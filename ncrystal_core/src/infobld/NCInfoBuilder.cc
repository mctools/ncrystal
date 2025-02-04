////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/infobld/NCInfoBuilder.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCIter.hh"
#include "NCrystal/internal/cfgutils/NCCfgManip.hh"

namespace NC = NCrystal;

struct NC::Info::InternalState {
  //NB: This struct only declared fully here, so only code in this .cc file will
  //be directly able to construct Info objects!
  shared_obj<const Info::Data> data;
  std::shared_ptr<const Info::OverrideableData> oData;
};

inline NC::Info::Info( InternalState&& state )
  : m_data(std::move(state.data)),
    m_oData(std::move(state.oData))
{
  nc_assert_always(!m_data->composition.empty());
}

NC::Info::InternalState NC::Info::copyInternalState() const
{
  return NC::Info::InternalState{ this->m_data, this->m_oData };
}

namespace NCRYSTAL_NAMESPACE {
  namespace InfoBuilder {
    namespace detail {

      std::size_t totalNumberOfAtomsInUnitCell( const AtomInfoList& atomlist )
      {
        //Derives number of atoms in a unit cell from atom list. Also
        //validatesthat atomlist nonempty, that all ai.numberPerUnitCell() > 0,
        //and returns the sum of those.
        std::size_t ntot(0);
        if ( atomlist.empty() )
          NCRYSTAL_THROW2(BadInput,"AtomInfoList must be non-empty if provided");
        for ( const auto& ai : atomlist ) {
          if ( ai.numberPerUnitCell() == 0 )
            NCRYSTAL_THROW(BadInput,"AtomInfo object should not have numberPerUnitCell()==0");
          nc_assert( ai.numberPerUnitCell() > 0 );
          ntot += ai.numberPerUnitCell();
        }
        nc_assert_always( ntot > 0 );
        return ntot;
      }

      template<class TList, class TGetFractionFct, class TSetFractionFct, class TGetDescrObjFct>
      void validateFractionListAndSnapToUnity( const char* fieldname, TList& alist,
                                               const TGetFractionFct& getfracfct,
                                               const TSetFractionFct& setfracfct,
                                               const TGetDescrObjFct& getdescrobjfct )
      {
        if ( alist.empty() )
          NCRYSTAL_THROW2(BadInput,"invalid "<<fieldname<<" : no entries!");

        StableSum ftot;
        for (const auto& e : alist ) {
          double fraction = getfracfct(e);
          if ( fraction<=0 || fraction > 1.0 ) {
            NCRYSTAL_THROW2(BadInput,"invalid "<<fieldname<<" fraction for "
                            <<getdescrobjfct(e)<<" : "<<fraction);
          }
          ftot.add( fraction );
        }
        const double ftotval = ftot.sum();
        if ( ftotval >= 1.000000001 || ftotval < 0.999999999 )
          NCRYSTAL_THROW2(BadInput,"invalid "<<fieldname<<" : fractions do not sum to unity");
        const double snapfact = 1.0 / ftotval;
        for (auto& e : alist )
          setfracfct(e,getfracfct(e)*snapfact);
      }

      namespace details {
        bool dhkl_compare( const NC::HKLInfo& rh, const NC::HKLInfo& lh )
        {
          nc_assert( !std::isnan(lh.dspacing) );
          nc_assert( !std::isnan(rh.dspacing) );
          nc_assert( !std::isinf(lh.dspacing) );
          nc_assert( !std::isinf(rh.dspacing) );
          if ( ncabs( lh.dspacing - rh.dspacing ) > 1e-6 )
            return lh.dspacing < rh.dspacing;
          double lh_fm = lh.fsquared*lh.multiplicity;
          double rh_fm = rh.fsquared*rh.multiplicity;
          nc_assert( !std::isnan(lh_fm) );
          nc_assert( !std::isnan(rh_fm) );
          nc_assert( !std::isinf(lh_fm) );
          nc_assert( !std::isinf(rh_fm) );
          if( ncabs( lh_fm - rh_fm ) > 1e-6 )
            return lh_fm < rh_fm;
          if (lh.multiplicity != rh.multiplicity)
            return lh.multiplicity < rh.multiplicity;
          return lh.hkl < rh.hkl;
        }
        bool atominfo_pos_compare( const NC::AtomInfo::Pos& rh, const NC::AtomInfo::Pos& lh )
        {
          if (rh[0]!=lh[0]) return rh[0] < lh[0];
          if (rh[1]!=lh[1]) return rh[1] < lh[1];
          return rh[2] < lh[2];
        }
        bool atominfo_pos_compare_yfirst( const NC::AtomInfo::Pos& rh, const NC::AtomInfo::Pos& lh )
        {
          if (rh[1]!=lh[1]) return rh[1] < lh[1];
          if (rh[0]!=lh[0]) return rh[0] < lh[0];
          return rh[2] < lh[2];
        }
        bool atominfo_pos_compare_zfirst( const NC::AtomInfo::Pos& rh, const NC::AtomInfo::Pos& lh )
        {
          if (rh[2]!=lh[2]) return rh[2] < lh[2];
          if (rh[1]!=lh[1]) return rh[1] < lh[1];
          return rh[0] < lh[0];
        }
        void atominfo_pos_remap( double& x ) {
          const double xorig = x;
          if ( x<0.0 )
            x += 1.0;
          else if ( x>=1.0 )
            x -= 1.0;
          if ( ! (x>=0 && x<1.0) )//must be in [0,1) and not NaN
            NCRYSTAL_THROW2(BadInput,"Invalid coordinate of atom position encountered (out of range or NaN): "<<xorig);
          if (x==0.0)
            x=0.0;//remaps -0 to 0
        }
        void detect_duplicate_positions(  const AtomInfo::AtomPositions& plist ) {
          const double pos_tolerance = 0.0001;//NB: Best if value matches the one in NCrystal/cifutils.py!
          for (std::size_t i = 1; i < plist.size(); ++i) {
            const AtomInfo::Pos& p1 = plist.at(i-1);
            const AtomInfo::Pos& p2 = plist.at(i);
            if ( ncabs(p1[0]-p2[0])<pos_tolerance && ncabs(p1[1]-p2[1])<pos_tolerance && ncabs(p1[2]-p2[2])<pos_tolerance )
              NCRYSTAL_THROW2(BadInput,"The same atom position used more than once: ("<<p1[0]<<", "<<p1[1]<<", "<<p1[2]<<")");
          }
        }
      }//details

      void validateAndCompleteUnitCell( InfoBuilder::UnitCell& uc )
      {
        auto& si = uc.structinfo;

        if ( uc.atomlist.has_value() ) {
          auto& uc_atomlist = uc.atomlist.value();
          //First validate total number of atoms is >0 and consistently specified:
          auto natompos = totalNumberOfAtomsInUnitCell( uc_atomlist );//throws if empty or have entries with 0 atoms
          nc_assert(natompos>0);
          if ( si.n_atoms != natompos )
            NCRYSTAL_THROW2(BadInput,"Inconsistent total number of atoms deduced from StructureInfo and AtomInfo ("
                            <<si.n_atoms<<" vs. "<<natompos<<")");

          //Sort for reproducibility:
          std::stable_sort(uc_atomlist.begin(),uc_atomlist.end(),
                           []( const NC::AtomInfo& rh, const NC::AtomInfo& lh )
                           { return rh.atomData() < lh.atomData(); } );

          //Check AtomInfo objects for consistency and standardise them a bit:
          AtomInfo::AtomPositions all_positions;
          all_positions.reserve(natompos);

          for ( auto& ai : uc_atomlist ) {
            //Map all atom positions to interval [0,1) (e.g. 1.0 becomes 0.0, -0.3
            //becomes 0.7, etc.). However at most move 1.0:
            for ( auto& pos : ai.detail_accessPos() ) {
              for ( auto& coord : pos )
                details::atominfo_pos_remap( coord );
              all_positions.push_back( pos );
            }

            //Sort positions, ensuring a well defined ordering:
            std::stable_sort( ai.detail_accessPos().begin(),
                              ai.detail_accessPos().end(),
                              details::atominfo_pos_compare );

            nc_assert( ai.numberPerUnitCell() > 0 );//was validated in totalNumberOfAtomsInUnitCell(..)

            //Debye temps and msd's should be specified for all or none of the AtomInfo objects:
            if ( ai.debyeTemp().has_value() != uc_atomlist.front().debyeTemp().has_value() )
              NCRYSTAL_THROW(BadInput,"Inconsistency: Debye temperatures specified for some but not all AtomInfo objects.");

            if ( ai.debyeTemp().has_value() )
              ai.debyeTemp().value().validate();

            if ( ai.msd().has_value() != uc_atomlist.front().msd().has_value() )
              NCRYSTAL_THROW(BadInput,"Inconsistency: MSD values specified for some but not all AtomInfo objects.");

          }

          //Ensure only one atom exists at a given position, within a tolerance. To make
          //sure this works we sort three times after x, y and z coordinates
          //respectively, and compare neighbouring elements each time:
          std::sort(all_positions.begin(),all_positions.end(),details::atominfo_pos_compare);
          details::detect_duplicate_positions(all_positions);
          std::sort(all_positions.begin(),all_positions.end(),details::atominfo_pos_compare_yfirst);
          details::detect_duplicate_positions(all_positions);
          std::sort(all_positions.begin(),all_positions.end(),details::atominfo_pos_compare_zfirst);
          details::detect_duplicate_positions(all_positions);
        }

        if ( ! ( si.n_atoms > 0 ) )
          NCRYSTAL_THROW2(BadInput,"StructureInfo n_atoms not a positive number: "<<si.n_atoms);

        checkAndCompleteLattice( si.spacegroup, si.lattice_a,
                                 si.lattice_b,  si.lattice_c );

        checkAndCompleteLatticeAngles( si.spacegroup, si.alpha, si.beta, si.gamma );


        if ( ! ( si.alpha > 0 && si.alpha < 180 &&
                 si.beta  > 0 && si.beta  < 180 &&
                 si.gamma > 0 && si.gamma < 180 ) ) {
          NCRYSTAL_THROW(BadInput,"Lattice angles must all be >0 and <180 degrees.");
        }

        if ( ! ( si.alpha > kPi || si.beta > kPi || si.gamma > kPi ) )
          NCRYSTAL_THROW(BadInput,"Lattice angles specified in radians where degrees were expected.");

        //Calculate unit cell volume:
        auto latRot = getLatticeRot( si.lattice_a, si.lattice_b, si.lattice_c,
                                     si.alpha*kDeg, si.beta*kDeg, si.gamma*kDeg );
        double expected_volume = latRot.colX().cross( latRot.colY() ).dot( latRot.colZ() );
        if ( si.volume > 0.0 && !floateq(si.volume,expected_volume,1e-3) )
          NCRYSTAL_THROW2(BadInput,"Provided ("<<si.volume<<"Aa3) versus calculated ("
                          <<expected_volume<<"Aa3) unit cell volume are incompatible!");
        si.volume = expected_volume;
        nc_assert( si.volume > 0.0 );
      }

      void validateAndCompleteDynamics( DynamicInfoList& dynamics )
      {
        if ( dynamics.empty() )
          NCRYSTAL_THROW2(BadInput,"If providing a DynamicInfoList, it must be non-empty");

        //Verify no duplicated entries in dynamics list and no links set up already:
        {
          std::set<IndexedAtomData> seen;
          for ( auto& di : dynamics ) {
            if ( !seen.insert(di->atom()).second )
              NCRYSTAL_THROW2(BadInput,"Multiple dynamic info sections for "<<di->atom());
            if ( di->correspondingAtomInfo() )
              NCRYSTAL_THROW2(BadInput,"Do not setup links on DynamicInfo objects before passing to InfoBuilder");
          }
        }
        //Validate+snap fractions:
        using DIPtr = DynamicInfoList::value_type;
        validateFractionListAndSnapToUnity( "DynamicInfoList", dynamics,
                                            [](const DIPtr& di){ nc_assert(di!=nullptr); return di->fraction(); },
                                            [](DIPtr& di, double f){ nc_assert(di!=nullptr); return di->changeFraction(f); },
                                            [](const DIPtr& di){ nc_assert(di!=nullptr); return di->atom(); } );

        //Sort for reproducibility:
        std::stable_sort(dynamics.begin(),dynamics.end(),
                         [](const std::unique_ptr<DynamicInfo>& a,const std::unique_ptr<DynamicInfo>& b)
                         {
                           nc_assert( !!a && !!b );
                           return ( a->atom() == b->atom()
                                    ? a->fraction() > b->fraction()
                                    : a->atom() < b->atom() );
                         });

        //Various:
        auto dblIsReallyInt = [](double x) { double dummy; return std::modf(x,&dummy) == 0.0; };

        for ( auto& di : dynamics ) {
          //Temperature and fractions validated elsewhere.

          //Nothing special for DI_Sterile or DI_FreeGas.
          auto di_sk = dynamic_cast<const DI_ScatKnl*>(di.get());
          if ( di_sk ) {
            auto egrid = di_sk->energyGrid();
            if ( egrid != nullptr ) {
              bool egrid_ok = egrid->size()>=3;
              if ( egrid_ok && egrid->size() == 3 ) {
                if (!dblIsReallyInt(egrid->at(2)))
                  egrid_ok = false;
                if ( egrid->at(1)!=0.0 && !( egrid->at(1) > egrid->at(0) ) )
                  egrid_ok = false;
              }
              if (!egrid_ok)
                NCRYSTAL_THROW2(BadInput,"Invalid energy grid returned from DI_ScatKnl.");
            }
            //Nothing further for DI_ScatKnlDirect.
            auto di_vdos = dynamic_cast<const DI_VDOS*>(di_sk);
            if ( di_vdos ) {
              if ( di_vdos->vdosOrigEgrid().empty() != di_vdos->vdosOrigDensity().empty() )
                NCRYSTAL_THROW2(BadInput,"DI_VDOS instance must provide none or both of vdosOrigEgrid and vdosOrigDensity");
              if ( ! floateq( di_vdos->atomData().averageMassAMU().dbl(),
                              di_vdos->vdosData().elementMassAMU().dbl() ) )
                NCRYSTAL_THROW2(BadInput,"DI_VDOS instance provides inconsistent mass values in atomData() and vdosData() fields.");
            }
            auto di_vdosdebye = dynamic_cast<const DI_VDOSDebye*>(di_sk);
            if ( di_vdosdebye )
              di_vdosdebye->debyeTemperature().validate();
          }//end di_ski
        }

      }

      void setupAtomInfoDynInfoLinks( AtomInfoList& atomlist, DynamicInfoList& dynamics )
      {
        nc_assert_always(!atomlist.empty());
        nc_assert_always(!dynamics.empty());
        if ( dynamics.size() != atomlist.size() )
          NCRYSTAL_THROW2(BadInput,"incompatible unit cell and dynamics info provided"
                          " (the two lists have a different number of atoms)");
        std::size_t nlinks(0);
        for ( auto& ai : atomlist ) {
          for ( auto& di : dynamics ) {
            if ( di->atom().index == ai.atom().index ) {
              ++nlinks;
              ai.detail_setupLink(di.get());//will throw if link has already been setup.
            }
          }
        }
        if ( nlinks != atomlist.size() )
          NCRYSTAL_THROW2(BadInput,"incompatible unit cell and dynamics info provided"
                          " (the two lists do not have the same IndexedAtomData fields present)");
      }

      void validateAndCompleteUnitCellAndDynamics( Optional<InfoBuilder::UnitCell>& unitcell,
                                                   Optional<DynamicInfoList>& dynamics )
      {
        if ( unitcell.has_value() )
          validateAndCompleteUnitCell( unitcell.value() );
        if ( dynamics.has_value() )
          validateAndCompleteDynamics( dynamics.value() );
        if ( unitcell.has_value() && unitcell.value().atomlist.has_value() && dynamics.has_value() ) {
          //Both dynamics and atomlist are present. First setup links:
          setupAtomInfoDynInfoLinks( unitcell.value().atomlist.value(), dynamics.value() );
          //Next, verify fractions are compatible and snap values to those
          //calculated from the unit cell (which can be calculated to full
          //machine precision). Also verify that Debye temperatures match up
          //where appropriate. This is easily done with the links once we know
          //the total number of atoms in the unit cell:
          const auto ntot = totalNumberOfAtomsInUnitCell( unitcell.value().atomlist.value() );
          nc_assert_always( ntot > 0 );
          for ( auto & di : dynamics.value() ) {
            auto ai = di->correspondingAtomInfo();
            nc_assert(ai!=nullptr);
            const double calc_frac = ai->numberPerUnitCell() / double(ntot);
            if ( ! floateq( calc_frac, di->fraction(), 1e-6, 1e-6 ) )
              NCRYSTAL_THROW2( BadInput,"Fractions specified in DynamicInfoList are incompatible with"
                               " those calculated from unit cell content for"<<di->atom()<<": "
                               <<  di->fraction()<<" vs " <<calc_frac );
            di->changeFraction(calc_frac);

            auto di_vdosdebye = dynamic_cast<const DI_VDOSDebye*>(di.get());
            if ( di_vdosdebye ) {
              if ( !ai->debyeTemp().has_value() )
                NCRYSTAL_THROW(BadInput,"AtomInfo object associated with DI_VDOSDebye object must have Debye temperature available!");
              if ( ai->debyeTemp().value() != di_vdosdebye->debyeTemperature() )
                NCRYSTAL_THROW(BadInput,"Associated AtomInfo and DI_VDOSDebye objects do not have the same Debye temperature specified!");
            }
          }
        }
      }

      void validateAtomIndexes( const Info::Composition& composition )
      {
        std::set<AtomIndex> seen;
        for ( auto& c : composition ) {
          if ( !seen.insert(c.atom.index).second )
            NCRYSTAL_THROW2(BadInput,"Invalid AtomIndex setup (repeated indices found in composition list)");
          if ( c.atom.index.get() >= composition.size() )
            NCRYSTAL_THROW2(BadInput,"Invalid AtomIndex setup (must be one of 0,...,ncomponents-1)");
        }
      }

      void validateAndCompleteComposition( Optional<Info::Composition>& composition,
                                           const Optional<UnitCell>& unitcell,
                                           Optional<DynamicInfoList>& dynamics )
      {
        const bool has_uc_atomlist = unitcell.has_value() && unitcell.value().atomlist.has_value();
        if ( composition.has_value() ) {
          //Composition is set, so just verify that dynamics/unitcell are not set, validate the provided fractions, and return.
          if ( has_uc_atomlist || dynamics.has_value() )
            NCRYSTAL_THROW(BadInput,"Do not set explicit composition on SinglePhaseBuilder when providing unitcell.atomlist or dynamics.");
          validateFractionListAndSnapToUnity( "composition list", composition.value(),
                                              [](const Info::CompositionEntry& e){ return e.fraction; },
                                              [](Info::CompositionEntry& e, double f){ e.fraction = f; },
                                              [](const Info::CompositionEntry& e){ return e.atom; } );
          return;
        }

        //Composition is absent, must be able to get from either unit cell or dynamics.
        if ( !has_uc_atomlist && !dynamics.has_value() ) {
          NCRYSTAL_THROW(BadInput,"SinglePhaseBuilder must have at least one of the following"
                         " pieces of information: composition, atomlist in unit cell, or dynamics.");
        }

        //Create composition based on either unit cell or dynamics (their
        //fractions have already been aligned but we anyway check the more
        //precise unitcell source first):

        Info::Composition new_composition;
        if ( has_uc_atomlist ) {
          //From unit cell:
          auto& atomlist = unitcell.value().atomlist.value();
          const auto ntot = totalNumberOfAtomsInUnitCell( atomlist );
          nc_assert_always( ntot > 0 );
          new_composition.reserve( atomlist.size() );
          for ( const auto& ai : atomlist ) {
            nc_assert( ai.numberPerUnitCell() > 0 );
            new_composition.emplace_back( ai.numberPerUnitCell() / double(ntot), ai.atom() );
          }
        } else {
          //From dynamics
          new_composition.reserve( dynamics.value().size() );
          for ( auto& di : dynamics.value() )
            new_composition.emplace_back( di->fraction(), di->atom() );
        }
        composition = std::move(new_composition);
      }

      void validateAndCompleteDSpacingRange( const PairDD& dspacingRange )
      {
        if ( ! (dspacingRange.second > dspacingRange.first )
             || !(dspacingRange.first>0.0) || !(dspacingRange.second>0.0) ) {
          NCRYSTAL_THROW2( BadInput, "Unvalid dspacingRange : [" << dspacingRange.first
                           << ", " << dspacingRange.second << "]." );
        }
      }

      void validateAndCompleteHKLList( HKLList& hkllist, const PairDD& dspacingRange )
      {
        validateAndCompleteDSpacingRange(dspacingRange);
        auto checkDSpacing = [&dspacingRange](double dsp)
        {
          if ( !(dsp>=dspacingRange.first) || !(dsp<=dspacingRange.second) )
            NCRYSTAL_THROW(BadInput,"Invalid HKL list produced. Some dspacing values are not in the requested range.");
        };

        std::stable_sort(hkllist.begin(),hkllist.end(),details::dhkl_compare);
        hkllist.shrink_to_fit();

        HKLInfoType hitype = ( hkllist.empty() ? HKLInfoType::Minimal : hkllist.front().type() );
        for ( auto& hkl : hkllist ) {
          if ( hkl.multiplicity < 2 || hkl.multiplicity > 99998 )
            NCRYSTAL_THROW(BadInput,"HKL multiplicity is not in range 2..99998");
          if ( hkl.multiplicity % 2 != 0 )
            NCRYSTAL_THROW(BadInput,"HKL multiplicity is not an even number.");
          if ( ! (hkl.fsquared >= 0.0 ) )
            NCRYSTAL_THROW(BadInput,"HKL fsquared is not a non-negative number");
          checkDSpacing(hkl.dspacing);
          if ( hkl.type() != hitype )
            NCRYSTAL_THROW(BadInput,"Inconsistency: HKLInfoType is not the same on all HKLInfo objects in the same list");
          if ( hitype == HKLInfoType::ExplicitHKLs ) {
            nc_assert(hkl.explicitValues != nullptr );
            auto& ll = hkl.explicitValues->list.get<std::vector<HKL>>();
            if ( ll.size()*2 != (std::size_t)hkl.multiplicity )
              NCRYSTAL_THROW(BadInput,"Explicit HKL values provided but number does not match multiplicity");
            ll.shrink_to_fit();
          } else if ( hitype == HKLInfoType::ExplicitNormals ) {
            auto& ll = hkl.explicitValues->list.get<std::vector<HKLInfo::Normal>>();
            if ( ll.size()*2 != (std::size_t)hkl.multiplicity )
              NCRYSTAL_THROW(BadInput,"Explicit HKL normals provided but number does not match multiplicity");
            ll.shrink_to_fit();
            for ( const auto& demi_normal : ll ) {
              if ( ! demi_normal.as<Vector>().isUnitVector() )
                NCRYSTAL_THROW(BadInput,"Provided demi_normals must have unit lengths");
            }
          }
        }
      }

      void validateTemperatures( const Optional<Temperature>& temperature,
                                 const Optional<DynamicInfoList>& dynamics )
      {
        Temperature tempval{-1.0};
        if ( temperature.has_value() ) {
          tempval = temperature.value();
          if ( ! ( tempval.get() >= Temperature::allowed_range.first
                   && tempval.get() <= Temperature::allowed_range.second ) )
            NCRYSTAL_THROW2(BadInput,"Invalid or out-of-range temperature value provided: "<<tempval);
        }
        if ( !dynamics.has_value() )
          return;
        if ( !temperature.has_value() )
          NCRYSTAL_THROW(BadInput,"Temperature is required whenever dynamic info is supplied");
        nc_assert(tempval.get()>0.0);
        tempval.validate();
        //
        //Check consistency of temperature fields (checking consistency of
        //DI_VDOSDebye Debye temperature fields elsewhere):
        for (auto& di : dynamics.value()) {

          if ( di->temperature() != tempval )
            NCRYSTAL_THROW(BadInput,"temperature info on DynamicInfo object"
                           " is different than temperature value of containing phase!");

          auto di_vdos = dynamic_cast<const DI_VDOS*>(di.get());
          if ( di_vdos && di_vdos->vdosData().temperature() != di_vdos->temperature() )
            NCRYSTAL_THROW(BadInput,"temperature info on VDOSData object provided by DI_VDOS object"
                           " is different than temperature on DI_VDOS object itself!");

          //NB: for DI_ScatKnlDirect we can't validate the temperature field on
          //the SABData object here, since it is built on-demand. However, we do
          //it upon construction in the ensureBuildThenReturnSAB() method.
        }
      }

      AtomMass calculateAverageAtomMass(const Info::Composition& c)
      {
        StableSum sum;
        for ( auto& di : c )
          sum.add( di.fraction * di.atom.data().averageMassAMU().dbl() );
        return AtomMass{ sum.sum() };
      }

      void validateDensities( Density d, NumberDensity nd )
      {
        if ( ! (d.dbl() >= 0.0 ) || std::isinf(d.dbl()) )
          NCRYSTAL_THROW2(BadInput,"Invalid density value: "<<d);
        if ( ! (nd.dbl() >= 0.0 ) || std::isinf(nd.dbl()) )
          NCRYSTAL_THROW2(BadInput,"Invalid number density value: "<<nd);
        if ( d.dbl() == 0.0 || nd.dbl() == 0.0 )
          NCRYSTAL_THROW(BadInput,"Materials with vanishing densities are not presently supported.");
        //General type validation:
        d.validate();
        nd.validate();
      }

      void validateAndCompleteDensities( AtomMass averageAtomMass,
                                         const Optional<UnitCell>& unitcell,
                                         Optional<Density>& density,
                                         Optional<NumberDensity>& numberDensity )
      {
        averageAtomMass.validate();
        nc_assert_always( averageAtomMass.get() > 0.0 );

        if ( density.has_value() && numberDensity.has_value() )
          NCRYSTAL_THROW(BadInput,"Do not supply both Density and NumberDensity on"
                         " SinglePhaseBuilder (supply at most one and the other will be calculated).");

        if ( unitcell.has_value() ) {
          //Always take final density/numberdensity values from unit cell when
          //possible - but check compatibility with any number already provided.
          const auto& si = unitcell.value().structinfo;
          nc_assert( si.volume > 0.0 && si.n_atoms > 0 );
          const auto numberDensityUC = NumberDensity{ si.n_atoms / si.volume };
          const auto densityUC = Density( numberDensityUC, averageAtomMass );
          if ( numberDensity.has_value() && !floateq( numberDensityUC.dbl(), numberDensity.value().dbl(), 1e-2 ) )
            NCRYSTAL_THROW2(BadInput,"Provided (" << numberDensity.value() << ") versus calculated-from-unit-cell"
                            " (" << numberDensityUC << ") number density values are incompatible!");
          if ( density.has_value() && !floateq( densityUC.dbl(), density.value().dbl(), 1e-2 ) )
            NCRYSTAL_THROW2(BadInput,"Provided (" << density.value() << ") versus calculated-from-unit-cell"
                            " (" << densityUC << ") density values are incompatible!");
          numberDensity = numberDensityUC;
          density = densityUC;
        } else {
          //If one of density/numberdensity is set, calculate the other:
          if ( density.has_value() && !numberDensity.has_value() )
            numberDensity = NumberDensity( density.value(), averageAtomMass );

          if ( numberDensity.has_value() && !density.has_value() )
            density = Density( numberDensity.value(), averageAtomMass );

          //Final requirement is always that a density must be available:
          if ( !numberDensity.has_value() || !density.has_value() )
            NCRYSTAL_THROW(BadInput,"Density/NumberDensity values must always be supplied directly or it"
                           " must be possible to deduce them (from each other or unit cell information).");
        }

        //Final sanity check:
        nc_assert_always( density.has_value() && numberDensity.has_value() );
        validateDensities( density.value(), numberDensity.value() );
      }

      void validateCustomData( const Info::CustomData& customData )
      {
        for (const auto& e : customData ) {
          if (e.first.empty() || !contains_only(e.first,"ABCDEFGHIJKLMNOPQRSTUVWXYZ"))
            NCRYSTAL_THROW2(BadInput,"invalid custom section name: \""<<e.first
                            <<"\" (must be non-empty and contain only capitalised letters A-Z)");
        }
      }

      void validateAndCompleteStateOfMatter( bool isCrystalline,
                                             const Optional<DynamicInfoList>& dynamics,
                                             Info::StateOfMatter& stateOfMatter )
      {
        //Verify and possibly fill out StateOfMatter:
        //
        //  Crystalline fields or presence of VDOS/VDOSDebye dynamics implies
        //  solid. In that case, verify that current state is either Unknown or
        //  Solid (and replace Unknown with Solid).
        //
        bool should_be_solid(isCrystalline);
        if ( dynamics.has_value() ) {
          for ( auto& di : dynamics.value() ) {
            if (should_be_solid)
              break;
            if (dynamic_cast<const DI_VDOSDebye*>(di.get())||dynamic_cast<const DI_VDOS*>(di.get()))
              should_be_solid = true;
          }
        }
        if ( should_be_solid ) {
          if ( !isOneOf( stateOfMatter, Info::StateOfMatter::Unknown, Info::StateOfMatter::Solid ) )
            NCRYSTAL_THROW2(BadInput,"Info objects that are crystalline or have at least one VDOS"
                            " (or VDOSDebye) can not be designated as \""<<Info::toString(stateOfMatter)<<"\"");
          stateOfMatter = Info::StateOfMatter::Solid;
        }
      }

      std::pair< std::vector<AtomDataSP>, VectS > createAtomDataSPAndLabelsLists( const Info::Composition& composition )
      {
        //TODO: We could combine the next two vectors into just:
        //std::vector<std::pair<IndexedAtomData,std::string>> That would be more
        //efficient, and we could even return const refs to IndexAtomData objects.
        std::pair< std::vector<AtomDataSP>, VectS > res;

        auto& displayLabels = res.second;
        displayLabels.resize(composition.size());

        //Order by atomdata, then index (so display labels won't depend on whether
        //or not the info factory did sensible sorting):
        std::vector<const IndexedAtomData*> v;
        v.reserve(composition.size());
        for ( const auto& e : composition )
          v.push_back( &e.atom );
        std::stable_sort(v.begin(),v.end(),
                         [](const IndexedAtomData* a,const IndexedAtomData* b)
                         {
                           if ( a->atomDataSP->getUniqueID() == b->atomDataSP->getUniqueID() )
                             return a->index < b->index;
                           return *a->atomDataSP < *b->atomDataSP;
                         });

        std::map<std::string,std::vector<AtomIndex>> lbl2indices;
        for ( const auto& e : v ) {
          const IndexedAtomData& iad = *e;
          const AtomData& ad = iad.atomDataSP;
          std::string lbl;
          if ( ad.isElement() ) {
            lbl = ad.elementName();
            if ( ad.isSingleIsotope() ) {
              unsigned A = ad.A();
              if (ad.Z()==1&&(A==2||A==3))
                lbl = ( A==2 ? "D" : "T" );
              else
                lbl += std::to_string(A);
            }
          } else {
            lbl = "Mix";
          }
          lbl2indices[lbl].push_back(iad.index);
        }
        auto idx_to_alphalbl = [](unsigned i) {
          static const std::string lc="abcdefghijklmnopqrstuvwxyz"_s;// a=0, b=1, c=2
          const unsigned nlc = static_cast<unsigned>(lc.size());
          std::string lbl;
          while (true) {
            lbl = lc.at(i%nlc)+lbl;
            if (i<nlc)
              break;
            i /= nlc;
            --i;
          }
          return lbl;
        };
        for (auto& e : lbl2indices ) {
          if (e.second.size()==1) {
            displayLabels.at(e.second.front().get()) = e.first;
          } else {
            for (auto&& ee: enumerate(e.second)) {
              displayLabels.at(ee.val.get()) = e.first + "("_s + idx_to_alphalbl(static_cast<unsigned>(ee.idx)) + ")"_s;
            }
          }
        }
        for (const auto& dl : displayLabels)
          nc_assert_always(!dl.empty());

        std::vector<OptionalAtomDataSP> tmp_atomDataSPs;
        tmp_atomDataSPs.resize(composition.size());

        for ( const auto& e : v ) {
          const IndexedAtomData& iad = *e;
          tmp_atomDataSPs.at(iad.index.get()) = iad.atomDataSP;
        }

        auto& atomDataSPs = res.first;
        atomDataSPs.reserve(tmp_atomDataSPs.size());
        for ( auto&& e : tmp_atomDataSPs )
          atomDataSPs.emplace_back(std::move(e));

        return res;
      }

      void validateDataSourceName( const DataSourceName& dsn )
      {
        if ( StrView(dsn.str().c_str(),dsn.str().size()).contains('\x0') )
          NCRYSTAL_THROW2(BadInput,"Null character encountered in data source name.");
      }

      void validateAndCompleteSinglePhaseInput( SinglePhaseBuilder& in )
      {
        validateDataSourceName(in.dataSourceName);

        //Unit cell, dynamics, and their connections:
        validateAndCompleteUnitCellAndDynamics( in.unitcell, in.dynamics );

        //Composition (taking from unitcell or dynamics if absent):
        validateAndCompleteComposition( in.composition, in.unitcell, in.dynamics );
        nc_assert(in.composition.has_value());
        validateAtomIndexes( in.composition.value() );

        //Check that temperature fields are valid and consistent:
        validateTemperatures( in.temperature, in.dynamics );

        if ( in.hklPlanes.has_value() ) {
          auto& hkl = in.hklPlanes.value();
          //Dspacing range:
          if ( !( hkl.dspacingRange.first < hkl.dspacingRange.second ) )
            NCRYSTAL_THROW2(BadInput,"Do not provide hklPlanes field with a dspacingRange of non-positive length");
          validateAndCompleteDSpacingRange(hkl.dspacingRange);

          //Validate and complete (sort etc.) HKL list if already available
          //(otherwise it will be delayed until construction):
          if ( hkl.source.has_value<HKLList>() )
            validateAndCompleteHKLList( hkl.source.get<HKLList>(), hkl.dspacingRange );
        }

        //Densities (based on average atomic mass from composition, calculating
        //from unit cell volume if needed and possible):
        nc_assert(in.composition.has_value());
        validateAndCompleteDensities( calculateAverageAtomMass( in.composition.value() ),
                                      in.unitcell, in.density, in.numberDensity );


        //Make sure that the "isCrystalline" concept is exactly as described in
        //NCInfo.hh, including which sub-fields must be available if others
        //are.: "Additionally, all crystalline phases will have HKL lists
        //available, and if Atom lists are available, Structure info must also
        //be available.". The latter requirements comes from the way we have set
        //up the InfoBuilder::UnitCell struct, so we just have to check that if
        //the unitcell field is set, hklPlanes must also be there.:

        if ( in.unitcell.has_value() && ! in.hklPlanes.has_value() )
          NCRYSTAL_THROW2(BadInput,"Info objects that have unit cell structure"
                          " available must always have hklPlanes available as well.");
        const bool isCrystalline = in.hklPlanes.has_value();

        //State of matter:
        validateAndCompleteStateOfMatter( isCrystalline, in.dynamics, in.stateOfMatter );

        //Custom data:
        if ( in.customData.has_value() )
          validateCustomData( in.customData.value() );
        //NB: We do not perform any validation of bkgdxsectprovider. We could in
        //principle query a few values and check they provide non-negative cross
        //sections.

      }

      void transferSinglePhaseData( SinglePhaseBuilder&& in, Info::Data& out )
      {
        //First derive some additional fields:
        nc_assert( in.composition.has_value() );
        std::tie(out.atomDataSPs,
                 out.displayLabels) = createAtomDataSPAndLabelsLists( in.composition.value() );

        //Non-moving transfers first for robustness:
        nc_assert_always(in.density.has_value());
        out.oData.fields.density.set( in.density.value() );
        nc_assert_always( in.numberDensity.has_value() );
        out.oData.fields.numberdensity.set( in.numberDensity.value() );
        out.temp.set( in.temperature );
        out.stateOfMatter = in.stateOfMatter;

        out.detail_braggthreshold = -1.0;
        out.detail_hklInfoType = Info::Data::hKLInfoTypeInt_unsetval;

        if ( in.hklPlanes.has_value() ) {
          auto& hkl = in.hklPlanes.value();
          out.hkl_dlower_and_dupper = hkl.dspacingRange;
          if ( hkl.source.has_value<HKLList>() ) {
            //List is provided directly:
            out.detail_hkllist_needs_init = false;
            out.detail_hklList = std::move( hkl.source.get<HKLList>() );
            out.detail_braggthreshold = ( out.detail_hklList.empty() ? 0.0 : out.detail_hklList.front().dspacing * 2.0 );
            out.detail_hklInfoType = enumAsInt( out.detail_hklList.empty()  ? HKLInfoType::Minimal : out.detail_hklList.front().type() );
          } else {
            //Delayed init via generator function. We wrap the genfct call to
            //take care of StructureInfo/AtomInfoList arguments + make sure we
            //always ends up calling validateAndCompleteHKLList.
            using HKLListGenFct = SinglePhaseBuilder::HKLPlanes::HKLListGenFct;
            nc_assert_always( hkl.source.has_value<HKLListGenFct>() );
            HKLListGenFct genfct = hkl.source.get<HKLListGenFct>();
            const Info::Data * dataptr = &out;//Ok since Info::Data is NoCopyMove
            std::function<HKLList(PairDD)> wrapped_genfct = [genfct,dataptr](PairDD dspacingRange) -> HKLList
            {
              nc_assert(dataptr!=nullptr);
              nc_assert(genfct!=nullptr);
              auto siptr = ( dataptr->structinfo.has_value() ? &dataptr->structinfo.value() : nullptr );
              auto aiptr = ( dataptr->atomlist.empty() ? nullptr : &dataptr->atomlist );
              auto hkllist = genfct( siptr, aiptr, dspacingRange );
              validateAndCompleteHKLList( hkllist, dspacingRange );
              return hkllist;
            };
            out.detail_hkllist_needs_init = true;
            out.hkl_ondemand_fct = wrapped_genfct;
          }
        }

        //All that is left now are simple moves:
        if ( in.dataSourceName.str().empty() ) {
          static const DataSourceName s_dsn_unknown = DataSourceName("<unknown>");
          out.dataSourceName = s_dsn_unknown;
        } else {
          out.dataSourceName = std::move(in.dataSourceName);
        }

        if ( in.unitcell.has_value() ) {
          out.structinfo = std::move( in.unitcell.value().structinfo );
          if ( in.unitcell.value().atomlist.has_value() )
            out.atomlist = std::move( in.unitcell.value().atomlist.value() );
        }
        if ( in.dynamics.has_value() )
          out.dyninfolist = std::move( in.dynamics.value() );

        out.composition = std::move( in.composition.value() );
        if ( in.customData.has_value() )
          out.custom = std::move( in.customData.value() );

        out.xsectprovider = std::move(in.bkgdxsectprovider);
      }

      void validateAndCompleteMultiPhaseInput( MultiPhaseBuilder& in )
      {
        //Ensure non-empty list and check/snap fractions:
        validateFractionListAndSnapToUnity( "phase list", in.phases,
                                            [](const Info::Phase& p){ return p.first; },
                                            [](Info::Phase& p, double f){ p.first = f;; },
                                            [](const Info::Phase&){ return "phase"; } );
      }

      void transferMultiPhaseData( MultiPhaseBuilder&& in, Info::Data& out )
      {
        nc_assert_always( in.phases.size() >= 2 );
        out.oData.fields.phases = makeSO<Info::PhaseList>( std::move(in.phases) );
        const auto& phases = *out.oData.fields.phases;

        static const DataSourceName s_dsn_mp = DataSourceName("<multiphase>");
        out.dataSourceName = s_dsn_mp;

        ///////////////////////////////////////////////
        // Calculate parameters derived from phases: //
        ///////////////////////////////////////////////

        //Set temperature only if present and identical on all phases:
        {
          Optional<Temperature> common_temperature;
          if ( phases.back().second->hasTemperature() ) {
            common_temperature = phases.back().second->getTemperature();
            for ( auto i : ncrange(phases.size()-1) ) {
              if ( ! phases.at(i).second->hasTemperature() ||
                   phases.at(i).second->getTemperature() != common_temperature.value() ) {
                common_temperature.reset();
                break;
              }
            }
          }
          out.temp.set( common_temperature );
        }

        //Set state of matter only if known and identical on all phases:
        {
          Info::StateOfMatter stateOfMatter = phases.back().second->stateOfMatter();
          if ( stateOfMatter != Info::StateOfMatter::Unknown ) {
            for ( auto i : ncrange(phases.size()-1) ) {
              if ( stateOfMatter != phases.at(i).second->stateOfMatter() ) {
                stateOfMatter = Info::StateOfMatter::Unknown;
                break;
              }
            }
          }
          out.stateOfMatter = stateOfMatter;
        }

        //Density/number density:
        {
          StableSum densitySum;
          StableSum numberDensitySum;
          for ( const auto& ph : phases ) {
            densitySum.add( ph.first * ph.second->getDensity().dbl() );
            numberDensitySum.add( ph.first * ph.second->getNumberDensity().dbl() );
          }
          out.oData.fields.density = Density{ densitySum.sum() };
          out.oData.fields.numberdensity = NumberDensity{ numberDensitySum.sum() };
        }

        //Composition:
        {
          //Go through phases and combine the AtomData and fraction information in
          //their composition. In case the same AtomData object is used in multiple
          //phases, we combine them here.
          //
          //Important note: The phase fractions indicate the volume-percentage
          //of each phase. For the composition we on the other hand need to
          //weight with that phase's relative contribution to the total number
          //density, which is the phase volume fraction multiplied by its number
          //density.

          nc_assert( out.oData.fields.numberdensity.dbl() > 0.0 );
          const double combined_numberdensity = out.oData.fields.numberdensity.dbl();

          std::map<AtomDataSP,StableSum> atomdata2frac;//collect fractions and associated AtomDataSP objects
          for ( const auto& ph : phases ) {
            const double phase_volume_fraction = ph.first;
            nc_assert(phase_volume_fraction>0.0&&phase_volume_fraction<=1.0);
            const double phase_numberdensity = ph.second->getNumberDensity().dbl();
            const double phase_atomcount_fraction = phase_volume_fraction * ( phase_numberdensity / combined_numberdensity );
            for ( auto & cc : ph.second->getComposition() ) {
              nc_assert(cc.fraction>0.0&&cc.fraction<=1.0);
              const double f = phase_atomcount_fraction * cc.fraction;
              nc_assert(f>=0.0&&f<=1.0);
              if ( f > 0.0 )
                atomdata2frac[cc.atom.atomDataSP].add(f);
            }
          }
          //Setup new composition list (sorting will be done later);
          Info::Composition comp;
          comp.reserve(atomdata2frac.size());
          for ( auto& e : atomdata2frac ) {
            comp.emplace_back( e.second.sum(),
                               IndexedAtomData{e.first,AtomIndex::createInvalidObject()});
          }
          out.composition = std::move(comp);
        }

      }

      void finalCommonValidateAndComplete( Info::Data& data )
      {
        //Sort composition for reproducibility:
        std::stable_sort( data.composition.begin(), data.composition.end(),
                          [](const Info::CompositionEntry& a,const Info::CompositionEntry& b)
                          {
                            return ( a.atom == b.atom
                                     ? a.fraction > b.fraction
                                     : a.atom < b.atom );
                          });
        //sanity check temperature:
        if (data.temp.has_value())
          data.temp.value().validate();
      }
    }
  }
}

NC::InfoPtr NC::InfoBuilder::buildInfoPtr( InfoPtr orig, Density density )
{
  detail::validateDensities( density, NumberDensity{1.0} );
  if ( orig->getDensity() == density )
    return orig;
  nc_assert_always( orig->getDensity().dbl()>0.0 );
  return buildInfoPtrWithScaledDensity( orig, density.dbl() / orig->getDensity().dbl() );
}

NC::InfoPtr NC::InfoBuilder::buildInfoPtr( InfoPtr orig, NumberDensity numberDensity )
{
  detail::validateDensities( Density{1.0}, numberDensity );
  if ( orig->getNumberDensity() == numberDensity )
    return orig;
  nc_assert_always( orig->getNumberDensity().dbl()>0.0 );
  return buildInfoPtrWithScaledDensity( orig, numberDensity.dbl() / orig->getNumberDensity().dbl() );
}

NC::Info NC::InfoBuilder::buildInfo( SinglePhaseBuilder&& input )
{
  InfoBuilder::detail::validateAndCompleteSinglePhaseInput(input);
  auto dataptr = makeSO<Info::Data>();
  Info::Data& data = *dataptr;
  InfoBuilder::detail::transferSinglePhaseData(std::move(input),data);
  InfoBuilder::detail::finalCommonValidateAndComplete(data);
  return Info( Info::InternalState{ std::move(dataptr), nullptr });
}

NC::Info NC::InfoBuilder::buildInfo( MultiPhaseBuilder&& input )
{
  //NB: Apart from recognising that all phases are the same, we do not in
  //general merge identical phases in the input.phases list. This is because it
  //could potentially screw up with SANS physics.

  InfoBuilder::detail::validateAndCompleteMultiPhaseInput(input);

  {
    //Check if all phases are the same, if so we do not actually have a
    //multi-phase material and we simply return that phase instead:
    nc_assert(!input.phases.empty());
    auto ph0 = input.phases.back();
    bool all_phases_same = true;
    for ( auto i : ncrange(input.phases.size()-1) ) {
      if ( input.phases.at(i).second != ph0.second ) {
        all_phases_same = false;
        break;
      }
    }
    if (all_phases_same)
      return Info( ph0.second->copyInternalState() );//aka "clone"
  }

  //Looks like true multiphase material:
  auto dataptr = makeSO<Info::Data>();
  Info::Data& data = *dataptr;
  InfoBuilder::detail::transferMultiPhaseData(std::move(input),data);
  InfoBuilder::detail::finalCommonValidateAndComplete(data);

  //One last thing to do is to see if any CfgData items are shared between all
  //the child phases. If so, we should register those on the mother object:
  class CfgDataIter {
    Info::PhaseList::const_iterator m_it,m_itE;
  public:
    CfgDataIter(const Info::PhaseList& pl) : m_it(pl.begin()), m_itE(pl.end()) {}
    const Cfg::CfgData* operator()() { return m_it == m_itE ? nullptr : &( (m_it++)->second->getCfgData() ); }
  };

  nc_assert(data.oData.fields.phases!=nullptr);
  auto common_entries = Cfg::CfgManip::findCommonEntries( CfgDataIter(*data.oData.fields.phases) );
  if (!common_entries.empty()) {
    auto varfilter_onlycommon = Cfg::CfgManip::createFilter( common_entries, Cfg::CfgManip::FilterType::OnlyListed );
    Cfg::CfgManip::apply( data.oData.fields.cfgData,
                          data.oData.fields.phases->front().second->getCfgData(),
                          varfilter_onlycommon );
  }

  return Info( Info::InternalState{ std::move(dataptr), nullptr });
}

NC::InfoPtr NC::InfoBuilder::buildInfoPtr( SinglePhaseBuilder&& bldr ) { return makeSO<const Info>( buildInfo( std::move(bldr) ) ); }
NC::InfoPtr NC::InfoBuilder::buildInfoPtr( MultiPhaseBuilder&& bldr )  { return makeSO<const Info>( buildInfo( std::move(bldr) ) ); }

#include "NCrystal/internal/utils/NCAtomUtils.hh"
#include "NCrystal/internal/atomdb/NCAtomDB.hh"

NC::Info::Composition NC::InfoBuilder::buildCompositionFromChemForm( const std::string& cf_str )
{
  auto cf = tryDecodeSimpleChemicalFormula( cf_str );
  if ( !cf.has_value() )
    NCRYSTAL_THROW2(BadInput,"Could not decode chemical formula (needed for composition): \""<<cf_str<<"\"");
  Info::Composition composition;
  uint64_t nelem_tot(0);
  for ( auto& n_smb : cf.value() )
    nelem_tot += n_smb.first;
  nc_assert_always(nelem_tot>0);
  for ( auto& n_smb : cf.value() ) {
    nc_assert_always(n_smb.second.isElement());
    auto atom = AtomDB::getNaturalElement( n_smb.second.Z() );
    if ( !atom )
      NCRYSTAL_THROW2(DataLoadError,"Does not have data for element with Z="<<n_smb.second.Z());
    composition.emplace_back( n_smb.first / double(nelem_tot),
                              IndexedAtomData{atom,AtomIndex{static_cast<AtomIndex::value_type>(composition.size())}} );
  }
  return composition;
}

namespace NCRYSTAL_NAMESPACE {
  namespace InfoBuilder {
    namespace {

      bool oDataFieldsEqual( const Info::OverrideableDataFields& a,
                             const Info::OverrideableDataFields& b )
      {
        nc_assert(!std::isnan(a.density.dbl()));
        nc_assert(!std::isnan(b.density.dbl()));
        nc_assert(!std::isnan(a.numberdensity.dbl()));
        nc_assert(!std::isnan(b.numberdensity.dbl()));
        if ( a.density != b.density )
          return false;
        if ( a.numberdensity != b.numberdensity )
          return false;
        if (!Cfg::CfgManip::equal(a.cfgData,b.cfgData))
          return false;
        if ( (a.phases == nullptr) != (b.phases == nullptr) )
          return false;
        if ( a.phases != nullptr ) {
          auto& a_pl = *a.phases;
          auto& b_pl = *b.phases;
          if ( a_pl.size() != b_pl.size() )
            return false;
          for ( auto i : ncrange(a_pl.size()) ) {
            auto& pa = a_pl.at(i);
            auto& pb = b_pl.at(i);
            nc_assert(!std::isnan(pa.first));
            nc_assert(!std::isnan(pb.first));
            if ( pa.first != pb.first )
              return false;
            if ( pa.second->getUniqueID() != pb.second->getUniqueID() )
              return false;
          }
        }
        return true;
      }

      InfoPtr overrideInfoFieldsWithCache( InfoPtr orig,
                                           const Info::OverrideableDataFields& fields )
      {
        if ( oDataFieldsEqual( orig->ovrData(), fields ) )
          return orig;
        auto state =  orig->copyInternalState();
        if ( oDataFieldsEqual( state.data->oData.fields, fields ) ) {
          //compatible with underlying Info object:
          state.oData.reset();
          return makeSO<const Info>( std::move(state) );
        }
        auto& cache = state.data->oDataCache;
        NCRYSTAL_LOCK_GUARD(cache.mtx);
        auto& v = cache.entries;
        for ( auto& e : v ) {
          if ( oDataFieldsEqual( e->fields, fields ) ) {
            //already in cache
            state.oData = e;
            return makeSO<const Info>( std::move(state) );
          }
        }
        auto newOData = makeSO<Info::OverrideableData>();
        newOData->fields = fields;
        state.oData = std::move(newOData);
        static constexpr std::size_t nlim = 1000;
        static constexpr std::size_t nlimpreservefrontback = 250;
        static_assert(nlimpreservefrontback*2<=nlim,"");
        if ( v.size() == nlim ) {
          //TODO: unit test this (was tested manually by lowering numbers above from 1000/250 to 3/1)
          //Protect against out-of-control caches while retaining some
          //repeatability. Preserve oldest + youngest and hope middle ones are most
          //unlikely to be requested again:
          std::decay<decltype(v)>::type v2;
          v2.reserve(nlimpreservefrontback*2+1);//+1 since we will do a final push_back below
          for ( auto i : ncrange(nlimpreservefrontback) )
            v2.push_back( std::move(v.at(i)) );
          for ( auto i : ncrange(nlimpreservefrontback) )
            v2.push_back( std::move(v.at((nlim-nlimpreservefrontback)+i)) );
          v = std::move(v2);
        }
        v.push_back( state.oData );
        return makeSO<const Info>( std::move(state) );
      }
    }
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace InfoBuilder {
    namespace {
      bool detail_phaseListsIdentical( const Info::PhaseList& a, const Info::PhaseList& b )
      {
        if ( a.size() != b.size() )
          return false;
        for ( auto i : ncrange(a.size()) ) {
          auto pa = a.at(i);
          auto pb = b.at(i);
          if ( ! (pa.first == pb.first ) )
            return false;
          if ( pa.second->getUniqueID() != pb.second->getUniqueID() )
            return false;
        }
        return true;
      }
    }
  }
}

NC::InfoPtr NC::InfoBuilder::recordCfgDataOnInfoObject( InfoPtr orig, const Cfg::CfgData& data )
{
  if ( Cfg::CfgManip::empty(data) )
    return orig;

  Cfg::CfgData newdata = orig->getCfgData();
  Cfg::CfgManip::apply( newdata, data,
                                  [](Cfg::VarId varid){ return Cfg::varGroup(varid) != Cfg::VarGroupId::Info; } );
  if ( Cfg::CfgManip::equal( newdata, orig->getCfgData() ) )
    return orig;//no changes (even if multiphase, since orig->getCfgData() contains the items common to all phases)

  Info::OverrideableDataFields fields;
  fields.density = orig->getDensity();
  fields.numberdensity = orig->getNumberDensity();
  fields.cfgData = std::move(newdata);
  if ( orig->isMultiPhase() ) {
    //Multiphase, we must record on all children as well. No way around breaking
    //up, recording, and reassembling the phase list.
    auto new_phases = makeSO<Info::PhaseList>();
    auto& phases = orig->getPhases();
    new_phases->reserve(phases.size());
    for ( auto ph : phases )
      new_phases->emplace_back( ph.first,
                                recordCfgDataOnInfoObject( ph.second, data ) );
    if ( detail_phaseListsIdentical( phases, new_phases ) )
      fields.phases = orig->detail_getPhasesSP();//prefer orig shared instance
    else
      fields.phases = std::move(new_phases);
  }
  return overrideInfoFieldsWithCache( orig, fields );
}

NC::InfoPtr NC::InfoBuilder::buildInfoPtrWithScaledDensity( InfoPtr orig, double scaleFactor )
{
  nc_assert_always(scaleFactor>=0.0);
  if ( scaleFactor == 1.0 )
    return orig;
  Info::OverrideableDataFields fields;
  fields.numberdensity = NumberDensity{ orig->getNumberDensity().dbl() * scaleFactor };
  fields.density = Density{ orig->getDensity().dbl() * scaleFactor };
  fields.cfgData = orig->getCfgData();
  if ( orig->isMultiPhase() ) {
    //Multiphase, we must scale all children as well. No way around breaking up,
    //scaling, and reassembling the phase list.
    auto new_phases = makeSO<Info::PhaseList>();
    auto& phases = orig->getPhases();
    new_phases->reserve(phases.size());
    for ( auto ph : phases )
      new_phases->emplace_back( ph.first,
                                buildInfoPtrWithScaledDensity( ph.second, scaleFactor ) );
    if ( detail_phaseListsIdentical( phases, new_phases ) )
      fields.phases = orig->detail_getPhasesSP();//prefer orig shared instance
    else
      fields.phases = std::move(new_phases);
  }
  return overrideInfoFieldsWithCache( orig, fields );
}
