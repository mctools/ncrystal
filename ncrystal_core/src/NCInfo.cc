////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
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

#include "NCrystal/NCInfo.hh"
#include "NCrystal/internal/NCLatticeUtils.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCIter.hh"
#include <algorithm>
#include <cstring>//for memcpy
#include <cstdlib>
namespace NC=NCrystal;

NC::Info::Info()
  : m_hkl_dlower(-1.0),
    m_hkl_dupper(-2.0),
    m_density(-1.0),
    m_numberdensity(-1.0),
    m_xsect_free(-1.0),
    m_xsect_absorption(-1.0),
    m_temp(-1.0),
    m_debyetemp(-1.0),
    m_lock(false)
{
  m_structinfo.spacegroup = 999999;//unset
}

NC::Info::~Info() = default;

namespace NCrystal {

  bool dhkl_compare( const NC::HKLInfo& rh, const NC::HKLInfo& lh )
  {
    if( ncabs(lh.dspacing-rh.dspacing) > 1.0e-6 )
      return lh.dspacing < rh.dspacing;
    else if( ncabs(lh.fsquared*lh.multiplicity - rh.fsquared*rh.multiplicity) > 1.0e-6 )
      return lh.fsquared*lh.multiplicity < rh.fsquared*rh.multiplicity;
    else if (lh.multiplicity != rh.multiplicity)
      return lh.multiplicity < rh.multiplicity;
    else if( lh.h!=rh.h )
      return rh.h < lh.h;
    else if( lh.k!=rh.k )
      return rh.k < lh.k;
    return rh.l < lh.l;
  }
  bool atominfo_compare( const NC::AtomInfo& rh, const NC::AtomInfo& lh )
  {
    return rh.atom < lh.atom;
  }
  bool atominfo_pos_compare( const NC::AtomInfo::Pos& rh, const NC::AtomInfo::Pos& lh )
  {
    if (rh.x!=lh.x) return rh.x < lh.x;
    if (rh.y!=lh.y) return rh.y < lh.y;
    return rh.z < lh.z;
  }
  bool atominfo_pos_compare_yfirst( const NC::AtomInfo::Pos& rh, const NC::AtomInfo::Pos& lh )
  {
    if (rh.y!=lh.y) return rh.y < lh.y;
    if (rh.x!=lh.x) return rh.x < lh.x;
    return rh.z < lh.z;
  }
  bool atominfo_pos_compare_zfirst( const NC::AtomInfo::Pos& rh, const NC::AtomInfo::Pos& lh )
  {
    if (rh.z!=lh.z) return rh.z < lh.z;
    if (rh.y!=lh.y) return rh.y < lh.y;
    return rh.x < lh.x;
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
  void detect_duplicate_positions(  std::vector<AtomInfo::Pos>& plist ) {
    const double pos_tolerance = 0.01;
    for (std::size_t i = 1; i < plist.size(); ++i) {
      const AtomInfo::Pos& p1 = plist.at(i-1);
      const AtomInfo::Pos& p2 = plist.at(i);
      if ( ncabs(p1.x-p2.x)<pos_tolerance && ncabs(p1.y-p2.y)<pos_tolerance && ncabs(p1.z-p2.z)<pos_tolerance )
        NCRYSTAL_THROW2(BadInput,"The same atom position used more than once: ("<<p1.x<<", "<<p1.y<<", "<<p1.z<<")");
    }
  }
}

void NC::Info::objectDone()
{
  //TODO: Throw LogicErrors or BadInput here?
  ensureNoLock();
  m_lock=true;

  //avoid excess memory usage in hkl list:
  nc_assert_always(m_hkllist.empty()||hasHKLInfo());
  m_hkllist.shrink_to_fit();

  //Check that no nullptr AtomDataSP were provided:
  for (const auto& e: m_composition) {
    if (!e.atom.atomDataSP)
      NCRYSTAL_THROW2(BadInput,"AtomData in provided composition is a nullptr!");
  }
  for (const auto& e : m_atomlist) {
    if (!e.atom.atomDataSP)
      NCRYSTAL_THROW2(BadInput,"AtomData in provided AtomInfo is a nullptr!");
  }

  //sort lists for reproducibility:
  std::stable_sort(m_hkllist.begin(),m_hkllist.end(),dhkl_compare);
  std::stable_sort(m_atomlist.begin(),m_atomlist.end(),atominfo_compare);
  std::stable_sort(m_composition.begin(),m_composition.end(),
                   [](const CompositionEntry& a,const CompositionEntry& b)
                   {
                     nc_assert(!!a.atom.atomDataSP && !!b.atom.atomDataSP);
                     return ( a.atom==b.atom
                              ? a.fraction > b.fraction
                              : a.atom < b.atom );
                   });
  std::stable_sort(m_dyninfolist.begin(),m_dyninfolist.end(),
                   [](const std::unique_ptr<const DynamicInfo>& a,const std::unique_ptr<const DynamicInfo>& b)
                   {
                     nc_assert( !!a && !!b );
                     return ( a->atom() == b->atom()
                              ? a->fraction() > b->fraction()
                              : a->atom() < b->atom() );
                   });

  //check that per-element debye temp, positions and MSD's are consistently
  //specified (either all or none must have):
  std::set<unsigned> z_seen;
  std::vector<AtomInfo::Pos> all_positions;
  all_positions.reserve(128);

  unsigned ntotatoms_from_atominfo(0);

  std::map<IndexedAtomData,double> composition_map_fromatompos;
  AtomList::iterator itAtm(m_atomlist.begin()), itAtmE(m_atomlist.end());
  for (;itAtm!=itAtmE;++itAtm)
  {
    nc_assert(itAtm->atom.atomDataSP);

    //Map all atom positions to interval [0,1) (e.g. 1.0 becomes 0.0, -0.3
    //becomes 0.7, etc.). However at most move 1.0:
    std::vector<AtomInfo::Pos>::iterator itPos(itAtm->positions.begin()), itPosE(itAtm->positions.end());
    for (;itPos!=itPosE;++itPos) {
      atominfo_pos_remap(itPos->x);
      atominfo_pos_remap(itPos->y);
      atominfo_pos_remap(itPos->z);
      all_positions.push_back(*itPos);
    }
    //Sort positions, ensuring a well defined ordering:
    std::stable_sort(itAtm->positions.begin(),itAtm->positions.end(),atominfo_pos_compare);

    if (!itAtm->number_per_unit_cell)
      NCRYSTAL_THROW(BadInput,"atominfo found for element with number_per_unit_cell=0!.");

    ntotatoms_from_atominfo += itAtm->number_per_unit_cell;

    //per-element debye temps should be 0 (indicating absence) or have a sane value between 1K and 100000K.
    if ( ! (itAtm->debye_temp==0.0 || (itAtm->debye_temp>1.0&&itAtm->debye_temp<100000.0) ) )
      NCRYSTAL_THROW(LogicError,"Inconsistency: per-element Debye temperatures with invalid values encountered.");
    if ( bool(itAtm->debye_temp>0) != bool(m_atomlist.front().debye_temp>0) )
      NCRYSTAL_THROW(LogicError,"Inconsistency: per-element Debye temperatures specified for some but not all elements.");
    if ( itAtm->positions.empty() != m_atomlist.front().positions.empty() )
      NCRYSTAL_THROW(LogicError,"Inconsistency: positions specified for some but not all elements.");
    if ( !itAtm->positions.empty() && itAtm->positions.size()!=itAtm->number_per_unit_cell)
      NCRYSTAL_THROW(LogicError,"Inconsistency: inconsistency between length of positions vector and number_per_unit_cell");
    if ( ! ( itAtm->mean_square_displacement >= 0.0 ) )//will catch NaN as well
      NCRYSTAL_THROW(LogicError,"Inconsistency: mean_square_displacement must be >= 0.0 and not NaN");
    if ( bool(itAtm->mean_square_displacement>0) != bool(m_atomlist.front().mean_square_displacement>0) )
      NCRYSTAL_THROW(LogicError,"Inconsistency: mean_square_displacement specified for some but not all elements.");

    composition_map_fromatompos[itAtm->atom] += itAtm->number_per_unit_cell;
  }

  for (auto& e: composition_map_fromatompos)
    e.second /= ntotatoms_from_atominfo;

  //Ensure only one atom exists at a given position, within a tolerance. To make
  //sure this works we sort three times after x, y and z coordinates
  //respectively, and compare neighbouring elements each time:
  std::sort(all_positions.begin(),all_positions.end(),atominfo_pos_compare);
  detect_duplicate_positions(all_positions);
  std::sort(all_positions.begin(),all_positions.end(),atominfo_pos_compare_yfirst);
  detect_duplicate_positions(all_positions);
  std::sort(all_positions.begin(),all_positions.end(),atominfo_pos_compare_zfirst);
  detect_duplicate_positions(all_positions);

  //Check that hkl normal information is self-consistent:
  nc_assert_always(sizeof(HKLInfo::Normal)==3*sizeof(double));

  HKLList::iterator ithkl = m_hkllist.begin();
  HKLList::iterator ithklE = m_hkllist.end();
  int has_demi_normals(-1);
  int has_eqv_hkl(-1);
  for (;ithkl!=ithklE;++ithkl) {
    if (! ithkl->demi_normals.empty() ) {
      if ( ithkl->demi_normals.size() != ithkl->demi_normals.capacity() ) {
        //Remove over-capacity:
        std::vector<HKLInfo::Normal>(ithkl->demi_normals.begin(),ithkl->demi_normals.end()).swap(ithkl->demi_normals);
      }
      if (has_demi_normals==0)
        NCRYSTAL_THROW(LogicError,"Inconsistency: Some but not all HKLInfo objects provide demi_normals");
      has_demi_normals=1;
      if (ithkl->multiplicity < 1 || ithkl->multiplicity > 99999)
        NCRYSTAL_THROW(LogicError,"HKL multiplicity is not in range 1..99999");
      if (ithkl->demi_normals.size()*2 != (size_t)ithkl->multiplicity)
        NCRYSTAL_THROW(LogicError,"HKL normals provided but number does not match multiplicity");

      if (has_eqv_hkl!=-1 && (ithkl->eqv_hkl?1:0)!=has_eqv_hkl )
        NCRYSTAL_THROW(LogicError,"Inconsistency: Some but not all HKLInfo objects provide eqv_hkl");
      has_eqv_hkl = (ithkl->eqv_hkl?1:0);
      //check demi-normals are normalised:
      std::vector<HKLInfo::Normal>::const_iterator itN, itNE = ithkl->demi_normals.end();
      for (itN = ithkl->demi_normals.begin();itN!=itNE;++itN) {
        double n2 = itN->x*itN->x + itN->y*itN->y + itN->z*itN->z;
        if (ncabs(n2-1.0)>1.0e-6)
          NCRYSTAL_THROW(BadInput,"Provided demi_normals must have unit lengths");
      }
    } else {
      if (has_demi_normals==1)
        NCRYSTAL_THROW(LogicError,"Inconsistency: Some but not all HKLInfo objects provide demi_normals");
      has_demi_normals=0;
      if (ithkl->eqv_hkl)
        NCRYSTAL_THROW(LogicError,"eqv_hkl provided although demi_normals are not!");
    }
    if ( ithkl->eqv_hkl && ithkl->multiplicity%2 != 0 )
      NCRYSTAL_THROW(LogicError,"Expanded HKL info (eqv_hkl) provided, but multiplicity is not an even number.");
  }

  if (hasStructureInfo()) {
    if ( ! (m_structinfo.volume > 0.0) )
      NCRYSTAL_THROW2(BadInput,"StructureInfo volume not a positive number: "<<m_structinfo.volume);
    if ( ! (m_structinfo.n_atoms > 0.0) )
      NCRYSTAL_THROW2(BadInput,"StructureInfo n_atoms not a positive number: "<<m_structinfo.n_atoms);
    if ( hasAtomInfo() && m_structinfo.n_atoms != ntotatoms_from_atominfo )
      NCRYSTAL_THROW2(BadInput,"Inconsistent total number of atoms deduced from StructureInfo and AtomInfo ("
                      <<m_structinfo.n_atoms<<" vs. "<<ntotatoms_from_atominfo<<")");

    checkAndCompleteLattice( m_structinfo.spacegroup, m_structinfo.lattice_a,
                             m_structinfo.lattice_b, m_structinfo.lattice_c );

    if ( ! ( m_structinfo.alpha > 0 && m_structinfo.alpha < 180 &&
             m_structinfo.beta  > 0 && m_structinfo.beta  < 180 &&
             m_structinfo.gamma > 0 && m_structinfo.gamma < 180 ) ) {
      NCRYSTAL_THROW(BadInput,"Lattice angles must all be >0 and <180 degrees.");
    }

    if ( ! ( m_structinfo.alpha > kPi || m_structinfo.beta >kPi || m_structinfo.gamma > kPi ) ) {
      NCRYSTAL_THROW(BadInput,"Lattice angles in radians where degrees was expected.");
    }

  }

  if (hasXSectFree()) {
    nc_assert_always(getXSectFree()>=0.0&&getXSectFree()<1e20);
  }

  if (hasXSectAbsorption()) {
    nc_assert_always(getXSectAbsorption()>=0.0&&getXSectAbsorption()<1e20);
  }

  if (hasTemperature()) {
    nc_assert_always(getTemperature()>0.0&&getTemperature()<1e6);
  }

  if (hasGlobalDebyeTemperature()) {
    nc_assert_always(getGlobalDebyeTemperature()>0.0&&getGlobalDebyeTemperature()<1e6);
  }

  if (hasHKLInfo()) {
    nc_assert_always( m_hkl_dlower < m_hkl_dupper);
    nc_assert_always( m_hkl_dlower > 0.0 );
    nc_assert_always( m_hkl_dupper > 0.0 );
  }

  if ( hasDensity() ) {
    nc_assert_always(m_density>0.0);
    //TODO: Now that we know the masses, we can check consistency with atominfo/composition and number density!
  }
  if ( hasNumberDensity() ) {
    nc_assert_always(m_numberdensity>0.0);
  }

  if (hasDynamicInfo()) {
    if ( !hasTemperature() )
      NCRYSTAL_THROW(BadInput,"temperature info is required whenever dynamic info is supplied");

    //Check consistency of temperature fields:
    for (auto& di : getDynamicInfoList()) {
      if ( di->temperature() != this->getTemperature() )
        NCRYSTAL_THROW(BadInput,"temperature info on DynamicInfo (DI_VDOS) object"
                       " is different than temperature on owning Info object!");

      //NB: for DI_ScatKnlDirect we can't validate temperature field on SABData
      //object here, since it is built on-demand. However, we do it above in
      //the ensureBuildThenReturnSAB() method.

      auto di_vdos = dynamic_cast<const DI_VDOS*>(di.get());
      if ( di_vdos && di_vdos->vdosData().temperature() != di_vdos->temperature() )
        NCRYSTAL_THROW(BadInput,"temperature info on VDOSData object provided by DI_VDOS object"
                       " is different than temperature on DI_VDOS object itself!");

      auto di_vdosdebye = dynamic_cast<const DI_VDOSDebye*>(di.get());
      if ( di_vdosdebye ) {
        if ( !hasAnyDebyeTemperature() )
          NCRYSTAL_THROW(BadInput,"DI_VDOSDebye added to Info object without Debye temperatures!");
        if ( di_vdosdebye->debyeTemperature() != getDebyeTemperatureByElement(di_vdosdebye->atom().index) )
          NCRYSTAL_THROW(BadInput,"Debye temperature on DI_VDOSDebye object different than the one provided on the owning Info object!");
      }
    }
  }

  if ( !isCrystalline() and !hasDynamicInfo() )
    NCRYSTAL_THROW(BadInput,"Non-crystalline materials must have dynamic info present.");

  std::map<IndexedAtomData,double> composition_map_fromdyninfo;
  if (hasDynamicInfo()) {
    for (auto& di : getDynamicInfoList()) {
      if (composition_map_fromdyninfo.count(di->atom()))
        NCRYSTAL_THROW2(BadInput,"Multiple dynamic info sections for \""<<di->atomData().description(false)<<"\" (AtomIndex "<<di->atom().index.value<<")");
      composition_map_fromdyninfo[ di->atom() ] = di->fraction();
    }
  }

  if (m_composition.empty() && !composition_map_fromdyninfo.empty()) {
    for ( auto& e : composition_map_fromdyninfo ) {
      CompositionEntry entry;
      entry.fraction = e.second;
      entry.atom = e.first;
      m_composition.push_back(std::move(entry));
    }
    composition_map_fromdyninfo.clear();
  }
  if (m_composition.empty() && !composition_map_fromatompos.empty()) {
    for ( auto& e : composition_map_fromatompos ) {
      CompositionEntry entry;
      entry.fraction = e.second;
      entry.atom = e.first;
      m_composition.push_back(std::move(entry));
    }
    composition_map_fromatompos.clear();
  }
  if ( hasComposition() ) {
    double ftot(0.0);
    for (const auto& e : getComposition()) {
      if (e.atom.atomDataSP==nullptr)
        NCRYSTAL_THROW2(BadInput,"Nullptr atomData provided");
      ftot += e.fraction;
      if ( e.fraction<=0 || e.fraction>1.0)
        NCRYSTAL_THROW2(BadInput,"invalid composition fraction for element \""<<e.atom.data().description()<<"\": "<<e.fraction);
    }
    if (ftot >= 1.000000001 || ftot<0.999999999)
      NCRYSTAL_THROW(BadInput,"invalid composition : fractions do not sum to unity");
    //Check consistency between sources of composition:
    auto checkComposition = [this](const std::map<IndexedAtomData,double>& comp2,const char* name2 )
                            {
                              if (comp2.empty())
                                return;
                              if (m_composition.size()!=comp2.size())
                                NCRYSTAL_THROW2(BadInput,"incompatible compositions specified in "<<name2<<" (different number of elements)");
                              for (const auto& e : m_composition) {
                                auto it = comp2.find(e.atom);
                                if (it==comp2.end())
                                  NCRYSTAL_THROW2(BadInput,"incompatible compositions specified in "<<name2<<" (element \""<<e.atom.data().description()
                                                  <<"\" not present everywhere [or specified via different AtomData instances!])");
                                if (!floateq(e.fraction,it->second))
                                  NCRYSTAL_THROW2(BadInput,"incompatible compositions specified in "<<name2<<" (fraction of element "<<e.atom.data().description()<<" not consistent)");
                              }
                            };
    checkComposition(composition_map_fromdyninfo,"DynInfo");
    checkComposition(composition_map_fromatompos,"Atomic Positions");
  }

  for (const auto& e : m_custom) {
    if (e.first.empty() || !contains_only(e.first,"ABCDEFGHIJKLMNOPQRSTUVWXYZ"))
      NCRYSTAL_THROW2(BadInput,"invalid custom section name: \""<<e.first
                      <<"\" (must be non-empty and contain only capitalised letters A-Z)");
  }

  //Setup display labels and m_atomDataSPs:
  if (!m_composition.empty()) {
    //TODO: We could combine the next two vectors into just:
    //std::vector<std::pair<IndexedAtomData,std::string>> That would be more
    //efficient, and we could even return const refs to IndexAtomData objects.
    m_displayLabels.clear();
    m_displayLabels.resize(m_composition.size());
    m_atomDataSPs.clear();
    m_atomDataSPs.resize(m_composition.size());

    //Order by atomdata, then index (so display labels won't depend on whether
    //or not the info factory did sensible sorting):
    std::vector<const IndexedAtomData*> v;
    v.reserve(m_composition.size());
    for ( const auto& e : m_composition )
      v.push_back( &e.atom );
    std::stable_sort(v.begin(),v.end(),
                     [](const IndexedAtomData* a,const IndexedAtomData* b)
                     {
                       if ( a->atomDataSP->getUniqueID() == b->atomDataSP->getUniqueID() )
                         return a->index.value < b->index.value;
                       return *a->atomDataSP < *b->atomDataSP;
                     });

    std::map<std::string,std::vector<AtomIndex>> lbl2indices;
    for ( const auto& e : v ) {
      const IndexedAtomData& iad = *e;
      const AtomData& ad = *iad.atomDataSP;
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
        m_displayLabels.at(e.second.front().value) = e.first;
      } else {
        for (auto&& ee: enumerate(e.second)) {
          m_displayLabels.at(ee.val.value) = e.first + "("_s + idx_to_alphalbl(static_cast<unsigned>(ee.idx)) + ")"_s;
        }
      }
    }
    for (const auto& dl : m_displayLabels)
      nc_assert_always(!dl.empty());

    for ( const auto& e : v ) {
      const IndexedAtomData& iad = *e;
      m_atomDataSPs.at(iad.index.value) = std::move(iad.atomDataSP);
    }

  }
}

void NC::Info::enableHKLInfo(double dlower, double dupper)
{
  ensureNoLock();
  m_hkl_dlower = dlower;
  m_hkl_dupper = dupper;
  nc_assert_always(hasHKLInfo());
}

void NC::Info::ensureNoLock()
{
  if (m_lock)
    NCRYSTAL_THROW(LogicError,"Modification of Info object after it is locked is forbidden");
}

NC::HKLList::const_iterator NC::Info::searchExpandedHKL(short h, short k, short l) const
{
  nc_assert_always(hasHKLInfo());
  nc_assert_always(hasExpandedHKLInfo());

  HKLList::const_iterator it(hklBegin()), itE(hklEnd());
  for(;it!=itE;++it)
  {
    for(unsigned i=0;i < it->multiplicity/2;i++)
    {
      nc_assert(it->eqv_hkl);
      if( (it->eqv_hkl[i*3]==h && it->eqv_hkl[i*3+1]==k && it->eqv_hkl[i*3+2]==l) ||
          (it->eqv_hkl[i*3]==-h && it->eqv_hkl[i*3+1]==-k && it->eqv_hkl[i*3+2]==-l) )
      {
        return it;
      }
    }
  }
  return itE;
}

//TODO: why not always provide eqv_hkl from .ncmat factories and remove
//demi_normals from the interface? The memory usage is 4 times lower and the
//added initialisation time is likely negligible.

double NC::Info::hklDMinVal() const
{
  nc_assert(hasHKLInfo());
  if (m_hkllist.empty())
    return kInfinity;
  return hklLast()->dspacing;
}

double NC::Info::hklDMaxVal() const
{
  nc_assert(hasHKLInfo());
  if (m_hkllist.empty())
    return kInfinity;
  return hklBegin()->dspacing;
}

double NC::Info::dspacingFromHKL( int h, int k, int l ) const
{
  if (!hasStructureInfo())
    NCRYSTAL_THROW(MissingInfo,"Info object lacks Structure information.");
  const StructureInfo & si = getStructureInfo();
  RotMatrix rec_lat = getReciprocalLatticeRot( si.lattice_a, si.lattice_b, si.lattice_c,
                                               si.alpha*kDeg, si.beta*kDeg, si.gamma*kDeg );
  return NC::dspacingFromHKL( h,k,l, rec_lat );
}

NC::DynamicInfo::DynamicInfo(double fr, NC::IndexedAtomData atom, double tt)
  : m_fraction(fr),
    m_atom(atom),
    m_temperature(tt)
{
  nc_assert_always(!!m_atom.atomDataSP);
}

bool NC::DI_ScatKnlDirect::hasBuiltSAB() const
{
  std::lock_guard<std::mutex> lock(m_mutex);
  return !! m_sabdata;
}

std::shared_ptr<const NC::SABData> NC::DI_ScatKnlDirect::ensureBuildThenReturnSAB() const
{
  std::lock_guard<std::mutex> lock(m_mutex);
  if ( ! m_sabdata ) {
    m_sabdata = buildSAB();
    nc_assert_always( !! m_sabdata );
    if ( m_sabdata->temperature() != this->temperature() )
        NCRYSTAL_THROW(BadInput,"temperature info on SABData object provided by DI_ScatKnlDirect object"
                       " is different than temperature on DI_ScatKnlDirect object itself!");
  }
  return m_sabdata;
}

double NC::Info::getDebyeTemperatureByElement( const AtomIndex& atomindex ) const
{
  if (m_debyetemp > 0.0)
    return m_debyetemp;//global
  if (!hasAnyDebyeTemperature() )
    NCRYSTAL_THROW2(BadInput,"getDebyeTemperatureByElement called but no Debye temperature is available");
  for ( const auto& ai : m_atomlist ) {
    if ( ai.atom.index == atomindex ) {
      nc_assert_always(ai.debye_temp > 0 );
      return ai.debye_temp;
    }
  }
  NCRYSTAL_THROW2(BadInput,"getDebyeTemperatureByElement called for AtomIndex \""
                  <<atomindex.value<<"\" which was not found in this material");
}

NC::DynamicInfo::~DynamicInfo() = default;//virtual destructors must be implemented despite being abstract!
NC::DI_Sterile::~DI_Sterile() = default;
NC::DI_FreeGas::~DI_FreeGas() = default;
NC::DI_ScatKnl::~DI_ScatKnl() = default;
NC::DI_ScatKnlDirect::~DI_ScatKnlDirect() = default;
NC::DI_VDOS::~DI_VDOS() = default;
NC::DI_VDOSDebye::~DI_VDOSDebye() = default;

unsigned NC::Info::countCustomSections(const NC::Info::CustomSectionName& sectionname ) const
{
  unsigned i = 0;
  for (const auto& e: m_custom) {
    if (e.first==sectionname)
      ++i;
  }
  return i;
}

const NC::Info::CustomSectionData& NC::Info::getCustomSection( const NC::Info::CustomSectionName& name,
                                                               unsigned index ) const
{
  unsigned i = 0;
  for (const auto& e: m_custom) {
    if (e.first!=name)
      continue;
    if (index==i)
      return e.second;
    ++i;
  }
  NCRYSTAL_THROW2(MissingInfo,"Call to Info::getCustomSectionData requested the section "<<name
                  <<" with index="<<index<<" but info does not have at least "<<index+1
                  <<" such entries. Check with countCustomSections(..) before calling this method.");
}
