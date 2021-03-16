////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2021 NCrystal developers                                   //
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
#include <cstring>//for memcpy
#include <cstdlib>
namespace NC=NCrystal;

namespace NCrystal {

  static_assert(std::is_nothrow_move_constructible<AtomInfo>::value,"");

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
    return rh.atomData() < lh.atomData();
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
    const double pos_tolerance = 0.01;
    for (std::size_t i = 1; i < plist.size(); ++i) {
      const AtomInfo::Pos& p1 = plist.at(i-1);
      const AtomInfo::Pos& p2 = plist.at(i);
      if ( ncabs(p1[0]-p2[0])<pos_tolerance && ncabs(p1[1]-p2[1])<pos_tolerance && ncabs(p1[2]-p2[2])<pos_tolerance )
        NCRYSTAL_THROW2(BadInput,"The same atom position used more than once: ("<<p1[0]<<", "<<p1[1]<<", "<<p1[2]<<")");
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

  //sort lists for reproducibility:
  std::stable_sort(m_hkllist.begin(),m_hkllist.end(),dhkl_compare);
  std::stable_sort(m_atomlist.begin(),m_atomlist.end(),atominfo_compare);
  std::stable_sort(m_composition.begin(),m_composition.end(),
                   [](const CompositionEntry& a,const CompositionEntry& b)
                   {
                     return ( a.atom==b.atom
                              ? a.fraction > b.fraction
                              : a.atom < b.atom );
                   });
  std::stable_sort(m_dyninfolist.begin(),m_dyninfolist.end(),
                   [](const std::unique_ptr<DynamicInfo>& a,const std::unique_ptr<DynamicInfo>& b)
                   {
                     nc_assert( !!a && !!b );
                     return ( a->atom() == b->atom()
                              ? a->fraction() > b->fraction()
                              : a->atom() < b->atom() );
                   });

  //Check AtomInfo objects for consistency and standardise them a bit:
  AtomInfo::AtomPositions all_positions;
  all_positions.reserve(128);

  unsigned ntotatoms_from_atominfo(0);

  std::map<IndexedAtomData,double> composition_map_fromatompos;
  AtomList::iterator itAtm(m_atomlist.begin()), itAtmE(m_atomlist.end());
  for (;itAtm!=itAtmE;++itAtm)
  {
    //Map all atom positions to interval [0,1) (e.g. 1.0 becomes 0.0, -0.3
    //becomes 0.7, etc.). However at most move 1.0:
    for ( auto& pos : itAtm->m_pos )
      for ( unsigned i = 0; i < 3; ++i )
        atominfo_pos_remap( pos[i] );

    //Sort positions, ensuring a well defined ordering:
    std::stable_sort(itAtm->m_pos.begin(),itAtm->m_pos.end(),atominfo_pos_compare);

    if (!itAtm->numberPerUnitCell())
      NCRYSTAL_THROW(BadInput,"atominfo found for element with number_per_unit_cell=0!.");

    ntotatoms_from_atominfo += itAtm->numberPerUnitCell();

    //Debye temps and msd's should be specified for all or none:
    if ( itAtm->debyeTemp().has_value() != m_atomlist.front().debyeTemp().has_value() )
      NCRYSTAL_THROW(LogicError,"Inconsistency: Debye temperatures specified for some but not all AtomInfo objects.");
    if ( itAtm->msd().has_value() != m_atomlist.front().msd().has_value() )
      NCRYSTAL_THROW(LogicError,"Inconsistency: MSD values specified for some but not all AtomInfo objects.");

    composition_map_fromatompos[itAtm->atom()] += itAtm->numberPerUnitCell();
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

  HKLList::iterator ithkl = m_hkllist.begin();
  HKLList::iterator ithklE = m_hkllist.end();
  int has_demi_normals(-1);
  int has_eqv_hkl(-1);
  for (;ithkl!=ithklE;++ithkl) {
    if (! ithkl->demi_normals.empty() ) {
      ithkl->demi_normals.shrink_to_fit();
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
        if ( ! itN->as<Vector>().isUnitVector() )
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

  if ( m_structinfo.has_value() ) {
    StructureInfo& structinfo = m_structinfo.value();
    if ( ! (structinfo.volume > 0.0) )
      NCRYSTAL_THROW2(BadInput,"StructureInfo volume not a positive number: "<<structinfo.volume);
    if ( ! (structinfo.n_atoms > 0.0) )
      NCRYSTAL_THROW2(BadInput,"StructureInfo n_atoms not a positive number: "<<structinfo.n_atoms);
    if ( hasAtomInfo() && structinfo.n_atoms != ntotatoms_from_atominfo )
      NCRYSTAL_THROW2(BadInput,"Inconsistent total number of atoms deduced from StructureInfo and AtomInfo ("
                      <<structinfo.n_atoms<<" vs. "<<ntotatoms_from_atominfo<<")");

    checkAndCompleteLattice( structinfo.spacegroup, structinfo.lattice_a,
                             structinfo.lattice_b, structinfo.lattice_c );

    if ( ! ( structinfo.alpha > 0 && structinfo.alpha < 180 &&
             structinfo.beta  > 0 && structinfo.beta  < 180 &&
             structinfo.gamma > 0 && structinfo.gamma < 180 ) ) {
      NCRYSTAL_THROW(BadInput,"Lattice angles must all be >0 and <180 degrees.");
    }

    if ( ! ( structinfo.alpha > kPi || structinfo.beta >kPi || structinfo.gamma > kPi ) ) {
      NCRYSTAL_THROW(BadInput,"Lattice angles in radians where degrees was expected.");
    }

  }

  if (hasXSectFree()) {
    nc_assert_always(getXSectFree().get()>=0.0&&getXSectFree().get()<1e20);
  }

  if (hasXSectAbsorption()) {
    nc_assert_always(getXSectAbsorption().get()>=0.0&&getXSectAbsorption().get()<1e20);
  }

  if (hasTemperature()) {
    nc_assert_always(getTemperature().get()>0.0&&getTemperature().get()<1e6);
  }

  if (hasHKLInfo()) {
    nc_assert_always( hklDLower() <= hklDUpper() );
    nc_assert_always( hklDLower() > 0.0 );
    nc_assert_always( hklDUpper() > 0.0 );
  }

  if ( hasDensity() ) {
    nc_assert_always(getDensity().dbl()>0.0);
    //TODO: Now that we know the masses, we can check consistency with atominfo/composition and number density!
  }
  if ( hasNumberDensity() ) {
    nc_assert_always(getNumberDensity().dbl()>0.0);
  }

  if (hasDynamicInfo()) {
    if ( !hasTemperature() )
      NCRYSTAL_THROW(BadInput,"temperature info is required whenever dynamic info is supplied");

    //Check consistency of temperature fields (checking consistency of
    //DI_VDOSDebye Debye temperature fields below:
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
    }
  }

  if ( !isCrystalline() and !hasDynamicInfo() )
    NCRYSTAL_THROW(BadInput,"Non-crystalline materials must have dynamic info present.");

  std::map<IndexedAtomData,double> composition_map_fromdyninfo;
  if (hasDynamicInfo()) {
    for (auto& di : getDynamicInfoList()) {
      if (composition_map_fromdyninfo.count(di->atom()))
        NCRYSTAL_THROW2(BadInput,"Multiple dynamic info sections for \""<<di->atomData().description(false)<<"\" (AtomIndex "<<di->atom().index<<")");
      composition_map_fromdyninfo[ di->atom() ] = di->fraction();
    }
  }

  if (m_composition.empty() && !composition_map_fromdyninfo.empty()) {
    for ( auto& e : composition_map_fromdyninfo )
      m_composition.emplace_back( e.second, e.first );
    composition_map_fromdyninfo.clear();
  }
  if (m_composition.empty() && !composition_map_fromatompos.empty()) {
    for ( auto& e : composition_map_fromatompos )
      m_composition.emplace_back( e.second, e.first );
    composition_map_fromatompos.clear();
  }
  if ( hasComposition() ) {
    double ftot(0.0);
    for (const auto& e : getComposition()) {
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
        m_displayLabels.at(e.second.front().get()) = e.first;
      } else {
        for (auto&& ee: enumerate(e.second)) {
          m_displayLabels.at(ee.val.get()) = e.first + "("_s + idx_to_alphalbl(static_cast<unsigned>(ee.idx)) + ")"_s;
        }
      }
    }
    for (const auto& dl : m_displayLabels)
      nc_assert_always(!dl.empty());

    std::vector<OptionalAtomDataSP> tmp_atomDataSPs;
    tmp_atomDataSPs.resize(m_composition.size());

    for ( const auto& e : v ) {
      const IndexedAtomData& iad = *e;
      tmp_atomDataSPs.at(iad.index.get()) = iad.atomDataSP;
    }

    m_atomDataSPs.reserve(tmp_atomDataSPs.size());
    for ( auto&& e : tmp_atomDataSPs )
      m_atomDataSPs.emplace_back(std::move(e));

  }
  //Setup AtomInfo<->DynamicInfo links:
  if ( !m_atomlist.empty() && !m_dyninfolist.empty()) {
    nc_assert_always( m_dyninfolist.size() == m_atomlist.size() );
    unsigned nlinks(0);
    for ( auto& ai : m_atomlist ) {
      for ( auto& di : m_dyninfolist ) {
        if ( di->atom().index == ai.atom().index ) {
          ++nlinks;
          di->m_atomInfo = &ai;
          ai.m_dyninfo = di.get();
        }
      }
    }
    nc_assert_always(nlinks==m_atomlist.size());
  }
  //Verify consistent Debye temperatures in DI_VDOSDebye objects:
  for ( auto& di : m_dyninfolist ) {
    auto di_vdosdebye = dynamic_cast<const DI_VDOSDebye*>(di.get());
    if ( !di_vdosdebye )
      continue;
    if ( di_vdosdebye->m_atomInfo == nullptr )
      NCRYSTAL_THROW(BadInput,"DI_VDOSDebye can only be added on an Info object when AtomInfo objects are also added!");
    if ( ! di_vdosdebye->m_atomInfo->debyeTemp().has_value() )
      NCRYSTAL_THROW(BadInput,"AtomInfo object associated with DI_VDOSDebye object must have Debye temperature available!");
    if ( di_vdosdebye->m_atomInfo->debyeTemp().value() != di_vdosdebye->debyeTemperature() )
      NCRYSTAL_THROW(BadInput,"Associated AtomInfo and DI_VDOSDebye objects do not have the same Debye temperature specified!");
  }
}

void NC::Info::enableHKLInfo(double dlower, double dupper)
{
  ensureNoLock();
  nc_assert_always( !m_hkl_dlower_and_dupper.has_value() );
  m_hkl_dlower_and_dupper = std::make_pair(dlower,dupper);
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

NC::DynamicInfo::DynamicInfo(double fr, NC::IndexedAtomData atom, Temperature tt)
  : m_fraction(fr),
    m_atom(atom),
    m_temperature(tt)
{
}

bool NC::DI_ScatKnlDirect::hasBuiltSAB() const
{
  NCRYSTAL_LOCK_GUARD(m_mutex);
  return !! m_sabdata;
}

std::shared_ptr<const NC::SABData> NC::DI_ScatKnlDirect::ensureBuildThenReturnSAB() const
{
  NCRYSTAL_LOCK_GUARD(m_mutex);
  if ( ! m_sabdata ) {
    m_sabdata = buildSAB();
    nc_assert_always( !! m_sabdata );
    if ( m_sabdata->temperature() != this->temperature() )
        NCRYSTAL_THROW(BadInput,"temperature info on SABData object provided by DI_ScatKnlDirect object"
                       " is different than temperature on DI_ScatKnlDirect object itself!");
  }
  return m_sabdata;
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

NC::AtomInfo::AtomInfo( IndexedAtomData iad,
                        AtomPositions&& pos,
                        Optional<DebyeTemperature> dt,
                        Optional<double> msd )
  : m_iad(std::move(iad)),
    m_dt(std::move(dt)),
    m_msd(std::move(msd)),
    m_pos(std::move(pos))
{
  nc_assert_always( m_pos.size()<100000 && m_pos.size() < std::numeric_limits<unsigned>::max() );
  if ( m_pos.empty() )
    NCRYSTAL_THROW(BadInput,"Empty position list passed to AtomInfo constructor.");
  if ( m_msd.has_value() && !(m_msd.value()>0.0&&m_msd.value()<1e20) )
    NCRYSTAL_THROW2(BadInput,"Invalid msd value passed to AtomInfo constructor:" << m_msd.value() );
  if ( m_dt.has_value() && ! ( m_dt.value().dbl() >= 0.1 && m_dt.value().dbl() <= 1.0e6 ) )
    NCRYSTAL_THROW2(LogicError,"Invalid debye temperature value passed to AtomInfo constructor: " << m_dt.value());
}

