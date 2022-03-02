////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2022 NCrystal developers                                   //
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
#include "NCrystal/internal/NCMath.hh"

namespace NC=NCrystal;

namespace NCrystal {
  static_assert(std::is_nothrow_move_constructible<AtomInfo>::value,"");
}

//TODO: why not always provide eqv_hkl from .ncmat factories and remove
//demi_normals from the interface? The memory usage is 4 times lower and the
//added initialisation time is likely negligible.

double NC::Info::hklDMinVal() const
{
  singlePhaseOnly(__func__);
  nc_assert(hasHKLInfo());
  auto& hklList = m_data->hklList();
  if (hklList.empty())
    return kInfinity;
  return std::prev(hklList.end())->dspacing;
}

double NC::Info::hklDMaxVal() const
{
  singlePhaseOnly(__func__);
  nc_assert(hasHKLInfo());
  auto& hklList = m_data->hklList();
  if (hklList.empty())
    return kInfinity;
  return hklList.front().dspacing;
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
  singlePhaseOnly(__func__);
  unsigned i = 0;
  for (const auto& e: m_data->custom) {
    if (e.first==sectionname)
      ++i;
  }
  return i;
}

const NC::Info::CustomSectionData& NC::Info::getCustomSection( const NC::Info::CustomSectionName& name,
                                                               unsigned index ) const
{
  singlePhaseOnly(__func__);
  unsigned i = 0;
  for (const auto& e: m_data->custom) {
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

std::string NC::Info::toString(StateOfMatter som)
{
  switch(som) {
  case StateOfMatter::Unknown: return "Unknown"_s;
  case StateOfMatter::Solid: return "Solid"_s;
  case StateOfMatter::Gas: return "Gas"_s;
  case StateOfMatter::Liquid: return "Liquid"_s;
  default:
    nc_assert_always(false);
    return ""_s;
  };
}

const NC::Info::PhaseList& NC::detail::getEmptyPL()
{
  static Info::PhaseList pl;
  return pl;
}

void NC::Info::singlePhaseOnlyRaiseError(const char*fctname) const
{
  nc_assert(isMultiPhase());
  NCRYSTAL_THROW2(LogicError,"Info::"<<fctname<<" should only be called on single-phase Info objects");
}


void NC::AtomInfo::detail_setupLink( DynamicInfo* di )
{
  nc_assert_always(di!=nullptr);
  //verify not set up already
  nc_assert_always(m_dyninfo==nullptr);
  nc_assert_always(di->m_atomInfo==nullptr);
  //set up:
  di->m_atomInfo = this;
  m_dyninfo = di;
}

NC::HKLList::const_iterator NC::Info::searchExpandedHKL(short h, short k, short l) const
{
  singlePhaseOnly(__func__);
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

NC::SigmaAbsorption NC::Info::getXSectAbsorption() const
{
  StableSum sum;
  for ( auto& di : getComposition() )
    sum.add( di.fraction * di.atom.data().captureXS().dbl() );
  return SigmaAbsorption{ DoValidate, sum.sum() };
}

NC::SigmaFree NC::Info::getXSectFree() const
{
  StableSum sum;
  for ( auto& di : getComposition() )
    sum.add( di.fraction * di.atom.data().freeScatteringXS().dbl() );
  return SigmaFree{ DoValidate, sum.sum() };
}

NC::ScatLenDensity NC::Info::getSLD() const
{
  StableSum sum;
  for ( auto& di : getComposition() )
    sum.add( di.fraction * di.atom.data().coherentScatLen() );
  return ScatLenDensity{ DoValidate, sum.sum() * getNumberDensity().dbl() * 100.0 };
}

NC::AtomMass NC::Info::getAverageAtomMass() const
{
  StableSum sum;
  for ( auto& di : getComposition() )
    sum.add( di.fraction * di.atom.data().averageMassAMU().dbl() );
  return AtomMass{ DoValidate, sum.sum() };
}

void NC::Info::Data::doInitHKLList() const
{
  nc_assert(hkl_dlower_and_dupper.has_value());
  nc_assert(hkl_ondemand_fct!=nullptr);
  //Do actual work (time consuming)
  auto res = hkl_ondemand_fct(hkl_dlower_and_dupper.value());
  //Now we must register the result - unless someone else beat us too it (using
  //global mtx is OK for the few instructions below):
  static std::mutex s_mtx;
  NCRYSTAL_LOCK_GUARD(s_mtx);
  if (!detail_hkllist_needs_init.load())
    return;//someone beat us to it - return (discarding our own result)
  detail_hklList = std::move(res);
  detail_hkllist_needs_init = false;
  if ( detail_braggthreshold.load() == -1.0 )
    detail_braggthreshold = detail_hklList.empty() ? 0.0 : detail_hklList.front().dspacing * 2.0;
}

NC::Optional<NC::NeutronWavelength> NC::Info::getBraggThreshold() const
{
  singlePhaseOnly(__func__);
  auto& data = *m_data;
  if( !data.hkl_dlower_and_dupper.has_value() )
    return NullOpt;
  double bt = data.detail_braggthreshold.load();
  if ( bt > 0.0 )
    return NeutronWavelength{bt};
  if ( bt == 0.0 )
    return NullOpt;
  //Needs init, try to avoid full init:
  for ( auto dlow : { 1.5, 0.75 } ) {
    if ( dlow >= data.hkl_dlower_and_dupper.value().second )
      continue;//trivially won't select any planes
    if ( data.detail_braggthreshold.load() >= 0.0 )
      return getBraggThreshold();//was filled (by ourselves in a previous loop or by other caller)
    if ( dlow <= data.hkl_dlower_and_dupper.value().first ) {
      //Break and do full init below
      break;
    } else {
      hklListPartialCalc( dlow );//will set data.detail_braggthreshold if dlow low enough
    }
  }

  //Fall back to full init if not done already:
  data.hklList();
  return getBraggThreshold();
}

NC::Optional<NC::HKLList> NC::Info::hklListPartialCalc( Optional<double> dlower,
                                                        Optional<double> dupper ) const
{
  singlePhaseOnly(__func__);
  auto& data = *m_data;
  nc_assert(data.hkl_dlower_and_dupper.has_value());
  if ( !data.hkl_ondemand_fct )
    return NullOpt;
  double dlow = ncmax( dlower.has_value() ? dlower.value() : data.hkl_dlower_and_dupper.value().first,
                       data.hkl_dlower_and_dupper.value().first );
  double dupp = ncmin( dupper.has_value() ? dupper.value() : data.hkl_dlower_and_dupper.value().second,
                       data.hkl_dlower_and_dupper.value().second );
  if ( ! ( dlow <= dupp ) || std::isnan(dlow) || std::isnan(dupp)  )
    NCRYSTAL_THROW2(BadInput,"hklListPartialCalc got invalid dspacing range request: ["<<dlow
                    <<", "<<dupp<<"] (once constrained to ["<<data.hkl_dlower_and_dupper.value().first
                    <<", "<<data.hkl_dlower_and_dupper.value().second<<"])");
  auto hklList = data.hkl_ondemand_fct( PairDD(dlow,dupp) );
  if ( !hklList.empty() && data.detail_braggthreshold.load() == -1.0 )
    data.detail_braggthreshold = hklList.front().dspacing * 2.0;
  return hklList;
}

#include "NCrystal/internal/NCLatticeUtils.hh"//Needed for dspacingFromHKL
double NC::Info::dspacingFromHKL( int h, int k, int l ) const
{
  singlePhaseOnly(__func__);
  if (!hasStructureInfo())
    NCRYSTAL_THROW(MissingInfo,"Info object lacks Structure information.");
  const StructureInfo & si = getStructureInfo();
  RotMatrix rec_lat = getReciprocalLatticeRot( si.lattice_a, si.lattice_b, si.lattice_c,
                                               si.alpha*kDeg, si.beta*kDeg, si.gamma*kDeg );
  return ::NC::dspacingFromHKL( h,k,l, rec_lat );
}

