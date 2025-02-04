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

#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/internal/utils/NCMath.hh"
namespace NC=NCrystal;

namespace NCRYSTAL_NAMESPACE {
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

namespace NCRYSTAL_NAMESPACE {
  namespace {
    template<class T>
    bool atomic_setValueIfHasValue(std::atomic<T>& av, T new_value, T value_required_for_set ) {
      //Set av to new_value if current value is value_required_for_set.
      T old_value = av.load();
      do {
        if ( old_value != value_required_for_set )
          return false;
        //bool compare_exchange_weak( T& expected, T desired, ... )
      } while ( !av.compare_exchange_weak( old_value, new_value ) );
      return true;
    }
  }
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

  //Take this chance to update Bragg threshold / HKLInfoType fields:
  double bt = detail_hklList.empty() ? 0.0 : detail_hklList.front().dspacing * 2.0;
  auto tp = enumAsInt(detail_hklList.empty() ? HKLInfoType::Minimal : detail_hklList.front().type() );
  atomic_setValueIfHasValue( detail_braggthreshold, bt, -1.0 );
  atomic_setValueIfHasValue( detail_hklInfoType, tp, hKLInfoTypeInt_unsetval );

  detail_hkllist_needs_init = false;
}

NC::HKLInfoType NC::Info::hklInfoType() const
{
  singlePhaseOnly(__func__);
  auto& data = *m_data;
  if( !data.hkl_dlower_and_dupper.has_value() ) {
    hklList();//This triggers exception
    return HKLInfoType::Minimal;//not used
  }
  auto tp = data.detail_hklInfoType.load();
  if ( tp != Data::hKLInfoTypeInt_unsetval )
    return static_cast<HKLInfoType>(tp);

  //trigger init and try again:
  getBraggThreshold();

  tp = data.detail_hklInfoType.load();
  nc_assert( tp != Data::hKLInfoTypeInt_unsetval );
  return static_cast<HKLInfoType>(tp);
}


NC::Optional<NC::NeutronWavelength> NC::Info::getBraggThreshold() const
{
  singlePhaseOnly(__func__);

  auto& data = *m_data;
  if( !data.hkl_dlower_and_dupper.has_value() )
    return NullOpt;

  auto retval = []( double bt )
  {
    Optional<NeutronWavelength> rv;
    if ( bt > 0 )
      rv.emplace( bt );
    return rv;
  };

  double bt = data.detail_braggthreshold.load();
  if ( bt >= 0.0 )
    return retval(bt);//already calculated

  //Needs init, try to avoid full init:
  for ( auto dlow : { 5.0, 1.5, 0.75 } ) {
    if ( (bt=data.detail_braggthreshold.load()) >= 0.0 )
      return retval(bt);//was filled (by ourselves in a previous loop or by concurrent call)
    if ( dlow > data.hkl_dlower_and_dupper.value().second )
      continue;//trivially won't select any planes
    if ( dlow <= data.hkl_dlower_and_dupper.value().first ) {
      //Better to break and do full init below
      break;
    } else {
      auto hl = hklListPartialCalc( dlow );//will set data.detail_braggthreshold if dlow low enough
    }
  }

  if ( (bt=data.detail_braggthreshold.load()) >= 0.0 )
    return retval(bt);//was filled (by ourselves above or by other caller)

  //Fall back to full init if not done already:
  data.hklList();
  nc_assert(data.detail_braggthreshold.load()>=0.0);
  return retval(data.detail_braggthreshold.load());
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
  if ( !hklList.empty() && !dupper.has_value() ) {
    //Take this chance to update Bragg threshold / HKLInfoType fields:
    double bt = hklList.front().dspacing * 2.0;
    auto tp = enumAsInt( hklList.front().type() );
    atomic_setValueIfHasValue( data.detail_braggthreshold, bt, -1.0 );
    atomic_setValueIfHasValue( data.detail_hklInfoType, tp, Data::hKLInfoTypeInt_unsetval );
  }
  return hklList;
}

#include "NCrystal/internal/utils/NCLatticeUtils.hh"//Needed for dspacingFromHKL
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

