////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2023 NCrystal developers                                   //
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

#include "NCrystal/NCFactRequests.hh"
#include "NCrystal/NCMatCfg.hh"
#include "NCrystal/NCFactImpl.hh"
#include "NCrystal/internal/NCCfgManip.hh"
#include "NCrystal/NCSCOrientation.hh"
namespace NC = NCrystal;
namespace NCF = NCrystal::FactImpl;

namespace NCrystal {
  namespace FactImpl {
    namespace {
      using Cfg::CfgManip;
      void validateMatCfgState( const MatCfg& cfg ) {
        //Validate MatCfg object ok for XXXRequest constructor:
        if ( !cfg.isTrivial() )
          NCRYSTAL_THROW(BadInput,"Only trivial MatCfg objects can be passed to constructors of Request objects.");
        if ( cfg.isThinned() )
          NCRYSTAL_THROW(BadInput,"Thinned MatCfg objects can not be passed to constructors of Request objects.");
        nc_assert( !cfg.isMultiPhase() );
        nc_assert( cfg.getPhaseChoices().empty() );
        nc_assert( cfg.get_density() == DensityState() );
      }
    }
  }
}

void NCF::InfoRequest::checkParamConsistency() const
{
  CfgManip::checkParamConsistency_Info( rawCfgData() );
}

void NCF::ScatterRequest::checkParamConsistency() const
{
  CfgManip::checkParamConsistency_ScatterBase( rawCfgData() );
  CfgManip::checkParamConsistency_ScatterExtra( rawCfgData() );
}

void NCF::AbsorptionRequest::checkParamConsistency() const
{
  CfgManip::checkParamConsistency_Absorption( rawCfgData() );
}

NCF::InfoRequest::InfoRequest( const MatCfg& cfg )
  : m_textDataSP( [&cfg]() {
    validateMatCfgState(cfg);
    auto tdsp = cfg.textDataSP();
    nc_assert(tdsp!=nullptr);
    return tdsp;
  }()),
    m_textDataUID( m_textDataSP->dataUID() ),
    m_dataSourceName( cfg.getDataSourceName() )
{
  CfgManip::apply( m_data,
                   cfg.rawCfgData(),
                   [](Cfg::detail::VarId varid){ return Cfg::varGroup(varid) == Cfg::VarGroupId::Info; } );
  checkParamConsistency();
}

bool NCF::InfoRequest::cmpDataLT( const InfoRequest& o ) const
{
  return ( m_dataSourceName.str() == o.m_dataSourceName.str()
           ? CfgManip::lessThan( m_data, o.m_data )
           : m_dataSourceName.str() < o.m_dataSourceName.str() );
}

bool NCF::InfoRequest::cmpDataEQ( const InfoRequest& o ) const
{
  return ( m_dataSourceName.str() == o.m_dataSourceName.str()
           ? CfgManip::equal( m_data, o.m_data )
           : false );
}

void NCF::InfoRequest::stream( std::ostream &os ) const
{
  os << m_dataSourceName;
  if ( !CfgManip::empty(m_data) ) {
    os << ';';
    CfgManip::stream(m_data,os);
  }
}

void NCF::InfoRequest::streamParamsOnly( std::ostream &os ) const
{
  if ( !CfgManip::empty(m_data) )
    CfgManip::stream(m_data,os);
}

NC::Temperature NCF::InfoRequest::get_temp() const { return CfgManip::get_temp(m_data); }
double NCF::InfoRequest::get_dcutoff() const { return CfgManip::get_dcutoff(m_data); }
double NCF::InfoRequest::get_dcutoffup() const { return CfgManip::get_dcutoffup(m_data); }
std::string NCF::InfoRequest::get_infofactory() const { return CfgManip::get_infofactory( m_data ).to_string(); }
std::string NCF::InfoRequest::get_atomdb() const { return CfgManip::get_atomdb( m_data ).to_string(); }
std::vector<NC::VectS> NCF::InfoRequest::get_atomdb_parsed() const { return CfgManip::get_atomdb_parsed( m_data ); }


inline bool NCF::ScatterRequest::varIsApplicable(Cfg::detail::VarId varid)
{
  auto gr = Cfg::varGroup(varid);
  return ( gr == Cfg::VarGroupId::ScatterBase || gr == Cfg::VarGroupId::ScatterExtra );
}

inline bool NCF::AbsorptionRequest::varIsApplicable(Cfg::detail::VarId varid)
{
  return Cfg::varGroup(varid) == Cfg::VarGroupId::Absorption;
}

template<typename TRequest>
NCF::ProcessRequestBase<TRequest>::ProcessRequestBase( internal_t, InfoPtr infoptr, const Cfg::CfgData* opt_data )
  //We want to record just the underlying part of infoptr along with the
  //cfgdata. The reason for this is is that processes provide
  //cross-section-per-barn and thus can be shared between materials which only
  //differ by a density override. Unfortunately the call to
  //Info::detail_copyUnderlying(infoptr) below is somewhat fragile, since it
  //also discards other fields on Info::OverrideableDataFields, namely
  //overridden cfgdata and phaselist (in case of multiphase material). Now, the
  //cfgdata is reapplied manually below, but we could potentially have buggy
  //behaviour if the underlying and overriden phaselists would not simply differ
  //by a general density scale. Ideally, this is guaranteed by having the code
  //in NCInfoBuilder.cc only override phaselists in such a way, but we
  //nonetheless perform a few sanity checks below.
  : m_infoPtr(Info::detail_copyUnderlying(infoptr)),//(NB: don't std::move in this line!!)
    m_infoUID(m_infoPtr->getUniqueID()),
    m_dataSourceName(m_infoPtr->getDataSourceName())
{
  //Sanity checks that underlying and overridden phase lists only differ in a density scale.
  nc_assert( m_infoPtr->isMultiPhase() == infoptr->isMultiPhase() );
  if ( m_infoPtr.get() != infoptr.get() && m_infoPtr->isMultiPhase() ) {
    auto &pl1 = infoptr->getPhases();
    auto &pl2 = m_infoPtr->getPhases();
    nc_assert_always( pl1.size() == pl2.size() );
    for ( auto i : ncrange(pl1.size()) ) {
      nc_assert_always( pl1.at(i).first == pl2.at(i).first );
      nc_assert_always( pl1.at(i).second->detail_getUnderlyingUniqueID()
                        == pl2.at(i).second->detail_getUnderlyingUniqueID() );
      //NB: We could go ahead and verify stuff like relative density and cfgdata...
    }
  }

  //We take the cfg data from infoptr, NOT m_infoPtr (due to the
  //detail_copyUnderlying call above it was discarded from m_infoPtr):
  CfgManip::apply( m_data, infoptr->getCfgData(), TRequest::varIsApplicable );
  if ( opt_data )
    CfgManip::apply( m_data, *opt_data, TRequest::varIsApplicable );
  static_cast<const TRequest*>(this)->checkParamConsistency();
}

template<typename TRequest>
NCF::ProcessRequestBase<TRequest>::ProcessRequestBase( const MatCfg& cfg )
  : ProcessRequestBase( internal_t(),
                        [&cfg]()
                        {
                          validateMatCfgState(cfg);
                          return FactImpl::create(InfoRequest(cfg));
                        }(),
                        (cfg.isTrivial()?&cfg.rawCfgData():nullptr)//of course, validateMatCfgState above ensures
                                                                   //cfg.isTrivial()==true, but it might get invoked
                                                                   //AFTER this line
                        )
{
}

template<typename TRequest>
TRequest NCF::ProcessRequestBase<TRequest>::createChildRequest( unsigned ichild ) const
{
  auto nchildren = isMultiPhase() ? info().getPhases().size() : 0;
  if ( ichild >= nchildren )
    NCRYSTAL_THROW2(BadInput,"createChildRequest index out of range (ichild="<<ichild<<", nchildren="<<nchildren<<")");
  auto info_child = info().getPhases().at(ichild).second;
  auto child_request = TRequest( info_child );
  CfgManip::apply( child_request.m_data, m_data );//cfg settings trickle down
                                                  //and override those on
                                                  //children [important to do
                                                  //this here due to the
                                                  //detail_copyUnderlying(..)
                                                  //usage)
  return child_request;
}

template<typename TRequest>
bool NCF::ProcessRequestBase<TRequest>::cmpDataLT( const ProcessRequestBase& o ) const
{
  return ( m_dataSourceName.str() == o.m_dataSourceName.str()
           ? CfgManip::lessThan( m_data, o.m_data )
           : m_dataSourceName.str() < o.m_dataSourceName.str() );
}

template<typename TRequest>
bool NCF::ProcessRequestBase<TRequest>::cmpDataEQ( const ProcessRequestBase& o ) const
{
  return ( m_dataSourceName.str() == o.m_dataSourceName.str()
           ? CfgManip::equal( m_data, o.m_data )
           : false );
}

template<typename TRequest>
void NCF::ProcessRequestBase<TRequest>::stream( std::ostream &os ) const
{
  os << m_dataSourceName << ";...";//;... represents the unknown Info params
  if ( !CfgManip::empty(m_data) ) {
    os << ';';
    CfgManip::stream(m_data,os);
  }
}

template<typename TRequest>
void NCF::ProcessRequestBase<TRequest>::streamParamsOnly( std::ostream &os ) const
{
  if ( !CfgManip::empty(m_data) )
    CfgManip::stream(m_data,os);
}

template<typename TRequest>
TRequest NCF::ProcessRequestBase<TRequest>::modified( internal_t, const char*strdata, std::size_t len ) const
{
  StrView sv(strdata,len);
  Cfg::CfgData tmpdata;
  auto toplvlvars = CfgManip::applyStrCfg( tmpdata, sv );
  auto varNotApplicaple = [](const Cfg::detail::VarId varid){ return !TRequest::varIsApplicable(varid); };
  if ( !toplvlvars.empty() || CfgManip::filterSelectsAny( tmpdata, varNotApplicaple ) )
    NCRYSTAL_THROW2(BadInput,"Invalid cfgstr passed to Request::modified function: \""<<sv
                    <<"\" (only settings applicable to the process type are allowed in this context)");
  auto res = TRequest( *static_cast<const TRequest*>(this) );
  CfgManip::apply( res.m_data, tmpdata );
  return res;
}

template<typename TRequest>
TRequest NCF::ProcessRequestBase<TRequest>::modified( const std::string& str ) const
{
  return modified( internal_t(), str.c_str(), str.size() );
}

template<typename TRequest>
TRequest NCF::ProcessRequestBase<TRequest>::modified( const char* cstr ) const
{
  StrView sv(cstr);
  return modified( internal_t(), sv.data(), sv.size() );
}

int NCF::ScatterRequest::get_vdoslux() const { return CfgManip::get_vdoslux(rawCfgData()); }
bool NCF::ScatterRequest::get_coh_elas() const { return CfgManip::get_coh_elas(rawCfgData()); }
bool NCF::ScatterRequest::get_incoh_elas() const { return CfgManip::get_incoh_elas(rawCfgData()); }
bool NCF::ScatterRequest::get_sans() const { return CfgManip::get_sans(rawCfgData()); }
std::string NCF::ScatterRequest::get_inelas() const { return CfgManip::get_inelas(rawCfgData()).to_string(); }
std::string NCF::ScatterRequest::get_scatfactory() const { return CfgManip::get_scatfactory(rawCfgData()).to_string(); }
NC::MosaicityFWHM NCF::ScatterRequest::get_mos() const { return CfgManip::get_mos(rawCfgData()); }
NC::OrientDir NCF::ScatterRequest::get_dir1() const { return CfgManip::get_dir1(rawCfgData()); }
NC::OrientDir NCF::ScatterRequest::get_dir2() const { return CfgManip::get_dir2(rawCfgData()); }
double NCF::ScatterRequest::get_mosprec() const { return CfgManip::get_mosprec(rawCfgData()); }
double NCF::ScatterRequest::get_sccutoff() const { return CfgManip::get_sccutoff(rawCfgData()); }
double NCF::ScatterRequest::get_dirtol() const { return CfgManip::get_dirtol(rawCfgData()); }
const NC::LCAxis& NCF::ScatterRequest::get_lcaxis() const { return CfgManip::get_lcaxis(rawCfgData()); }
std::int_least32_t NCF::ScatterRequest::get_lcmode() const { return CfgManip::get_lcmode(rawCfgData()); }
NC::StrView NCF::ScatterRequest::get_ucnmode_str() const { return CfgManip::get_ucnmode_str(rawCfgData()); }
NC::Optional<NC::UCNMode> NCF::ScatterRequest::get_ucnmode() const { return CfgManip::get_ucnmode(rawCfgData()); }

bool NCF::ScatterRequest::isSingleCrystal() const
{
  return CfgManip::isSingleCrystal( rawCfgData() );
}

bool NCF::ScatterRequest::isLayeredCrystal() const
{
  return CfgManip::isLayeredCrystal( rawCfgData() );
}

NC::SCOrientation NCF::ScatterRequest::createSCOrientation() const
{
  return CfgManip::createSCOrientation<SCOrientation>(rawCfgData());
}

std::string NCF::AbsorptionRequest::get_absnfactory() const { return CfgManip::get_absnfactory(rawCfgData()).to_string(); }

//Explicit instantiation (needed since we wanted templated code in this non-header file):
template class NCF::ProcessRequestBase<NCF::ScatterRequest>;
template class NCF::ProcessRequestBase<NCF::AbsorptionRequest>;
