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

#include "NCrystal/factories/NCFactRequestsImpl.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/cfgutils/NCCfgManip.hh"

namespace NC = NCrystal;
namespace NCFD = NCrystal::FactImpl::detail;

namespace NCRYSTAL_NAMESPACE {
  namespace FactImpl {
    namespace {
      using Cfg::CfgManip;
    }
  }
}

void NCFD::validateMatCfgState( const MatCfg& cfg ) {
  //Validate MatCfg object ok for XXXRequest constructor:
  if ( !cfg.isTrivial() )
    NCRYSTAL_THROW(BadInput,"Only trivial MatCfg objects can be passed to constructors of Request objects.");
  if ( cfg.isThinned() )
    NCRYSTAL_THROW(BadInput,"Thinned MatCfg objects can not be passed to constructors of Request objects.");
  nc_assert_always( !cfg.isMultiPhase() );
  nc_assert_always( cfg.getPhaseChoices().empty() );
  nc_assert_always( cfg.get_density() == DensityState() );
}

NCFD::ProcessRequestData NCFD::ProcessRequestData::cloneThinned() const
{
  ProcessRequestData res{ no_init };
  res.m_data = m_data;
  res.m_infoUID = m_infoUID;
  res.m_dataSourceName = m_dataSourceName;
  return res;
}

const NC::Info& NCFD::ProcessRequestData::info() const
{
  //Todo CalcError is not really a great type for this.
  if ( m_infoPtr == nullptr )
    NCRYSTAL_THROW(CalcError,"Do not use thinned ScatterRequest or"
                   " Absorptionrequest objects to access Info objects.");
  return *m_infoPtr;
}

NC::InfoPtr NCFD::ProcessRequestData::infoPtr() const {
  //Todo CalcError is not really a great type for this.
  if ( m_infoPtr == nullptr )
    NCRYSTAL_THROW(CalcError,"Do not use thinned ScatterRequest or"
                   " Absorptionrequest objects to access Info objects.");
  return m_infoPtr;
}


NCFD::ProcessRequestData::ProcessRequestData( internal_t,
                                              InfoPtr infoptr,
                                              const Cfg::CfgData* opt_data,
                                              ParamDefs pd )
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
    m_dataSourceName(m_infoPtr->getDataSourceName()),
    m_paramDefs( pd )
{
  nc_assert( m_paramDefs.varFilter != nullptr );
  nc_assert( m_paramDefs.checkParamConsistency != nullptr );

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
  CfgManip::apply( m_data, infoptr->getCfgData(), m_paramDefs.varFilter );
  if ( opt_data )
    CfgManip::apply( m_data, *opt_data, m_paramDefs.varFilter );
  m_paramDefs.checkParamConsistency( rawCfgData() );
}

NCFD::ProcessRequestData::ProcessRequestData( const MatCfg& cfg,
                                              ParamDefs pd )
  : ProcessRequestData( internal_t(),
                        [&cfg]()
                        {
                          validateMatCfgState(cfg);
                          return FactImpl::create(InfoRequest(cfg));
                        }(),
                        ( cfg.isTrivial()
                          ? &cfg.rawCfgData()
                          : nullptr),//of course, validateMatCfgState above
                                     //ensures cfg.isTrivial()==true, but it
                                     //might get invoked AFTER this line
                        pd
                        )
{
}

NCFD::ProcessRequestData NCFD::ProcessRequestData::createChildRequest( unsigned ichild ) const
{
  auto nchildren = isMultiPhase() ? info().getPhases().size() : 0;
  if ( ichild >= nchildren )
    NCRYSTAL_THROW2(BadInput,"createChildRequest index out of range (ichild="<<ichild<<", nchildren="<<nchildren<<")");
  auto info_child = info().getPhases().at(ichild).second;
  auto child_request = ProcessRequestData( info_child, m_paramDefs );
  CfgManip::apply( child_request.m_data, m_data );//cfg settings trickle down
                                                  //and override those on
                                                  //children [important to do
                                                  //this here due to the
                                                  //detail_copyUnderlying(..)
                                                  //usage).
  return child_request;
}

bool NCFD::ProcessRequestData::cmpDataLT( const ProcessRequestData& o ) const
{
  return ( m_dataSourceName.str() == o.m_dataSourceName.str()
           ? CfgManip::lessThan( m_data, o.m_data )
           : m_dataSourceName.str() < o.m_dataSourceName.str() );
}

bool NCFD::ProcessRequestData::cmpDataEQ( const ProcessRequestData& o ) const
{
  return ( m_dataSourceName.str() == o.m_dataSourceName.str()
           ? CfgManip::equal( m_data, o.m_data )
           : false );
}

void NCFD::ProcessRequestData::stream( std::ostream &os ) const
{
  os << m_dataSourceName << ";...";//;... represents the unknown Info params
  if ( !CfgManip::empty(m_data) ) {
    os << ';';
    CfgManip::stream(m_data,os);
  }
}

void NCFD::ProcessRequestData::streamParamsOnly( std::ostream &os ) const
{
  if ( !CfgManip::empty(m_data) )
    CfgManip::stream(m_data,os);
}

NCFD::ProcessRequestData NCFD::ProcessRequestData::modified( internal_t,
                                                             const char* strdata,
                                                             std::size_t len ) const
{
  StrView sv(strdata,len);
  Cfg::CfgData tmpdata;
  auto toplvlvars = CfgManip::applyStrCfg( tmpdata, sv );
  auto varNotApplicaple = [this](Cfg::detail::VarId varid){ return !m_paramDefs.varFilter(varid); };
  if ( !toplvlvars.empty() || CfgManip::filterSelectsAny( tmpdata, varNotApplicaple ) )
    NCRYSTAL_THROW2(BadInput,"Invalid cfgstr passed to Request::modified function: \""<<sv
                    <<"\" (only settings applicable to the process type are allowed in this context)");
  ProcessRequestData res( *this );
  CfgManip::apply( res.m_data, tmpdata );
  return res;
}

NCFD::ProcessRequestData NCFD::ProcessRequestData::modified( const std::string& str ) const
{
  return modified( internal_t(), str.c_str(), str.size() );
}

NCFD::ProcessRequestData NCFD::ProcessRequestData::modified( const char* cstr ) const
{
  StrView sv(cstr);
  return modified( internal_t(), sv.data(), sv.size() );
}
