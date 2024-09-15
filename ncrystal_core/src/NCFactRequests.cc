////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
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

namespace NCRYSTAL_NAMESPACE {
  namespace FactImpl {
    namespace {
      using Cfg::CfgManip;
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
    detail::validateMatCfgState(cfg);
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


bool NCF::ScatterRequest::varIsApplicable(Cfg::detail::VarId varid)
{
  auto gr = Cfg::varGroup(varid);
  return ( gr == Cfg::VarGroupId::ScatterBase || gr == Cfg::VarGroupId::ScatterExtra );
}

bool NCF::AbsorptionRequest::varIsApplicable(Cfg::detail::VarId varid)
{
  return Cfg::varGroup(varid) == Cfg::VarGroupId::Absorption;
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
