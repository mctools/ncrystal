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

#include "NCrystal/factories/NCFactRequests.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/cfgutils/NCCfgManip.hh"
#include "NCrystal/interfaces/NCSCOrientation.hh"
namespace NC = NCrystal;
namespace NCF = NCrystal::FactImpl;

namespace NCRYSTAL_NAMESPACE {
  namespace FactImpl {
    namespace {
      using Cfg::CfgManip;

      bool sr_varFilter(Cfg::detail::VarId varid)
      {
        auto gr = Cfg::varGroup(varid);
        return ( gr == Cfg::VarGroupId::ScatterBase
                 || gr == Cfg::VarGroupId::ScatterExtra );
      }

      bool ar_varFilter(Cfg::detail::VarId varid)
      {
        return Cfg::varGroup(varid) == Cfg::VarGroupId::Absorption;
      }

      void sr_checkParamConsistency( const Cfg::CfgData& data )
      {
        CfgManip::checkParamConsistency_ScatterBase( data );
        CfgManip::checkParamConsistency_ScatterExtra( data );
      }

      void ar_checkParamConsistency( const Cfg::CfgData& data )
      {
        CfgManip::checkParamConsistency_Absorption( data );
      }

      // static bool def_varFilter(Cfg::detail::VarId);
      // static void def_checkParamConsistency(const Cfg::CfgData&);

    }
  }
}

NCF::detail::ProcessRequestData::ParamDefs NCF::ScatterRequest::paramDefs()
{
  detail::ProcessRequestData::ParamDefs pd;
  pd.varFilter = sr_varFilter;
  pd.checkParamConsistency = sr_checkParamConsistency;
  return pd;
}

NCF::detail::ProcessRequestData::ParamDefs NCF::AbsorptionRequest::paramDefs()
{
  detail::ProcessRequestData::ParamDefs pd;
  pd.varFilter = ar_varFilter;
  pd.checkParamConsistency = ar_checkParamConsistency;
  return pd;
}

void NCF::ScatterRequest::checkParamConsistency() const
{
  sr_checkParamConsistency( rawCfgData() );
}

void NCF::AbsorptionRequest::checkParamConsistency() const
{
  ar_checkParamConsistency( rawCfgData() );
}

void NCF::InfoRequest::checkParamConsistency() const
{
  CfgManip::checkParamConsistency_Info( rawCfgData() );
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
                   [](Cfg::detail::VarId varid)
                   { return ( Cfg::varGroup(varid)
                              == Cfg::VarGroupId::Info ); } );
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

std::string NCF::AbsorptionRequest::get_absnfactory() const
{
  return CfgManip::get_absnfactory(rawCfgData()).to_string();
}

NCF::ScatterRequest NCF::ScatterRequest::createChildRequest( unsigned ichild ) const
{
  return { internal_t(), m_data.createChildRequest( ichild ) };
}

NCF::ScatterRequest
NCF::ScatterRequest::modified( const std::string& str) const
{
  return { internal_t(), m_data.modified( str ) };
}

NCF::ScatterRequest
NCF::ScatterRequest::modified( const char* cstr ) const
{
  return { internal_t(), m_data.modified( cstr ) };
}

NCF::AbsorptionRequest
NCF::AbsorptionRequest::createChildRequest( unsigned ichild ) const
{
  return { internal_t(), m_data.createChildRequest( ichild ) };
}

NCF::AbsorptionRequest
NCF::AbsorptionRequest::modified( const std::string& str) const
{
  return { internal_t(), m_data.modified( str ) };
}

NCF::AbsorptionRequest
NCF::AbsorptionRequest::modified( const char* cstr ) const
{
  return { internal_t(), m_data.modified( cstr ) };
}

NCF::ScatterRequest NCF::ScatterRequest::cloneThinned() const
{
  return { internal_t(), m_data.cloneThinned() };
}

NCF::AbsorptionRequest NCF::AbsorptionRequest::cloneThinned() const
{
  return { internal_t(), m_data.cloneThinned() };
}

