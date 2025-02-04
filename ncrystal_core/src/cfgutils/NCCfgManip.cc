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

#include "NCrystal/internal/cfgutils/NCCfgManip.hh"

namespace NC = NCrystal;

NC::Cfg::CfgData NC::Cfg::CfgManip::filter( const CfgData& data, VarIdFilter filter )
{
  CfgData filtered;
  for ( auto& e : data() ) {
    if ( filter( static_cast<VarId>(e.metaData()) ) )
      filtered().push_back( e );
  }
  return filtered;
}

void NC::Cfg::CfgManip::apply( CfgData& data, const CfgData& o, VarIdFilter filter )
{
  for ( auto& e : o() ) {
    if ( filter!=nullptr && !filter( static_cast<VarId>(e.metaData()) ) )
      continue;
    auto cloneVarBuf = [&e]() -> VarBuf { return e; };
    detail_setVar( data, e.metaData(), cloneVarBuf );
  }
}

NC::Cfg::TopLvlVarList NC::Cfg::CfgManip::applyStrCfg( CfgData& data, StrView str )
{
  nc_assert( str.has_value() );
  const char *errintro = "Syntax error in parameter string: ";
  TopLvlVarList toplvlVars;
  auto badchar_strrep = findForbiddenChar( str, forbidden_chars_non_multiphase,
                                           ExtraForbidOpt::RequireSimpleASCII );
  if ( badchar_strrep.has_value() )
    NCRYSTAL_THROW2(BadInput,errintro<<"forbidden "<<badchar_strrep.value()
                    <<" character encountered.");

  auto it = str.begin();
  auto itE = str.end();
  auto skipWS = [&it,&itE]() { while ( it!=itE && isWhiteSpace(*it) ) ++it; };
  auto skipWSAndSemicolons = [&it,&itE]() { while ( it!=itE && ( *it==';' || isWhiteSpace(*it)) ) ++it; };
  while (it!=itE) {
    skipWSAndSemicolons();//skipping leading semicolons is the same as skipping empty parts
    if ( it == itE )
      break;//done.
    //Expect variable name at lhs, ending at '=' or whitespace:
    auto it2 = it;
    while ( it2 != itE && *it2 != '=' && !isWhiteSpace(*it2) )
      ++it2;
    //Collect lhs parameter name (validate it below):
    auto par_str = StrView(it,std::distance(it,it2));
    //Advance iterator to where rhs value begins:
    it = it2;
    skipWS();//spaces before '=' sign
    if ( it==itE || *it!='=' ) {
      if ( par_str.trimmed()=="ignorefilecfg" ) {
        NCRYSTAL_THROW2(BadInput,errintro<<"The \"ignorefilecfg\" keyword is no longer supported.");
      } else {
        NCRYSTAL_THROW2(BadInput,errintro<<"missing \"=\" character.");
      }
    }
    nc_assert(*it=='=');
    ++it;//skip '=' sign
    skipWS();//spaces after '=' sign
    //Collect rhs value. Note that empty strings are allowed here.
    it2 = it;
    while ( it2 != itE && *it2 != ';' ) {
      if (*it2=='=')
        NCRYSTAL_THROW2(BadInput,errintro<<"unexpected \"=\" character.");
      ++it2;
    }
    //remove trailing spaces:
    auto it3 = it2;
    while ( it2 > it && isWhiteSpace(*std::prev(it2)) )
      --it2;
    auto val_str = StrView(it,std::distance(it,it2));
    nc_assert(!val_str.contains_any("="));
    nc_assert(!val_str.contains_any(";"));
    it = it3;
    //Decode parameter index (for efficiency we do no preliminary validation first):
    auto opt_varid = varIdFromName( par_str );
    if ( opt_varid.has_value() ) {
      //Real parameter:
      setVarByStr( data, opt_varid.value(), val_str );
      continue;
    }

    //Pseudo parameters:

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!! IMPORTANT                            !!!
    // !!! Any modification to pseudo or toplvl !!!
    // !!! parameters must be reflected in the  !!!
    // !!! docs in NCCfgVars.cc                 !!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    static constexpr auto sv_bragg = StrView::make("bragg");
    static constexpr auto sv_coh_elas = StrView::make("coh_elas");
    static constexpr auto sv_incoh_elas = StrView::make("incoh_elas");
    static constexpr auto sv_elas = StrView::make("elas");
    static constexpr auto sv_inelas = StrView::make("inelas");
    static constexpr auto sv_sans = StrView::make("sans");
    static constexpr auto sv_comp = StrView::make("comp");
    static constexpr auto sv_bkgd = StrView::make("bkgd");
    static constexpr auto sv_none = StrView::make("none");
    static constexpr auto sv_0 = StrView::make("0");
    if ( par_str == sv_bragg ) {
      setVarByStr( data, VarId::coh_elas, val_str );
      continue;
    } else if ( par_str == sv_elas ) {
      setVarByStr( data, VarId::coh_elas, val_str );
      setVarByStr( data, VarId::incoh_elas, val_str );
      setVarByStr( data, VarId::sans, val_str );
      continue;
    } else if ( par_str == sv_comp ) {
      bool comp_coh_elas = false;
      bool comp_incoh_elas = false;
      bool comp_inelas = false;
      bool comp_sans = false;
      for ( auto& e : val_str.splitTrimmedNoEmpty(',') ) {
        if ( e == sv_elas ) {
          comp_coh_elas = true;
          comp_incoh_elas = true;
        } else if ( e == sv_bragg || e == sv_coh_elas ) {
          comp_coh_elas = true;
        } else if ( e == sv_incoh_elas ) {
          comp_incoh_elas = true;
        } else if ( e == sv_inelas ) {
          comp_inelas = true;
        } else if ( e == sv_sans ) {
          comp_sans = true;
        } else {
          NCRYSTAL_THROW2(BadInput,"Unsupported component requested via the \"comp\" cfg-variable: \""<<e<<"\"");
        }
      }
      //Finally, turn off all comp_xxx that are now false (don't *enable*
      //anything, or we might suddenly introduce components that have been
      //already turned off):
      if ( !comp_coh_elas )
        setVarByStr( data, VarId::coh_elas, sv_0 );
      if ( !comp_incoh_elas )
        setVarByStr( data, VarId::incoh_elas, sv_0 );
      if ( !comp_inelas )
        setVarByStr( data, VarId::inelas, sv_0 );
      if ( !comp_sans )
        setVarByStr( data, VarId::sans, sv_0 );
      continue;
    } else if ( par_str == sv_bkgd ) {
      if ( !isOneOf(val_str,sv_none,sv_0) )
        NCRYSTAL_THROW(BadInput,"The \"bkgd\" parameter is obsolete and is available for backwards compatibility "
                       "only with the values \"0\" or \"none\". For control of inelastic, incoherent-elastic, or SANS "
                       "scattering, one must now instead use the parameters \"incoh_elas\", \"inelas\", and \"sans\".");
      setVarByStr( data, VarId::incoh_elas, sv_0 );
      setVarByStr( data, VarId::inelas, sv_0 );
      setVarByStr( data, VarId::sans, sv_0 );
      continue;
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!! IMPORTANT                            !!!
    // !!! Any modification to pseudo or toplvl !!!
    // !!! parameters must be reflected in the  !!!
    // !!! docs in NCCfgVars.cc                 !!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //Top-level parameters:
    static constexpr auto sv_density = StrView::make("density");
    static constexpr auto sv_phasechoice = StrView::make("phasechoice");
    if ( par_str == sv_density ) {
      Optional<double> val;
      Optional<DensityState::Type> dt;
      if ( val_str.endswith( StrView::make("x") ) ) {
        val = val_str.substr(0,val_str.size()-1).toDbl();
        dt = DensityState::Type::SCALEFACTOR;
      } else if ( val_str.endswith( StrView::make("gcm3") ) ) {
        val = val_str.substr(0,val_str.size()-4).toDbl();
        dt = DensityState::Type::DENSITY;
      } else if ( val_str.endswith( StrView::make("kgm3") ) ) {
        val = val_str.substr(0,val_str.size()-4).toDbl();
        if ( val.has_value() )
          val.value() *= 0.001;
        dt = DensityState::Type::DENSITY;
      } else if ( val_str.endswith( StrView::make("perAa3") ) ) {
        val = val_str.substr(0,val_str.size()-6).toDbl();
        dt = DensityState::Type::NUMBERDENSITY;
      }
      if ( !val.has_value() || !dt.has_value() )
        NCRYSTAL_THROW2(BadInput,"Bad value syntax in \"density\" parameter: \""<<val_str
                        <<"\". Expects value with unit \"gcm3\", \"kgm3\", \"perAa3\", or \"x\"");
      toplvlVars.emplace_back(DensityState{dt.value(),val.value()});
      continue;
    } else if ( par_str == sv_phasechoice ) {
      auto val = val_str.toInt();
      if ( !val.has_value() )
        NCRYSTAL_THROW2(BadInput,"Invalid phase choice specification: "<<val_str);
      if ( val.value()<0 || val.value()>10000 )
        NCRYSTAL_THROW2(BadInput,"Invalid phase choice index (out of range): "<<val.value());
      toplvlVars.emplace_back( TopLvlVar::phasechoice_t(), static_cast<unsigned>( val.value() ) );
      continue;
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!! IMPORTANT                            !!!
    // !!! Any modification to pseudo or toplvl !!!
    // !!! parameters must be reflected in the  !!!
    // !!! docs in NCCfgVars.cc                 !!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    //Unknown parameter! Try to produce the most helpful error message:
    if ( par_str.empty() )
      NCRYSTAL_THROW2(BadInput,errintro<<"missing parameter name.");
    auto semicolonpos = par_str.find(';');
    if ( semicolonpos != StrView::npos ) {
      //This means we encountered an ';' char before the '=' char, e.g. in "foo;bar=val".
      auto pn = par_str.substr(0,semicolonpos).trimmed();
      if ( !pn.empty() ) {
        NCRYSTAL_THROW2(BadInput,errintro<<"missing parameter value for parameter \""<<pn<<"\".");
      } else {
        NCRYSTAL_THROW2(BadInput,errintro<<"syntax error around the part \""<<par_str<<"\".");
      }
    }
    static constexpr auto sv_packfact = StrView::make("packfact");
    if ( par_str == sv_packfact )
      NCRYSTAL_THROW2(BadInput,errintro<<"The \"packfact\" parameter is obsolete. Previous"
                      " \"packfact=0.8\" is now \"density=0.8x\" (but note that the density scaling"
                      " of multiple such statements is now commulative!). Alternatively, it is"
                      " now also possible to set an absolute density, e.g. \"density=2.3gcm3\".");
    NCRYSTAL_THROW2(BadInput,errintro<<"unknown parameter \""<<par_str<<"\".");
  }
  return toplvlVars;
}

void NC::Cfg::CfgManip::stream( const CfgData& data, std::ostream& os, const VarIdFilter& filter )
{
  bool not_first = false;
  for ( auto& e : data() ) {
    auto varid = e.metaData();
    if ( filter != nullptr && !filter(varid) )
      continue;
    if ( not_first )
      os << ';';
    not_first = true;
    os << varName(varid);
    os << '=';
    varInfo(varid).stream( os , e );
  }
}

NC::Cfg::VarIdFilter NC::Cfg::CfgManip::createFilter(const VarIdList& vl, FilterType fltype )
{
  VarIdFilter result;
  if ( vl.empty() ) {
    //Trivial:
    if ( fltype == FilterType::ExcludeListed )
      result = [](VarId){ return true; };
    else
      result = [](VarId){ return false; };
    return result;
  }
  //Less trivial:
  class Filter {
    VarIdList m_local;
    bool m_exclude;
  public:
    Filter(const VarIdList& vl, bool exclude)
      : m_local(SVAllowCopy,vl),
        m_exclude(exclude)
    {
      std::sort(m_local.begin(),m_local.end());
    }
    bool operator()( VarId varid ) const
    {
      auto it = std::lower_bound(m_local.begin(),m_local.end(),varid);
      const bool found = it != m_local.end() && *it == varid;
      return found != m_exclude;
    }
  };
  result = Filter(vl, fltype == FilterType::ExcludeListed );
  return result;
}

bool NC::Cfg::CfgManip::lessThan( const CfgData& a, const CfgData& b )
{
  if ( &a == &b )
    return false;

  //First we simply do the cheapest thing - compare which variables are set at
  //all. Remember the m_data arrays are always sorted by metaData (aka varId)
  //value.
  if ( a().size() != b().size() )
    return a().size() < b().size();

  for ( auto i : ncrange( a().size() ) )
    if ( a()[i].metaData() != b()[i].metaData() )
      return a()[i].metaData() < b()[i].metaData();

  //Ok, the two have the same variables set - now we must check each variable
  //separately.
  for ( auto i : ncrange( a().size() ) ) {
    auto varid = static_cast<VarId>(a()[i].metaData());
    nc_assert( varid == static_cast<VarId>(b()[i].metaData()) );
    auto& varinfo = varInfo( varid );
    auto cmpval = varinfo.bufCmp( a()[i], b()[i] );
    if ( cmpval != 0 )
      return cmpval < 0;
  }

  return false;
}

bool NC::Cfg::CfgManip::equal( const CfgData& a, const CfgData& b )
{
  if ( &a == &b )
    return true;
  //Like lessThan, cheapest checks first.
  if ( a().size() != b().size() )
    return false;
  for ( auto i : ncrange( a().size() ) )
    if ( a()[i].metaData() != b()[i].metaData() )
      return false;
  for ( auto i : ncrange( a().size() ) ) {
    auto varid = static_cast<VarId>(a()[i].metaData());
    nc_assert( varid == static_cast<VarId>(b()[i].metaData()) );
    auto& varinfo = varInfo( varid );
    auto cmpval = varinfo.bufCmp( a()[i], b()[i] );
    if ( cmpval != 0 )
      return false;
  }
  return true;
}

bool NC::Cfg::CfgManip::isSingleCrystal(const CfgData& data)
{
  for ( auto& e : data() ) {
    if ( isOneOf(static_cast<VarId>(e.metaData()),VarId::mos,VarId::dir1,VarId::dir2,VarId::dirtol) )
      return true;
  }
  return false;
}

NC::Cfg::VarIdList NC::Cfg::CfgManip::findCommonEntries( std::function<const CfgData*()> cfgDataIter )
{
  VarIdList result;
  auto cfg0 = cfgDataIter();
  if (!cfg0)
    return result;//empty range
  //Candidates are any variable in cfg0, if they also fit the rest.
  SmallVector<std::pair<VarId,const VarBuf*>,16> v;
  for ( auto& e : (*cfg0)() )
    v.emplace_back( e.metaData(), &e );
  //Test the rest:
  while ( auto cfg = cfgDataIter() ) {
    for ( auto& e : v ) {
      if ( !e.second )
        continue;//already disqualified
      auto buf = searchBuf( *cfg, e.first );
      if ( !buf || varInfo(e.first).bufCmp(*e.second,*buf) != 0 ) {
        e.second = nullptr;
        continue;
      }
    }
  }
  for ( auto& e : v )
    if ( e.second )
      result.push_back(e.first);
  return result;
}

void NC::Cfg::CfgManip::streamJSON( const CfgData& data, std::ostream& os )
{
  os << '[';
  unsigned i(0);
  for ( auto& e : data() ) {
    if ( i++ )
      os << ',';
    auto varinfo = varInfo( e.metaData() );
    os<<'[';
    ::NCrystal::streamJSON(os,varinfo.nameSV());
    os<<',';
    varinfo.streamAsJSON(os,e);
    os<<']';
  }
  os << ']';
}


std::vector<NC::VectS> NC::Cfg::CfgManip::get_atomdb_parsed(const CfgData& data)
{
  auto sv_atomdb = get_atomdb(data);
  nc_assert(sv_atomdb.has_value());
  std::vector<NC::VectS> v;
  if ( sv_atomdb.empty() )
    return v;//usual case
  auto lines = sv_atomdb.splitTrimmedNoEmpty('@');
  v.reserve(lines.size());
  for ( const auto& line : lines ) {
    auto parts = line.split_any<8,StrView::SplitKeepEmpty::No,StrView::SplitTrimParts::Yes>(" \t\n\r:");
    if ( parts.empty() )
      continue;
    v.emplace_back();
    v.back().reserve(parts.size());
    for ( auto& e : parts )
      v.back().push_back(e.to_string());
  }
  v.shrink_to_fit();
  return v;
}

void NC::Cfg::CfgManip::checkParamConsistency_Info( const CfgData& data )
{
  auto buf_dcutoff = searchBuf( data, VarId::dcutoff );
  auto buf_dcutoffup = searchBuf( data, VarId::dcutoffup );
  if ( buf_dcutoff || buf_dcutoffup ) {
    auto parval_dcutoff = getValueFromBufPtr<vardef_dcutoff>(buf_dcutoff);
    auto parval_dcutoffup = getValueFromBufPtr<vardef_dcutoffup>(buf_dcutoffup);
    if (parval_dcutoff>=parval_dcutoffup)
      NCRYSTAL_THROW(BadInput,"dcutoff must be less than dcutoffup");
  }
}

void NC::Cfg::CfgManip::checkParamConsistency_ScatterBase( const CfgData& )
{
  //Currently all ScatterBase parameters can be set independently without
  //consistency constraints.
}

#include "NCrystal/internal/utils/NCLatticeUtils.hh"

void NC::Cfg::CfgManip::checkParamConsistency_ScatterExtra( const CfgData& data )
{
  auto buf_mos = searchBuf( data, VarId::mos );
  auto buf_dir1 = searchBuf( data, VarId::dir1 );
  auto buf_dir2 = searchBuf( data, VarId::dir2 );
  auto buf_dirtol = searchBuf( data, VarId::dirtol );

  int nOrient = (buf_dir1?1:0) + (buf_dir2?1:0) + (buf_mos?1:0);
  if (nOrient!=0 && nOrient<3)
    NCRYSTAL_THROW(BadInput,"Must set all or none of mos, dir1 and dir2 parameters");

  if ( !nOrient ) {
    if ( buf_dirtol )
      NCRYSTAL_THROW(BadInput,"mos, dir1 and dir2 parameters must all be set when dirtol is set");
    return;
  }

  auto parval_dir1 = getValueFromBufPtr<vardef_dir1>(buf_dir1);
  auto parval_dir2 = getValueFromBufPtr<vardef_dir2>(buf_dir2);
  auto parval_dirtol = getValueFromBufPtr<vardef_dirtol>(buf_dirtol);

  precheckLatticeOrientDef( parval_dir1, parval_dir2, parval_dirtol );
}

void NC::Cfg::CfgManip::checkParamConsistency_Absorption( const CfgData& )
{
  //Currently all Absorption parameters can be set independently without
  //consistency constraints.
}
