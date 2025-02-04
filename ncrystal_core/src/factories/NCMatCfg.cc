
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

#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/interfaces/NCSCOrientation.hh"
#include "NCrystal/internal/cfgutils/NCCfgManip.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <sstream>
namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    using Cfg::CfgManip;

    using DensityState = DensityState;
    Optional<DensityState> accumulateDensityState( const Optional<DensityState>& current, const DensityState& ds ) {
      ds.validate();
      if ( ds.type == DensityState::Type::SCALEFACTOR && ds.value == 1.0 )
        return current;//no change
      if ( !current.has_value() || ds.type != DensityState::Type::SCALEFACTOR )
        return ds;//override existing
      //scale existing, keeping type (not just a single return statement here due to c++11 compat.):
      DensityState newds;
      newds.type = current.value().type;
      newds.value = current.value().value * ds.value;
      newds.validate();
      return newds;
    }

    std::pair<StrView,StrView> splitCfgOnFilename( const StrView& sv )
    {
      //Essentially splits once on ';' and returns trimmed parts. If ';' is not
      //found, the second returned part is (i.e. only a filename was
      //provided. Also checks for the obsolete ignorefilecfg keyword and throws
      //a clear error if it is used.
      auto i = sv.find(';');
      if ( i == StrView::npos )
        return { sv.trimmed(), StrView::make("") };
      auto sv_fn = sv.substr(0,i).trimmed();
      auto sv_cfg = sv.substr(i+1).trimmed();
      if ( sv_cfg.startswith( StrView::make("ignorefilecfg") ) )
        NCRYSTAL_THROW2(BadInput,"Usage of the \"ignorefilecfg\" keyword in cfg strings is obsolete and no longer supported.");
      return { sv_fn, sv_cfg };
    }
  }//end anon namespace

  struct MatCfg::constructor_args : private MoveOnly {
    struct MultiPhase {
      Cfg::TopLvlVarList toplvlvars;
      PhaseList phases;
    };
    struct SinglePhase {
      TextDataSP td;
      StrView paramstr;
      StrView dataname;//Optional, for better reconstruction of original string
                       //(although the TextDataSP of course has a fall-back
                       //option, it might include factory name even though the
                       //original did not).
    };
    Variant<MultiPhase,SinglePhase> cfg;
  };
}

class NC::MatCfg::Impl2 {
  //Variables which are routinely discarded by the factory infrastructure are
  //kept on a separate COWPimpl object.
public:
  PhaseChoices m_phaseChoices;
  Optional<DensityState> m_densityState;

  static void checkPhaseChoiceRange( unsigned i )
  {
    //Just to keep it reasonable and catch obvious errors:
    if ( i > Cfg::TopLvlVar::phasechoice_max )
      NCRYSTAL_THROW2(BadInput,"Invalid phase choice index (too high): "<<i);
  }

  static void apply( const Cfg::TopLvlVarList& toplvlpars, COWPimpl<Impl2>& impl2, COWPimpl<Impl2>::Modifier * modptr = nullptr )
  {
    //Collect:
    MatCfg::PhaseChoices collected_phaseChoiceList;
    Optional<DensityState> collected_densityState;
    for ( auto& e : toplvlpars ) {
      if ( e.isPhaseChoice() ) {
        checkPhaseChoiceRange( e.getPhaseChoice() );
        collected_phaseChoiceList.push_back( e.getPhaseChoice() );
      } else {
        nc_assert( e.isDensity() );
        collected_densityState = accumulateDensityState( collected_densityState, e.getDensity() );
      }
    }
    //apply:
    if ( collected_phaseChoiceList.empty() && !collected_densityState.has_value() )
      return;
    Optional<COWPimpl<Impl2>::Modifier> modholder;
    if (!modptr) {
      modholder = impl2.modify();
      modptr = &modholder.value();
    }
    auto& mod = *modptr;
    if ( collected_densityState.has_value() )
      mod->m_densityState.set( accumulateDensityState( mod->m_densityState, collected_densityState.value() ) );

    for ( auto i : collected_phaseChoiceList )
      mod->m_phaseChoices.push_back(i);
  }
};

class NC::MatCfg::Impl {
public:
  TextDataUID m_textDataUID;
  std::string m_textDataType;
  DataSourceName m_dataSourceName;//data name as passed to MatCfg constructor or an appropriate "<anonymousXXX-data>"
  std::shared_ptr<PhaseList> m_phases;//set if multiphase
  Cfg::CfgData m_cfgdata;//Variables (if single phase)

  static PhaseList clonePhaseList( const PhaseList& pl )
  {
    PhaseList out;
    out.reserve( pl.size() );
    for ( auto& ph : pl )
      out.emplace_back( ph.first, ph.second.clone() );
    return out;
  }

  static std::string extractEmbeddedCfgStr( const DataSourceName&, const TextData& );
  static PhaseList cleanupAndCheckPhases( PhaseList&&, unsigned& );
  bool isMultiPhase() const noexcept { return m_phases != nullptr; }

  template<class TVar, class TSetVarMemberFunction>
  void setVar( const TVar& val, TSetVarMemberFunction setfct )
  {
    if ( m_phases==nullptr ) {
      //Single phase:
      setfct(m_cfgdata,val);
    } else {
      //Multiphase. To only parse once and ensure VarBuf remote buffer sharing,
      //we do it like this:
      Cfg::CfgData tmp;
      setfct( tmp, val );
      for ( auto& ph : *m_phases )
        CfgManip::apply( ph.second.m_impl.modify()->m_cfgdata, tmp );
    }
  }

  const Cfg::CfgData* tryReadVar( Cfg::VarId varid ) const
  {
    if (!isMultiPhase())
      return &m_cfgdata;
    const PhaseList& pl = *m_phases;
#ifndef NDEBUG
    for ( auto& ph : pl )
      nc_assert(!ph.second.isMultiPhase());
#endif
    nc_assert( pl.size()>=2 );
    const auto& cfg0 = pl.front().second.m_impl->m_cfgdata;
    auto phaselist_other = Span<const Phase>{pl}.subspan(1);
    if ( !CfgManip::hasValueSet(cfg0,varid) ) {
      //OK if all unset (optimise for this, the usual case):
      for ( auto& ph : phaselist_other )
        if ( CfgManip::hasValueSet( ph.second.m_impl->m_cfgdata,varid) )
          return nullptr;
      return &cfg0;
    } else {
      //OK if all set to same value
      for ( auto& ph : phaselist_other )
        if ( !CfgManip::hasSameValue( cfg0, ph.second.m_impl->m_cfgdata, varid ) )
          return nullptr;
      return &cfg0;
    }
  }

  const Cfg::CfgData& readVar( Cfg::VarId varid ) const
  {
    if (!isMultiPhase())
      return m_cfgdata;
    auto cfgptr = tryReadVar(varid);
    if (!cfgptr)
      NCRYSTAL_THROW2(CalcError,"Could not determine unique value of parameter \""<<Cfg::varName(varid)
                      <<"\" on multiphase MatCfg object (different values found in different phases).");
    return *cfgptr;
  }

  static Optional<constructor_args> decodeAndInitMultiPhaseCfg( StrView input )
  {
    //Decode "phases<FRAC1*CFG1&..&FRACN*CFGN>[;COMMONCFG]"
    Optional<constructor_args> result;
    input = input.trimmed();//redundant, but just to be safe
    auto badSyntax = [&input]()
    {
      NCRYSTAL_THROW2(BadInput,"Invalid syntax in multiphase configuration "
                      "string: \""<<input<<"\"");
    };

    ///////////////////////////////////////////////////////////////
    //Extract phases (frac1*cfg1&...&fracn*cfgn) and commoncfg parts:

    StrView commoncfgstr;
    StrView phases_str;
    {
      auto tmp = input.splitTrimmed<2>('<');
      if ( tmp.size() < 2 ) {
        //Had no '<'.
        if ( !input.contains_any("*>&") ) {
          //String started with "phases" but had no "<>*&" chars. This is OK,
          //might just be a single phase config string like "phases.ncmat".
          return result;
        }
        badSyntax();
      }
      if ( tmp.size() != 2 || tmp.front()!="phases" )
        badSyntax();
      auto tmp2 = tmp.back().splitTrimmed('>');
      if ( tmp2.size()!=2 )
        badSyntax();
      phases_str = tmp2.front();
      commoncfgstr = tmp2.back().trimmed();
    }

    if ( !commoncfgstr.empty() && !commoncfgstr.startswith(';') )
      NCRYSTAL_THROW2(BadInput,"Invalid syntax in multiphase configuration "
                      "string: \""<<input<<"\" (the part after the \">\""
                      " character must begin with a \";\" character).");

    //Create result like this to (hopefully) avoid redundant copies:
    result.emplace();//create value
    result.value().cfg = constructor_args::MultiPhase();
    auto& res_mp = result.value().cfg.get<constructor_args::MultiPhase>();

    ///////////////////////////////////////////////////////////////
    //Parse common cfg (top lvl vars stay at top level, the rest will be applied
    //to individual phases below):
    Cfg::CfgData common_cfgdata;
    res_mp.toplvlvars = CfgManip::applyStrCfg( common_cfgdata, commoncfgstr );

    ///////////////////////////////////////////////////////////////
    //Individual phase components:
    auto& phaselist = res_mp.phases;

    StableSum fracsum;
    for ( auto&& phasestr : phases_str.splitTrimmed('&') ) {
      auto tmp = phasestr.splitTrimmed('*');
      if ( tmp.size() != 2 ) {
        //Nice error if user used "&&" instead of "&":
        if ( phases_str.contains(StrView::make("&&")) ) {
          NCRYSTAL_THROW2(BadInput,"Invalid delimiter (\"&&\" instead of the"
                          " correct \"&\") in multiphase configuration "
                          "string: \""<<input<<"\"");
        } else {
          badSyntax();
        }
      }

      {
        auto badchar_strrep = findForbiddenChar( tmp.front(), Cfg::forbidden_chars_non_multiphase,ExtraForbidOpt::RequireSimpleASCII );
        if ( badchar_strrep.has_value() )
          NCRYSTAL_THROW2(BadInput,"Forbidden character "<<badchar_strrep.value()
                          <<" in configuration string! Problem found in string: "<<input);
      }
      double fraction;
      if ( !safe_str2dbl(tmp.front(),fraction) )
        fraction = -1.0;
      if ( fraction <= 0.0 || fraction > 1.0 )
        NCRYSTAL_THROW2(BadInput,"Invalid value of component fraction specified in "
                        "multiphase configuration string: \""<<tmp.front()<<"\"");
      fracsum.add(fraction);
      phaselist.emplace_back( Phase{ fraction, MatCfg(tmp.back().to_string()) } );
      nc_assert_always( !phaselist.back().second.isMultiPhase() );

      if ( !CfgManip::empty(common_cfgdata) ) {
        //Transfer common parameters to each daughter object (Note that this
        //does NOT transfer the topvll NB: this does on purpose NOT transfer the
        //common_toplvlvars parameters:
        CfgManip::apply( phaselist.back().second.m_impl.modify()->m_cfgdata,  common_cfgdata );
      }
    }

    if ( phaselist.empty() )
      badSyntax();

    ///////////////////////////////////////////////////////////////
    //Check/fix total fraction:
    const double fractot = fracsum.sum();
    if ( ! ( ncabs( fractot - 1.0 ) <= 1e-6 ) )
      NCRYSTAL_THROW2( BadInput, "Component fractions in multiphase configuration"
                       " string do not add up to unity: \""<<input<<"\"" );
    if ( fractot != 1.0 ) {
      //tiny correction, snap to unity (although redundant, will be done later as well):
      for ( auto& ph : phaselist )
        ph.first /= fractot;
    }

    return result;
  }

  void dump( const MatCfg *, std::ostream& out, bool add_endl, bool includePhaseChoice ) const;
  bool compareIgnoringTextDataUID(const MatCfg& ) const;
  std::string toStrCfg( const MatCfg&, bool include_datafile, const Cfg::VarIdFilter& only_pars, bool includePhaseChoice ) const;

};

std::string NC::MatCfg::toEmbeddableCfg() const
{
  if ( isMultiPhase()  )
    NCRYSTAL_THROW(BadInput,"MatCfg::toEmbeddableCfg() can not be"
                   " called for multiphase configurations");
  if ( m_impl2->m_densityState.has_value()
       && m_impl2->m_densityState.value().type == DensityState::Type::SCALEFACTOR
       && m_impl2->m_densityState.value().value != 1.0 ) {
    NCRYSTAL_THROW(BadInput,"MatCfg::toEmbeddableCfg() can not be"
                   " called with configurations where the density state is a scale factor");
  }
  std::stringstream out;
  out << "NCRYSTALMATCFG[" << m_impl->toStrCfg(*this,false,nullptr,true) << ']';
  return out.str();
}

std::string NC::MatCfg::Impl::toStrCfg( const MatCfg& cfgobj,
                                        bool include_datafile,
                                        const Cfg::VarIdFilter& only_pars,
                                        bool includePhaseChoice ) const
{
  std::ostringstream out;
  auto addSeparatorIfNeeded = [&out](){ if (!out.str().empty()) out<<';'; };

  auto appendPhaseChoices = [includePhaseChoice,&addSeparatorIfNeeded](std::ostringstream& os, const MatCfg& cfg)
  {
    const auto & pc = cfg.getPhaseChoices();
    if ( !includePhaseChoice || pc.empty() )
      return;
    bool first = true;
    for ( auto i : pc ) {
      if ( first ) {
        first = false;
        addSeparatorIfNeeded();
      } else {
        os << ';';
      }
      os << "phasechoice=" << std::to_string(i);
    }
  };
  auto appendDensity = [&addSeparatorIfNeeded](std::ostringstream&os, const MatCfg& cfg)
  {
    if ( !cfg.m_impl2->m_densityState.has_value() )
      return;
    auto& ds = cfg.m_impl2->m_densityState.value();
    if ( ds.type == DensityState::Type::SCALEFACTOR && ds.value == 1.0 )
      return;
    addSeparatorIfNeeded();
    os << "density=" << ds;
  };

  if ( isMultiPhase() ) {
    const PhaseList& phlist = *m_phases;
    class CfgDataIter {
      PhaseList::const_iterator m_it,m_itE;
    public:
      CfgDataIter(const PhaseList& pl) : m_it(pl.begin()), m_itE(pl.end()) {}
      const Cfg::CfgData* operator()() { return m_it == m_itE ? nullptr : &(m_it++)->second.m_impl->m_cfgdata; }
    };
    auto common_entries = CfgManip::findCommonEntries( CfgDataIter(phlist) );
    if ( only_pars ) {
      //Remove entries not in only_pars from common_entries:
      Cfg::VarIdList tmp;
      std::copy_if(common_entries.begin(), common_entries.end(), std::back_inserter(tmp), only_pars );
      std::swap(common_entries,tmp);
    }

    auto varfilter_onlycommon = CfgManip::createFilter( common_entries, CfgManip::FilterType::OnlyListed );
    auto varfilter_nocommon = CfgManip::createFilter( common_entries, CfgManip::FilterType::ExcludeListed );
    if ( only_pars != nullptr ) {
      class Filter_A_and_B {
        Cfg::VarIdFilter m_a, m_b;
      public:
        Filter_A_and_B( Cfg::VarIdFilter a, Cfg::VarIdFilter b ) : m_a(std::move(a)), m_b(std::move(b)) {}
        bool operator()(Cfg::VarId v) const { return ( !m_a || m_a(v) ) && ( !m_b || m_b(v) ); }
      };
      varfilter_onlycommon = Filter_A_and_B( only_pars, std::move(varfilter_onlycommon) );
      varfilter_nocommon = Filter_A_and_B( only_pars, std::move(varfilter_nocommon) );
    }
    nc_assert_always(m_phases!=nullptr);
    nc_assert_always(include_datafile=true);
    out<<"phases<";
    nc_assert(phlist.size()>=2);

    bool first=true;
    for ( const auto& ph : *m_phases ) {
      nc_assert( ph.second.isSinglePhase() );
      if (!first)
        out << '&';
      first = false;
      out << fmt(ph.first)
          << '*'
          << ph.second.m_impl->toStrCfg(ph.second,true,varfilter_nocommon,includePhaseChoice);
    }

    out << '>';
    if (!common_entries.empty()) {
      const Cfg::CfgData& phase0_cfg_data = phlist.front().second.m_impl->m_cfgdata;
      if ( !CfgManip::empty(phase0_cfg_data,varfilter_onlycommon) ) {
        addSeparatorIfNeeded();
        CfgManip::stream(phase0_cfg_data,out,varfilter_onlycommon);
      }
    }

    appendDensity(out,cfgobj);
    appendPhaseChoices(out,cfgobj);
    return out.str();
  }

  //single phase:
  if ( include_datafile )
    out << m_dataSourceName.str();

  if ( !CfgManip::empty(m_cfgdata,only_pars) ) {
    addSeparatorIfNeeded();
    CfgManip::stream(m_cfgdata,out,only_pars);
  }

  appendDensity(out,cfgobj);
  appendPhaseChoices(out,cfgobj);
  return out.str();
}

std::string NC::MatCfg::toStrCfg( bool include_datafile ) const
{
  if ( isMultiPhase() && !include_datafile )
    NCRYSTAL_THROW(BadInput,"MatCfg::toStrCfg can not be called"
                   " with include_datafile=false for multiphase"
                   " configurations");
  return m_impl->toStrCfg( *this, include_datafile, nullptr, true );
}

bool NC::MatCfg::isSingleCrystal() const
{
  if ( isMultiPhase() )
    NCRYSTAL_THROW(CalcError,"MatCfg::isSingleCrystal() should not be called for multiphase materials");
  return CfgManip::isSingleCrystal( m_impl->m_cfgdata );
}

bool NC::MatCfg::isLayeredCrystal() const
{
  if ( isMultiPhase() )
    NCRYSTAL_THROW(CalcError,"MatCfg::isLayeredCrystal() should not be called for multiphase materials");
  return CfgManip::isLayeredCrystal( m_impl->m_cfgdata );
}

void NC::MatCfg::checkConsistency() const
{
  if ( m_impl2->m_densityState.has_value() )
    m_impl2->m_densityState.value().validate();//redundant, done when setting
  if ( m_impl->m_phases!=nullptr ) {
    for ( auto& ph : *(m_impl->m_phases) )
      ph.second.checkConsistency();
    return;
  }

  nc_assert(isSinglePhase());
  CfgManip::checkParamConsistency_Info( m_impl->m_cfgdata );
  CfgManip::checkParamConsistency_ScatterBase( m_impl->m_cfgdata );
  CfgManip::checkParamConsistency_ScatterExtra( m_impl->m_cfgdata );
  CfgManip::checkParamConsistency_Absorption( m_impl->m_cfgdata );
}

void NC::MatCfg::set_dir1( const OrientDir& od )
{
  m_impl.modify()->setVar( od, &CfgManip::set_dir1 );
}

void NC::MatCfg::set_dir2( const OrientDir& od )
{
  m_impl.modify()->setVar( od, &CfgManip::set_dir2 );
}

void NC::MatCfg::set_dir1( const HKLPoint& c, const LabAxis& l )
{
  m_impl.modify()->setVar( OrientDir{ c, l }, &CfgManip::set_dir1 );
}

void NC::MatCfg::set_dir1( const CrystalAxis& c, const LabAxis& l )
{
  m_impl.modify()->setVar( OrientDir{ c, l }, &CfgManip::set_dir1 );
}

void NC::MatCfg::set_dir2( const HKLPoint& c, const LabAxis& l )
{
  m_impl.modify()->setVar( OrientDir{c, l}, &CfgManip::set_dir2 );
}

void NC::MatCfg::set_dir2( const CrystalAxis& c, const LabAxis& l )
{
  m_impl.modify()->setVar( OrientDir{c, l}, &CfgManip::set_dir2 );
}

NC::OrientDir NC::MatCfg::get_dir1() const
{
  return CfgManip::get_dir1( m_impl->readVar(Cfg::VarId::dir1) );
}

NC::OrientDir NC::MatCfg::get_dir2() const
{
  return CfgManip::get_dir2( m_impl->readVar(Cfg::VarId::dir2) );
}

NC::SCOrientation NC::MatCfg::createSCOrientation() const
{
  if ( isMultiPhase() )
    NCRYSTAL_THROW(CalcError,"MatCfg::createSCOrientation() should not be called for multiphase materials");
  if ( !isSingleCrystal() )
    NCRYSTAL_THROW(MissingInfo,"Can only create SCOrientation object for single crystals (must set dir1, dir2, and mos parameters)");


  //Check that we are single phase or that none of the four relevant parameters
  //have different values in different phases:
  const auto& cfg_a = m_impl->readVar(Cfg::VarId::dir1);
  const auto& cfg_b = m_impl->readVar(Cfg::VarId::dir2);
  const auto& cfg_c = m_impl->readVar(Cfg::VarId::mos);
  const auto& cfg = m_impl->readVar(Cfg::VarId::mos);
  nc_assert_always( &cfg == &cfg_a && &cfg == &cfg_b && &cfg == &cfg_c );

  return CfgManip::createSCOrientation<SCOrientation>(cfg);
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    void doSetSCOrient( Cfg::CfgData& data, const SCOrientation::Data& orient )
    {
      CfgManip::set_dir1( data, orient.dir1 );
      CfgManip::set_dir2( data, orient.dir2 );
      CfgManip::set_dirtol( data, orient.dirtol );
    }
  }
}

void NC::MatCfg::setOrientation( const SCOrientation& sco )
{
  if (!sco.isComplete())
    NCRYSTAL_THROW(BadInput,"setOrientation called with incomplete SCOrientation object");
  m_impl.modify()->setVar( sco.getData(), &doSetSCOrient );
}

void NC::MatCfg::apply( const Cfg::CfgData& cfgData )
{
  if ( CfgManip::empty(cfgData) )
    return;
  auto mod = m_impl.modify();
  if ( isMultiPhase() ) {
    for ( auto& ph : *mod->m_phases )
      ph.second.apply(cfgData);
  } else {
    CfgManip::apply( mod->m_cfgdata, cfgData);
  }
}

bool NC::MatCfg::isTrivial() const
{
  return m_impl2->m_phaseChoices.empty() && m_impl->m_phases==nullptr && !hasDensityOverride();
}

void NC::MatCfg::applyStrCfg( const std::string& str )
{
  Cfg::CfgData cfgData;
  auto toplvlpars = CfgManip::applyStrCfg(cfgData,str);
  apply(cfgData);
  Impl2::apply(toplvlpars,m_impl2);
}

NC::MatCfg::MatCfg( const char* datafile_and_parameters )
  : MatCfg(std::string(datafile_and_parameters))
{
}

NC::MatCfg NC::MatCfg::createFromRawData( std::string&& data, std::string pars, std::string ext )
{
  return MatCfg( from_raw_t(), std::move(data), std::move(pars), std::move(ext) );
}

NC::MatCfg::MatCfg( PhaseList&& phases )
  : MatCfg( [&phases]() -> constructor_args
  {
    constructor_args args;
    args.cfg = constructor_args::MultiPhase{ {}, std::move(phases) };
    return args;
  }())
{
}

NC::MatCfg::MatCfg( const PhaseList& phases )
  : MatCfg( [&phases]() -> constructor_args
  {
    constructor_args args;
    args.cfg = constructor_args::MultiPhase{ {}, Impl::clonePhaseList(phases) };
    return args;
  }())
{
}

NC::MatCfg::MatCfg( const std::string& datafile_and_parameters )
  : MatCfg( [&datafile_and_parameters]() -> constructor_args
  {
    auto input = StrView(datafile_and_parameters);

    {
      auto badchar_strrep = findForbiddenChar( input, Cfg::forbidden_chars_multiphase );
      if ( badchar_strrep.has_value() )
        NCRYSTAL_THROW2(BadInput,"Forbidden character "<<badchar_strrep.value()
                        <<" in configuration string! Problem found in string: "<<input);
    }
    input = input.trimmed();

    //Check for multiphase syntax:
    auto pos_MPOnlyChar = input.find_first_of("<>&*");
    if ( pos_MPOnlyChar != StrView::npos ) {
      if (!input.startswith("phases"))
        NCRYSTAL_THROW2(BadInput,"Invalid syntax in cfg-string "
                        "(\""<<input.at(pos_MPOnlyChar)<< "\" char only "
                        <<"allowed in multi-phase cfgs): \""<<input<<"\"");

      auto opt_multiphase_constructor_args = Impl::decodeAndInitMultiPhaseCfg(input);
      if ( opt_multiphase_constructor_args.has_value() )
        return std::move(opt_multiphase_constructor_args.value());
    }
    //Standard single-phase cfg string.
    StrView dataname,cfgstr;
    std::tie(dataname,cfgstr) = splitCfgOnFilename( input );
    nc_assert( dataname.has_value() );
    if ( dataname.empty() )
      NCRYSTAL_THROW2(BadInput,"Missing data name in \""<<input<<'"');

    {
      auto badchar_strrep = findForbiddenChar( cfgstr, Cfg::forbidden_chars_non_multiphase,ExtraForbidOpt::RequireSimpleASCII );
      if ( badchar_strrep.has_value() )
        NCRYSTAL_THROW2(BadInput,"Forbidden character "<<badchar_strrep.value()
                        <<" in configuration string! Problem found in string: "<<input);
      auto badchar_strrep2 = findForbiddenChar( dataname, Cfg::forbidden_chars_non_multiphase );
      if ( badchar_strrep2.has_value() )
        NCRYSTAL_THROW2(BadInput,"Forbidden character "<<badchar_strrep2.value()
                        <<" in configuration string! Problem found in string: "<<input);
    }
    nc_assert( dataname.has_value() );
    nc_assert( !dataname.empty() );
    auto td = FactImpl::createTextData( dataname.to_string() );
    constructor_args args;
    args.cfg = constructor_args::SinglePhase{ std::move(td), cfgstr, dataname};
    return args;
  }())
{
}

NC::MatCfg::MatCfg( TextDataSP td, std::string pars )
  : MatCfg( [&td,&pars]() -> constructor_args
  {
    constructor_args args;
    args.cfg = constructor_args::SinglePhase{ std::move(td), pars, NullOpt };
    return args;
  }())
{
}

NC::MatCfg::MatCfg( from_raw_t, std::string&& data, std::string pars, std::string dataType )
  : MatCfg( [&data,&pars,&dataType]() -> constructor_args
  {
    RawStrData rawdata(std::move(data));
    if ( dataType.empty() )
      dataType = FactImpl::guessDataType(rawdata);
    if ( dataType.empty() )
      NCRYSTAL_THROW2(BadInput,"Can not determine format of anonymous data (must be specified explicitly in this case):");
    auto td = makeSO<const TextData>( std::move(rawdata), TextData::DataType{std::move(dataType)} );
    constructor_args args;
    args.cfg = constructor_args::SinglePhase{ std::move(td), pars, NullOpt };
    return args;
  }())
{
}

NC::MatCfg::PhaseList NC::MatCfg::Impl::cleanupAndCheckPhases( PhaseList&& phlist, unsigned& recursion_lvl )
{
  if (++recursion_lvl==10000)
    NCRYSTAL_THROW(BadInput,"Self-referencing (or insanely deep) MatCfg::PhaseList detected");
  PhaseList out;
  out.reserve(phlist.size());
  //Flatten to a list of single-phase phases (could trigger if user constructed
  //PhaseList and used OO MatCfg constructor):

  for ( auto& ph : phlist ) {
    if ( ph.second.isSinglePhase() ) {
      out.push_back( std::move(ph) );
    } else {
      const double fraction = ph.first;
      PhaseList pl2 = clonePhaseList(*ph.second.m_impl->m_phases);
      pl2 = cleanupAndCheckPhases(std::move(pl2),recursion_lvl);
      for ( auto&& ph2 : pl2 )
        out.push_back( Phase{ fraction * ph2.first, std::move(ph2.second) } );
    }
  }
#ifndef NDEBUG
  for ( auto& ph : out )
    nc_assert_always( ph.second.isSinglePhase() );
#endif

  //Merge duplicate entries:
  auto cfgsEqual = [](const MatCfg& c1, const MatCfg& c2)
  {
    return ( c1.m_impl->m_textDataUID == c2.m_impl->m_textDataUID
             && !(c1<c2) && !(c2<c1) );
  };

  {
    PhaseList tmp;
    std::swap(tmp,out);
    out.reserve( tmp.size() );
    for ( auto i : ncrange(tmp.size()) ) {
      double fr(tmp.at(i).first);
      MatCfg cfg(std::move(tmp.at(i).second));
      if ( fr == 0.0 )
        continue;//ignore rather than throw error for now

      //Check if there is already an entry with same cfg:
      for ( auto& out_entry : out ) {
        if ( cfgsEqual(out_entry.second,cfg) ) {
          //consume on prior item
          out_entry.first += fr;
          fr = 0.0;
          break;
        }
      }
      if ( fr > 0.0 )
        out.push_back( Phase{fr, std::move(cfg) } );
    }
  }

  //Sanity check and snap-to-unity-sum of fractions (all in (0,1], unit sum):
  StableSum fracsum;
  for ( auto & ph : out ) {
    double fraction = ph.first;
    if ( fraction <= 0.0 || fraction > 1.0 )
      NCRYSTAL_THROW2(BadInput,"Invalid value of multiphase component fraction: "<<fraction<<"\"");
    fracsum.add(fraction);
  }

  const double fractot = fracsum.sum();
  if ( ! ( ncabs( fractot - 1.0 ) <= 1e-6 ) )
    NCRYSTAL_THROW2( BadInput, "Multiphase component fractions do not add up to unity!" );
  if ( fractot != 1.0 ) {
    //tiny correction, snap to unity:
    for ( auto& ph : out )
      ph.first /= fractot;
  }

  out.shrink_to_fit();
  return out;
}

NC::MatCfg::MatCfg( constructor_args&& args )
{
  //we just constructed COWPimpl objects from scratch, no need to lock as no one
  //else can refer to them:
  auto mod = m_impl.modifyWithoutLocking();
  auto mod2 = m_impl2.modifyWithoutLocking();

  if (args.cfg.has_value<constructor_args::MultiPhase>()) {
    auto& mpcfg = args.cfg.get<constructor_args::MultiPhase>();
    //This is the top-level object in a multi-phase material. All individual
    //phases have already been configured, and all we need is to adopt the phase
    //list.

    //Sanity check, flatten, snap fractions:
    unsigned recursion_lvl(0);
    PhaseList phases = Impl::cleanupAndCheckPhases( std::move(mpcfg.phases), recursion_lvl );
    if ( phases.size() == 1 ) {
      nc_assert(ncabs(phases.at(0).first-1.0)<1e-6);

      //Special case, not actually a multiphase cfg! We must become identical to
      //the first and only phase (and then of course apply the toplvlvars
      //below).

      //For safety, first remove all refs to present m_impl/m_impl2 before
      //replacing (although we didn't acquire a lock).

      mod.reset();
      mod2.reset();

      //Now become the first phase (which of course might be shared with other objects):
      *this = phases.at(0).second;
      //Apply top level variables::
      Impl2::apply(mpcfg.toplvlvars,m_impl2);
    } else {
      //Simply adopt the phase list:
      nc_assert_always(phases.size()>=2);
      mod->m_phases = std::make_shared<PhaseList>(std::move(phases));
      //Apply top level variables::
      Impl2::apply(mpcfg.toplvlvars,m_impl2,&mod2);
    }
  } else {
    nc_assert(args.cfg.has_value<constructor_args::SinglePhase>());
    auto& spcfg = args.cfg.get<constructor_args::SinglePhase>();
    //Standard single phase object.

    nc_assert(spcfg.td!=nullptr);
    m_textDataSP = std::move(spcfg.td);
    const TextData& textData = *m_textDataSP;

    //Cache TextData meta data which should remain even if thinned:
    mod->m_textDataUID = textData.dataUID();
    mod->m_textDataType = textData.dataType();
    const auto& textDataType = mod->m_textDataType;

    //Original file name (if any) is needed for more consistent re-serialisation:
    if ( spcfg.dataname.has_value() ) {
      mod->m_dataSourceName = spcfg.dataname.to_string();
    } else {
      if ( textDataType.empty() || textDataType == "unknown" ) {
        static DataSourceName s_anon("<anonymous>");
        mod->m_dataSourceName = s_anon;
      } else if ( textDataType == "ncmat" ) {
        static DataSourceName s_anon_ncmat("<anonymous-ncmat-data>");
        mod->m_dataSourceName = s_anon_ncmat;
      } else {
        std::ostringstream sss;
        sss << "<anonymous-"<< textDataType <<"-data>";
        mod->m_dataSourceName = sss.str();
      }
    }
    //Now look for embedded cfg str:
    std::string embedded_cfgstr = trim2(Impl::extractEmbeddedCfgStr(m_impl->m_dataSourceName,textData));
    if ( !embedded_cfgstr.empty() ) {
      auto embedded_toplvlvars = CfgManip::applyStrCfg(mod->m_cfgdata, embedded_cfgstr );
      Impl2::apply(embedded_toplvlvars,m_impl2,&mod2);
      //It is not allowed to have any variables with a commulative (as opposed
      //to overriding) effect in the embedded_cfgstr. Because if there would be
      //such variables, then reserialising the object as a cfg str would yield a
      //cfg string with those variables applied twice!
      //
      //Specifically it currently means we must disallow scaling densities and
      //phasechoice values.
      if ( !m_impl2->m_phaseChoices.empty() )
        NCRYSTAL_THROW2(BadInput,"phasechoice parameters are not allowed in embedded"
                        " cfg strings. Seen in data: "<<m_impl->m_dataSourceName);
      if  ( m_impl2->m_densityState.has_value()
            && m_impl2->m_densityState.value().type == DensityState::Type::SCALEFACTOR )
        NCRYSTAL_THROW2(BadInput,"density parameters with scale factors are not allowed"
                        " in embedded cfg strings. Seen in data: "<<m_impl->m_dataSourceName);
    }

    //Parse and apply cfg str:
    auto toplvlvars = CfgManip::applyStrCfg( mod->m_cfgdata, spcfg.paramstr );
    Impl2::apply(toplvlvars,m_impl2,&mod2);
  }
}

const NC::Cfg::CfgData& NC::MatCfg::rawCfgData() const
{
  if ( m_impl->isMultiPhase() )
    NCRYSTAL_THROW(LogicError,"MatCfg::rawCfgData called for multiphase object");
  return m_impl->m_cfgdata;
}

std::string NC::MatCfg::Impl::extractEmbeddedCfgStr( const DataSourceName& dataSourceName, const TextData& input )
{
  std::string res;
  std::string pattern="NCRYSTALMATCFG";
  for ( const std::string& line : input ) {
    std::size_t pos = line.find(pattern);
    if ( pos == std::string::npos )
      continue;
    if (!res.empty())
      NCRYSTAL_THROW2(BadInput,"Input data contains more than one "<<pattern<<" specification: "<<dataSourceName.str());
    std::string s = line.substr(pos+pattern.size());
    if (s.empty()||s.at(0)!='[')
      NCRYSTAL_THROW2(BadInput,"Input data contains "<<pattern<<" which is not followed by a '[' character: "<<dataSourceName.str());
    if (s.find(pattern)!=std::string::npos)
      NCRYSTAL_THROW2(BadInput,"Input data contains more than one "<<pattern<<" specification on a single line: "<<dataSourceName.str());
    s = s.substr(1);
    pos = s.find(']');
    if ( pos == std::string::npos )
      NCRYSTAL_THROW2(BadInput,"Input data contains "<<pattern<<" without a closing ']' character: "<<dataSourceName.str());
    res = s.substr(0,pos);
    if (res.empty())
      res = " ";//for detection of multiple occurances
  }
  trim(res);
  return res;
}

//All this nice MT-safe shallow copying stuff is generated by the COWPimpl and shared_ptr classes:
NC::MatCfg::MatCfg(const MatCfg&) = default;
NC::MatCfg& NC::MatCfg::operator=(const MatCfg&) = default;
NC::MatCfg::MatCfg( MatCfg&& ) = default;
NC::MatCfg& NC::MatCfg::operator=(MatCfg&&) = default;
NC::MatCfg::~MatCfg() = default;

NC::Temperature NC::MatCfg::get_temp() const { return CfgManip::get_temp( m_impl->readVar(Cfg::VarId::temp) ); }
double NC::MatCfg::get_dcutoff() const { return CfgManip::get_dcutoff( m_impl->readVar(Cfg::VarId::dcutoff) ); }
double NC::MatCfg::get_dcutoffup() const { return CfgManip::get_dcutoffup( m_impl->readVar(Cfg::VarId::dcutoffup) ); }
NC::MosaicityFWHM NC::MatCfg::get_mos() const { return CfgManip::get_mos( m_impl->readVar(Cfg::VarId::mos) ); }
double NC::MatCfg::get_mosprec() const { return CfgManip::get_mosprec( m_impl->readVar(Cfg::VarId::mosprec) ); }
double NC::MatCfg::get_sccutoff() const { return CfgManip::get_sccutoff( m_impl->readVar(Cfg::VarId::sccutoff) ); }
double NC::MatCfg::get_dirtol() const { return CfgManip::get_dirtol( m_impl->readVar(Cfg::VarId::dirtol) ); }
bool NC::MatCfg::get_coh_elas() const { return CfgManip::get_coh_elas( m_impl->readVar(Cfg::VarId::coh_elas) ); }
bool NC::MatCfg::get_incoh_elas() const { return CfgManip::get_incoh_elas( m_impl->readVar(Cfg::VarId::incoh_elas) ); }
bool NC::MatCfg::get_sans() const { return CfgManip::get_sans( m_impl->readVar(Cfg::VarId::sans) ); }

std::string NC::MatCfg::get_inelas() const { return CfgManip::get_inelas( m_impl->readVar(Cfg::VarId::inelas) ).to_string(); }
std::string NC::MatCfg::get_infofactory() const { return CfgManip::get_infofactory( m_impl->readVar(Cfg::VarId::infofactory) ).to_string(); }
std::string NC::MatCfg::get_scatfactory() const { return CfgManip::get_scatfactory( m_impl->readVar(Cfg::VarId::scatfactory) ).to_string(); }
std::string NC::MatCfg::get_absnfactory() const { return CfgManip::get_absnfactory( m_impl->readVar(Cfg::VarId::absnfactory) ).to_string(); }
const NC::LCAxis& NC::MatCfg::get_lcaxis() const { return CfgManip::get_lcaxis( m_impl->readVar(Cfg::VarId::lcaxis) ); }

void NC::MatCfg::set_temp( Temperature v ) { m_impl.modify()->setVar( v, &CfgManip::set_temp ); }
void NC::MatCfg::set_dcutoff( double v ) { m_impl.modify()->setVar( v, &CfgManip::set_dcutoff ); }
void NC::MatCfg::set_dcutoffup( double v ) { m_impl.modify()->setVar( v, &CfgManip::set_dcutoffup ); }
void NC::MatCfg::set_mos( MosaicityFWHM v ) { m_impl.modify()->setVar( v, &CfgManip::set_mos ); }
void NC::MatCfg::set_mosprec( double v ) { m_impl.modify()->setVar( v, &CfgManip::set_mosprec ); }
void NC::MatCfg::set_sccutoff( double v ) { m_impl.modify()->setVar( v, &CfgManip::set_sccutoff ); }
void NC::MatCfg::set_dirtol( double v ) { m_impl.modify()->setVar( v, &CfgManip::set_dirtol ); }
void NC::MatCfg::set_coh_elas( bool v ) { m_impl.modify()->setVar( v, &CfgManip::set_coh_elas ); }
void NC::MatCfg::set_incoh_elas( bool v ) { m_impl.modify()->setVar( v, &CfgManip::set_incoh_elas ); }
void NC::MatCfg::set_sans( bool v ) { m_impl.modify()->setVar( v, &CfgManip::set_sans ); }
void NC::MatCfg::set_inelas( const std::string& v ) { m_impl.modify()->setVar( v, &CfgManip::set_inelas_stdstr ); }
void NC::MatCfg::set_infofactory( const std::string& v ) { m_impl.modify()->setVar( v, &CfgManip::set_infofactory_stdstr ); }
void NC::MatCfg::set_scatfactory( const std::string& v ) { m_impl.modify()->setVar( v, &CfgManip::set_scatfactory_stdstr ); }
void NC::MatCfg::set_absnfactory( const std::string& v ) { m_impl.modify()->setVar( v, &CfgManip::set_absnfactory_stdstr ); }
void NC::MatCfg::set_lcmode( std::int_least32_t v ) { m_impl.modify()->setVar( v, &CfgManip::set_lcmode ); }
void NC::MatCfg::set_vdoslux( int v ) { m_impl.modify()->setVar( v, &CfgManip::set_vdoslux ); }
void NC::MatCfg::set_lcaxis( const LCAxis& axis ) { m_impl.modify()->setVar( axis, &CfgManip::set_lcaxis ); }
void NC::MatCfg::set_atomdb( const std::string& v ) { m_impl.modify()->setVar( v, &CfgManip::set_atomdb_stdstr ); }
std::int_least32_t NC::MatCfg::get_lcmode() const { return CfgManip::get_lcmode( m_impl->readVar(Cfg::VarId::lcmode) ); }
int NC::MatCfg::get_vdoslux() const { return CfgManip::get_vdoslux( m_impl->readVar(Cfg::VarId::vdoslux) ); }
std::string NC::MatCfg::get_atomdb() const { return CfgManip::get_atomdb( m_impl->readVar(Cfg::VarId::atomdb) ).to_string(); }
std::vector<NC::VectS> NC::MatCfg::get_atomdb_parsed() const { return CfgManip::get_atomdb_parsed( m_impl->readVar(Cfg::VarId::atomdb) ); }

void NC::MatCfg::Impl::dump( const MatCfg * self,
                             std::ostream& out,
                             bool add_endl,
                             bool includePhaseChoice ) const
{

  out << "MatCfg(\"";
  if ( isMultiPhase() ) {
    out << this->toStrCfg(*self,true,nullptr,includePhaseChoice);
  } else  {
    std::string strcfg = self->m_impl->toStrCfg( *self, false, nullptr,includePhaseChoice );
    out << m_dataSourceName.str();
    if (!strcfg.empty())
      out << (strcfg[0]==';'?"":";") << strcfg;
  }
  out<<"\")";
  if (add_endl)
    out<<std::endl;
}

void NC::MatCfg::dump( std::ostream& out, bool add_endl ) const
{
  m_impl->dump(this,out,add_endl,true);
}

const NC::TextDataUID NC::MatCfg::textDataUID() const
{
  if ( m_impl->isMultiPhase() )
    NCRYSTAL_THROW(LogicError,"MatCfg::textDataUID called for multiphase object");
  return m_impl->m_textDataUID;
}

const std::string& NC::MatCfg::getDataType() const
{
  if ( m_impl->isMultiPhase() )
    NCRYSTAL_THROW(LogicError,"MatCfg::getDataType called for multiphase object");
  return m_impl->m_textDataType;
}

bool NC::MatCfg::isThinned() const {
  if (!isMultiPhase())
    return m_textDataSP == nullptr;
  nc_assert(m_impl->m_phases!=nullptr);
  for ( const auto& ph : *m_impl->m_phases )
    if ( ! ph.second.isThinned() )
      return false;
  return true;
}

NC::MatCfg NC::MatCfg::cloneThinned() const
{
  MatCfg cfg(*this);
  if ( cfg.isThinned() )
    return cfg;
  cfg.m_textDataSP.reset();
  if ( cfg.m_impl->m_phases != nullptr ) {
    //must detach in order to thin phases.
    auto mod = cfg.m_impl.modify();
    nc_assert(mod->m_phases!=nullptr);
    for ( auto& ph : *mod->m_phases )
      ph.second = ph.second.cloneThinned();
  }
  return cfg;
}

const NC::TextData& NC::MatCfg::textData() const {
  //NB: Sync checks and messages between this fct and the next
  if ( m_impl->isMultiPhase() )
    NCRYSTAL_THROW(LogicError,"MatCfg::textDataSP called for multiphase object");
  if ( m_textDataSP == nullptr )
    NCRYSTAL_THROW(LogicError,"MatCfg::textDataSP/textData methods should not be"
                   " used in a MatCfg object which was thinned or moved-from.");
  return *m_textDataSP;
}

NC::TextDataSP NC::MatCfg::textDataSP() const
{
  //NB: Sync checks and messages between this fct and the previous
  if ( m_impl->isMultiPhase() )
    NCRYSTAL_THROW(LogicError,"MatCfg::textDataSP called for multiphase object");
  if ( m_textDataSP == nullptr )
    NCRYSTAL_THROW(LogicError,"MatCfg::textDataSP/textData methods should not be"
                   " used in a MatCfg object which was thinned or moved-from.");
  return m_textDataSP;
}

bool NC::MatCfg::Impl::compareIgnoringTextDataUID(const MatCfg& o ) const
{
  const Impl * oimpl = &*o.m_impl;
  nc_assert(!this->isMultiPhase());
  nc_assert(!oimpl->isMultiPhase());
  if ( this == oimpl )
    return false;//same internal data instance, must be equal

  if ( this->m_dataSourceName.str() != oimpl->m_dataSourceName.str() )
    return this->m_dataSourceName.str() < oimpl->m_dataSourceName.str();

  //Ok, just the actual parameters left:
  return CfgManip::lessThan( this->m_cfgdata, oimpl->m_cfgdata );
}

bool NC::MatCfg::operator<( const MatCfg& o ) const
{
  if ( this == &o )
    return false;

  //First compare m_phaseChoices if set on at least one:
  if ( !m_impl2->m_phaseChoices.empty() || !o.m_impl2->m_phaseChoices.empty() ) {
    if ( ! ( m_impl2->m_phaseChoices == o.m_impl2->m_phaseChoices ) )
      return m_impl2->m_phaseChoices < o.m_impl2->m_phaseChoices;
  }
  //Then number of (MatCfg-level) phases:
  auto nphase = !m_impl->m_phases ? 1 : m_impl->m_phases->size();
  auto o_nphase = !o.m_impl->m_phases ? 1 : o.m_impl->m_phases->size();
  if ( nphase != o_nphase )
    return nphase < o_nphase;

  if ( ! ( m_impl2->m_densityState == o.m_impl2->m_densityState ) ) {
    if ( m_impl2->m_densityState.has_value() != o.m_impl2->m_densityState.has_value() )
      return m_impl2->m_densityState.has_value();
    nc_assert_always( m_impl2->m_densityState.has_value() && o.m_impl2->m_densityState.has_value() );
    auto type_int   = enumAsInt(   m_impl2->m_densityState.value().type );
    auto o_type_int = enumAsInt( o.m_impl2->m_densityState.value().type );
    if ( type_int != o_type_int )
      return type_int < o_type_int;
    return m_impl2->m_densityState.value().value < o.m_impl2->m_densityState.value().value;
  }

  //Same MatCfg-level phases.
  if ( nphase == 1 ) {
    //Both singlephase (at MatCfg-level).
    if ( m_impl->m_textDataUID != o.m_impl->m_textDataUID )
      return m_impl->m_textDataUID < o.m_impl->m_textDataUID;
    //Phase choices:
    return m_impl->compareIgnoringTextDataUID(o);
  } else {
    //multiphase (at MatCfg-level) with equal number of phases
    nc_assert(m_impl->m_phases!=nullptr);
    nc_assert(o.m_impl->m_phases!=nullptr);
    const PhaseList& pl1 = *(m_impl->m_phases);
    const PhaseList& pl2 = *(o.m_impl->m_phases);

    //First we try to look for different fractions or text uids (because it is cheap):
    for ( auto i : ncrange(nphase) ) {
      if ( pl1.at(i).first != pl2.at(i).first )
        return pl1.at(i).first < pl2.at(i).first;
      if ( pl1.at(i).second.m_impl->m_textDataUID != pl2.at(i).second.m_impl->m_textDataUID )
        return pl1.at(i).second.m_impl->m_textDataUID < pl2.at(i).second.m_impl->m_textDataUID;
    }

    //More full-blown comparison:
    for ( auto i : ncrange(nphase) ) {
      if ( pl1.at(i).second < pl2.at(i).second )
        return true;//not equal and pl1<pl2
      else if ( pl2.at(i).second < pl1.at(i).second )
        return false;//not equal and pl2>pl1
      //pl1==pl2, continue
    }
    //Identical:
    return false;
  }
}

const NC::MatCfg::PhaseList& NC::MatCfg::phases() const
{
  if ( m_impl->m_phases != nullptr )
    return *m_impl->m_phases;
  static NC::MatCfg::PhaseList s_empty_pl;
  return s_empty_pl;
}

bool NC::MatCfg::isSinglePhase() const { return m_impl->m_phases == nullptr; }
bool NC::MatCfg::isMultiPhase() const { return m_impl->m_phases != nullptr; }

std::string NC::MatCfg::toJSONCfg() const
{
  std::ostringstream os;
  os << "{\"format\":\"NCrystal-MatCfg-v1\", \"ismultiphase\":"
     <<(isMultiPhase()?"true":"false");
  if ( !isMultiPhase() ) {
    os << ",\"data_name\":";
    streamJSON(os,m_impl->m_dataSourceName.str());
    os << ",\"textdata_uid\": \""<<m_impl->m_textDataUID.value()<<"\"";
    os << ",\"textdata_type\":";
    streamJSON(os,m_impl->m_textDataType);
    os << ",\"pars\":";
    CfgManip::streamJSON( m_impl->m_cfgdata, os );
  } else {
    //multiphase
    nc_assert(m_impl->m_phases!=nullptr);
    os <<",\"phases\":[";
    for ( auto i : ncrange( m_impl->m_phases->size() ) ) {
      const auto& ph = m_impl->m_phases->at(i);
      os << ( i?",":"");
      os << "[";
      streamJSON(os,ph.first);
      os<<','<<ph.second.toJSONCfg()<<"]";
    }
    os <<']';
  }

  os <<",\"phasechoices\":[";
  auto pcs = getPhaseChoices();
  for ( auto ipc : ncrange(pcs.size()) )
    os << (ipc?",":"") << pcs.at(ipc);
  os<<"]";
  auto ds = get_density();
  os <<",\"density\":{";
  os << "\"type\":\"";
  switch ( ds.type ) {
  case DensityState::Type::SCALEFACTOR:
    os << "scalefactor";
    break;
  case DensityState::Type::DENSITY:
    os << "density_gcm3";
    break;
  case DensityState::Type::NUMBERDENSITY:
    os << "numberdensity_perAa3";
    break;
  default:
    nc_assert_always(false);
    break;
  };
  os << "\",\"value\":";
  streamJSON(os,ds.value);
  os << "}}";
  return os.str();
}

const NC::MatCfg::PhaseChoices& NC::MatCfg::getPhaseChoices() const
{
  return m_impl2->m_phaseChoices;
}

void NC::MatCfg::appendPhaseChoices( const PhaseChoices& pc )
{
  if ( pc.empty() )
    return;
  for ( auto i : pc )
    Impl2::checkPhaseChoiceRange(i);
  auto mod = m_impl2.modify();
  for ( auto i : pc )
    mod->m_phaseChoices.push_back(i);
}

void NC::MatCfg::appendPhaseChoice( unsigned ipc )
{
  Impl2::checkPhaseChoiceRange(ipc);
  m_impl2.modify()->m_phaseChoices.push_back(ipc);
}

NC::MatCfg NC::MatCfg::cloneWithoutPhaseChoices() const
{
  if ( m_impl2->m_phaseChoices.empty() )
    return *this;
  MatCfg cfg_new(*this);
  cfg_new.m_impl2.modify()->m_phaseChoices.clear();
  return cfg_new;
}

NC::MatCfg NC::MatCfg::cloneWithoutDensityState() const
{
  if ( !m_impl2->m_densityState.has_value() )
    return *this;
  MatCfg cfg_new(*this);
  cfg_new.m_impl2.modify()->m_densityState.reset();
  return cfg_new;
}

bool NC::MatCfg::hasDensityOverride() const
{
  return ( m_impl2->m_densityState.has_value()
           && !( m_impl2->m_densityState.value() == DensityState() ) );
}

NC::DensityState NC::MatCfg::get_density() const
{
  return ( m_impl2->m_densityState.has_value()
           ?  m_impl2->m_densityState.value()
           : DensityState() );
}

void NC::MatCfg::set_density( const DensityState& ds ) {
  ds.validate();
  auto newds = accumulateDensityState( m_impl2->m_densityState, ds );
  nc_assert(newds.has_value());
  newds.value().validate();
  if ( ! ( newds == m_impl2->m_densityState ) )
    m_impl2.modify()->m_densityState.set( newds );
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    //C++11 is annoying, hence this helper:
    void set_density_helper( MatCfg& cfg, DensityState::Type type, double val )
    {
      DensityState ds;
      ds.type = type;
      ds.value = val;
      cfg.set_density(ds);
    }
  }
}

void NC::MatCfg::set_density( const Density& val ) { set_density_helper(*this,DensityState::Type::DENSITY, val.dbl()); }
void NC::MatCfg::set_density( const NumberDensity& val ) { set_density_helper( *this, DensityState::Type::NUMBERDENSITY, val.dbl() ); }

const NC::DataSourceName& NC::MatCfg::getDataSourceName() const
{
  if ( m_impl->isMultiPhase() )
    NCRYSTAL_THROW(LogicError,"MatCfg::getDataSourceName called for multiphase object");
  return m_impl->m_dataSourceName;
}

void NC::MatCfg::genDoc( std::ostream& os, GenDocMode gdm )
{
  nc_assert_always(isOneOf(gdm,GenDocMode::TXT_SHORT,GenDocMode::TXT_FULL,GenDocMode::JSON));
  auto cfgdumpmode = ( gdm==GenDocMode::TXT_SHORT ? Cfg::CfgVarListMode::TXT_SHORT
                       : ( gdm==GenDocMode::TXT_FULL ? Cfg::CfgVarListMode::TXT_FULL
                           : Cfg::CfgVarListMode::JSON ) );
  Cfg::dumpCfgVarList( os, cfgdumpmode );
}

void NC::MatCfg::set_ucnmode( const Optional<UCNMode>& v ) { m_impl.modify()->setVar( v, &CfgManip::set_ucnmode ); }
NC::StrView NC::MatCfg::get_ucnmode_str() const { return CfgManip::get_ucnmode_str( m_impl->readVar(Cfg::VarId::ucnmode) ); }
NC::Optional<NC::UCNMode> NC::MatCfg::get_ucnmode() const { return CfgManip::get_ucnmode( m_impl->readVar(Cfg::VarId::ucnmode) ); }

namespace {
  static_assert(std::is_move_constructible<NC::MatCfg>::value, "");
  static_assert(std::is_move_assignable<NC::MatCfg>::value, "");
  static_assert(std::is_copy_constructible<NC::MatCfg>::value, "");
  static_assert(std::is_copy_assignable<NC::MatCfg>::value, "");
}
