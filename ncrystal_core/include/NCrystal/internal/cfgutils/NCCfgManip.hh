#ifndef NCrystal_CfgManip_hh
#define NCrystal_CfgManip_hh

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

#include "NCrystal/internal/cfgutils/NCCfgVars.hh"

namespace NCRYSTAL_NAMESPACE {

  /////////////////////////////////////////////////////////////////////////////////
  //                                                                             //
  // CfgManip: Internal helper class which is responsible for parsing parameter  //
  // strings and maintaining lists of set variables in a CfgDataobject (cf.      //
  // NCCfgVars.hh for variable definitions).                                     //
  //                                                                             //
  // The internal data in CfgData objects is essentially just a heterogenic      //
  // sorted SBO vector of data of enabled variable's (again using SBO buffers).  //
  // Thus, the need for malloc's is kept low, and no storage is used for the     //
  // unset variables (which is usually the vast majority of them). It is also    //
  // efficient to create filtered objecs with only variables from certain groups //
  // (e.g. Info, Scatter, ...) present.                                          //
  //                                                                             //
  // Parameter string parsing happens via the applyStrCfg function, which also   //
  // parses recognised-but-unhandled top-level parameters (density and           //
  // phasechoice) to the calling code.                                           //
  //                                                                             //
  // The design allows the MatCfg class to construct the objects, apply string   //
  // configs, and handle top level parameters. From the MatCfg interface it is   //
  // then possible to create request objects relevant for specific factories     //
  // (InfoRequest, ScatterRequest, AbsorptionRequest).                           //
  //                                                                             //
  /////////////////////////////////////////////////////////////////////////////////

  namespace Cfg {

    class TopLvlVar {
    public:
      //Contains either density modification or phase choice:
      static_assert(std::is_same<unsigned,std::underlying_type<DensityState::Type>::type>::value,"");
      bool isDensity() const noexcept;
      DensityState getDensity() const;
      bool isPhaseChoice() const noexcept;
      unsigned getPhaseChoice() const;

      //Constructors:
      TopLvlVar( DensityState );
      struct phasechoice_t {};
      TopLvlVar( phasechoice_t, unsigned );

      static constexpr unsigned phasechoice_max = 10000;
    private:
      double m_dbl;//<0 means phase choice, otherwise it is density value;
      unsigned m_uint;//phase choice or density type
    };

    using TopLvlVarList = SmallVector<TopLvlVar,6>;
    using VarIdFilter = std::function<bool(VarId)>;
    using VarIdList = SmallVector_IC<VarId,8>;

    class CfgManip {
    public:
      static_assert(std::is_same<std::underlying_type<VarId>::type,std::uint32_t>::value,"");
      using size_type = std::size_t;

      //Create new objects with all except selected variables unset (removed):
      static CfgData filter( const CfgData&, VarIdFilter );

      static bool filterSelectsAny( const CfgData&, VarIdFilter );

      enum class FilterType { ExcludeListed, OnlyListed };
      static VarIdFilter createFilter( const VarIdList& vl, FilterType );

      //Cfg-string i/o. The applyStrCfg method parses strings like like
      //"par1=val1;par2=val2;..." (with optional whitespace at either sides of
      //pari and vali), and calls setVarByStr(..) accordingly. Both methods
      //throw BadInput in case of syntax errors, unknown variable names, or
      //forbidden characters. Note that the applyStrCfg returns a list of any
      //top level parameters encountered (it otherwise ignores them).
      ncnodiscard17 static TopLvlVarList applyStrCfg( CfgData&, StrView );

      //Output stream support:
      static void stream( const CfgData&, std::ostream&, const VarIdFilter& filter );
      static void stream( const CfgData& cd, std::ostream& os ) { stream(cd, os, nullptr); }

      //Stream set variables as json list, e.g.: [ ["temp", 293.15], ["vdoslus", 3] ]
      static void streamJSON( const CfgData&, std::ostream& );

      //Miscellaneous:
      static VarIdList findCommonEntries( std::function<const CfgData*()> cfgDataIter );

      //Comparisons:
      static bool equal( const CfgData&, const CfgData& );
      static bool lessThan( const CfgData&, const CfgData& );//unspecified ordering

      //Applies variables contained in other CfgData object to object:
      static void apply( CfgData&, const CfgData&, VarIdFilter = nullptr );

      //Check how many variables are set:
      ncnodiscard17 static ncconstexpr17 size_type size(const CfgData&) noexcept;
      ncnodiscard17 static ncconstexpr17 bool empty(const CfgData&) noexcept;
      ncnodiscard17 static bool empty( const CfgData&, const VarIdFilter& );//check if empty after filtering

      //Check if given variable has a value set:
      static bool hasValueSet( const CfgData&, VarId );

      //Check if a given variable has same value in two objects. This is true
      //only when both are either unset or set to the same value. If one is
      //unset and the other is set to a value which happens to be the default
      //value, this method will still return false:
      static bool hasSameValue( const CfgData&, const CfgData&, VarId );

    private:
      static void setVarByStr( CfgData&, VarId, StrView );
      template<class TVarDef>
      static decltype(TVarDef::get_val(VarBuf(NullOpt))) getValue(const CfgData&);
      template<class TVarDef, class TValType = typename TVarDef::value_type>
      static void setValue( CfgData&, const TValType& val );
      template<class TVarDef>
      static decltype(TVarDef::get_val(VarBuf(NullOpt))) getValueFromBufPtr(const VarBuf*);
    public:

      //Access individual variables programatically (for C++ API):
      static Temperature get_temp(const CfgData& data) { return Temperature{ getValue<vardef_temp>(data) }; }
      static void set_temp( CfgData& data, Temperature val ) { setValue<vardef_temp>(data,val.dbl()); }

      static bool get_coh_elas(const CfgData& data) { return getValue<vardef_coh_elas>(data); }
      static void set_coh_elas( CfgData& data, bool val ) { setValue<vardef_coh_elas>(data,val); }

      static bool get_incoh_elas(const CfgData& data) { return getValue<vardef_incoh_elas>(data); }
      static void set_incoh_elas( CfgData& data, bool val ) { setValue<vardef_incoh_elas>(data,val); }

      static bool get_sans(const CfgData& data) { return getValue<vardef_sans>(data); }
      static void set_sans( CfgData& data, bool val ) { setValue<vardef_sans>(data,val); }

      static double get_dcutoff(const CfgData& data) { return getValue<vardef_dcutoff>(data); }
      static void set_dcutoff( CfgData& data, double val ) { setValue<vardef_dcutoff>(data,val); }

      static double get_dcutoffup(const CfgData& data) { return getValue<vardef_dcutoffup>(data); }
      static void set_dcutoffup( CfgData& data, double val ) { setValue<vardef_dcutoffup>(data,val); }

      static double get_dirtol(const CfgData& data) { return getValue<vardef_dirtol>(data); }
      static void set_dirtol( CfgData& data, double val ) { setValue<vardef_dirtol>(data,val); }

      static MosaicityFWHM get_mos(const CfgData& data) { return MosaicityFWHM{ getValue<vardef_mos>(data) }; }
      static void set_mos( CfgData& data, MosaicityFWHM val ) { setValue<vardef_mos>(data,val.dbl()); }

      static OrientDir get_dir1(const CfgData& data) { return getValue<vardef_dir1>(data); }
      static void set_dir1( CfgData& data, const OrientDir& val ) { setValue<vardef_dir1>(data,val); }

      static OrientDir get_dir2(const CfgData& data) { return getValue<vardef_dir2>(data); }
      static void set_dir2( CfgData& data, const OrientDir& val ) { setValue<vardef_dir2>(data,val); }

      static double get_mosprec(const CfgData& data) { return getValue<vardef_mosprec>(data); }
      static void set_mosprec( CfgData& data, double val ) { setValue<vardef_mosprec>(data,val); }

      static double get_sccutoff(const CfgData& data) { return getValue<vardef_sccutoff>(data); }
      static void set_sccutoff( CfgData& data, double val ) { setValue<vardef_sccutoff>(data,val); }

      static int get_vdoslux(const CfgData& data) { return static_cast<int>( getValue<vardef_vdoslux>(data) ); }
      static void set_vdoslux( CfgData& data, int val ) { setValue<vardef_vdoslux>( data, static_cast<std::int64_t>(val) ); }

      static std::int_least32_t get_lcmode(const CfgData& data) { return static_cast<std::int_least32_t>( getValue<vardef_lcmode>(data) ); }
      static void set_lcmode( CfgData& data, std::int_least32_t val ) { setValue<vardef_lcmode>( data,static_cast<std::int_least32_t>(val) ); }

      static StrView get_ucnmode_str(const CfgData& data) { return getValue<vardef_ucnmode>(data); }
      static Optional<UCNMode> get_ucnmode( const CfgData& data ) { return vardef_ucnmode::decode_value( get_ucnmode_str(data) ); }
      static void set_ucnmode( CfgData& data, const Optional<UCNMode>& val ) {
        if ( !val.has_value() ) {
          setValue<vardef_ucnmode>(data,"");
        } else {
          std::ostringstream ss;
          ss << val.value();
          setValue<vardef_ucnmode>(data,ss.str());
        }
      }

      static const LCAxis& get_lcaxis(const CfgData& data) { return getValue<vardef_lcaxis>(data).as<LCAxis>(); }
      static void set_lcaxis( CfgData& data, const LCAxis& val ) { setValue<vardef_lcaxis>(data, val.as<Vector>()); }

      static StrView get_infofactory(const CfgData& data) { return getValue<vardef_infofactory>(data); }
      static void set_infofactory( CfgData& data, StrView val ) { setValue<vardef_infofactory>(data,val); }
      static void set_infofactory_stdstr( CfgData& data, const std::string& val ) { setValue<vardef_infofactory,std::string>(data,val); }
      static void set_infofactory_cstr( CfgData& data, const char * val ) { setValue<vardef_infofactory,const char *>(data,val); }

      static StrView get_scatfactory(const CfgData& data) { return getValue<vardef_scatfactory>(data); }
      static void set_scatfactory( CfgData& data, StrView val ) { setValue<vardef_scatfactory>(data,val); }
      static void set_scatfactory_stdstr( CfgData& data, const std::string& val ) { setValue<vardef_scatfactory,std::string>(data,val); }
      static void set_scatfactory_cstr( CfgData& data, const char * val ) { setValue<vardef_scatfactory,const char *>(data,val); }

      static StrView get_absnfactory(const CfgData& data) { return getValue<vardef_absnfactory>(data); }
      static void set_absnfactory( CfgData& data, StrView val ) { setValue<vardef_absnfactory>(data,val); }
      static void set_absnfactory_stdstr( CfgData& data, const std::string& val ) { setValue<vardef_absnfactory,std::string>(data,val); }
      static void set_absnfactory_cstr( CfgData& data, const char * val ) { setValue<vardef_absnfactory,const char *>(data,val); }

      static StrView get_atomdb(const CfgData& data) { return getValue<vardef_atomdb>(data); }
      static std::vector<VectS> get_atomdb_parsed(const CfgData&);
      static void set_atomdb_parsed( CfgData& data, std::vector<VectS> val)
      {
        std::string tmp;
        constexpr StrView svSC = StrView::make(":");
        for (auto& e : val ) {
          if ( !tmp.empty() )
            tmp += '@';
          tmp.append(joinstr(e,svSC));
        }
        set_atomdb_stdstr(data,std::move(tmp));
      }
      static void set_atomdb( CfgData& data, StrView val ) { setValue<vardef_atomdb>(data,val); }
      static void set_atomdb_stdstr( CfgData& data, const std::string& val ) { setValue<vardef_atomdb,std::string>(data,val); }
      static void set_atomdb_cstr( CfgData& data, const char * val ) { setValue<vardef_atomdb,const char *>(data,val); }

      static StrView get_inelas(const CfgData& data) { return getValue<vardef_inelas>(data); }
      static void set_inelas( CfgData& data, StrView val ) { setValue<vardef_inelas>(data,val); }
      static void set_inelas_stdstr( CfgData& data, const std::string& val ) { setValue<vardef_inelas,std::string>(data,val); }
      static void set_inelas_cstr( CfgData& data, const char * val ) { setValue<vardef_inelas,const char *>(data,val); }

      //Various utilities:
      static bool isSingleCrystal(const CfgData&);
      static bool isLayeredCrystal(const CfgData&);

      //Info, ScatterBase, ScatterExtra, Absorption

      //Check consistency of parameters in various groups. This implements
      //checks affecting more than a single parameter
      //(e.g. "dcutoff<=dcutoffup", but not "dcutoff>=0"), and results in
      //BadInput exceptions in case of inconsistencies.
      static void checkParamConsistency_Info( const CfgData& );
      static void checkParamConsistency_ScatterBase( const CfgData& );
      static void checkParamConsistency_ScatterExtra( const CfgData& );
      static void checkParamConsistency_Absorption( const CfgData& );

      //Validate orientation parameters (throws BadInput exception if not valid,
      //or if not oriented) and return SCOrientation object. Using template to
      //avoid direct dependencies here:
      template<class TSCOrientation>
      ncnodiscard17 static TSCOrientation createSCOrientation(const CfgData&);

    private:
      template<class TBufCreateFct>
      static void detail_setVar( CfgData& data, VarId varid, TBufCreateFct& doCreate );
      static const VarBuf* searchBuf( const CfgData& data, VarId );
      static VarBuf* searchBuf( CfgData& data, VarId );
    };

    std::ostream& operator<<( std::ostream&, const CfgData& );

  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace Cfg {

    ncnodiscard17 inline ncconstexpr17 CfgManip::size_type CfgManip::size( const CfgData& data ) noexcept
    {
      return data().size();
    }

    ncnodiscard17 inline ncconstexpr17 bool CfgManip::empty( const CfgData& data ) noexcept
    {
      return data().empty();
    }

    inline bool CfgManip::filterSelectsAny( const CfgData& data, VarIdFilter filter )
    {
      for ( auto& e : data() ) {
        if ( filter( static_cast<VarId>(e.metaData()) ) )
          return true;
      }
      return false;
    }

    ncnodiscard17 inline bool CfgManip::empty( const CfgData& data, const VarIdFilter& filter )
    {
      if ( data().empty() || !filter )
        return data().empty();
      for ( auto& e : data() ) {
        if ( filter( static_cast<VarId>(e.metaData()) ) )
          return false;
      }
      return true;
    }

    inline std::ostream& operator<<( std::ostream& os, const CfgData& cfgData )
    {
      CfgManip::stream(cfgData,os);
      return os;
    }

    inline const VarBuf* CfgManip::searchBuf( const CfgData& data, VarId varid )
    {
      auto it = std::lower_bound( data().begin(), data().end(), varid,
                                  []( const VarBuf& a, const VarId& b ) { return a.metaData() < b; } );
      return ( it != data().end() && it->metaData() == varid ) ? it : nullptr;
    }

    inline VarBuf* CfgManip::searchBuf( CfgData& data, VarId varid )
    {
      auto it = std::lower_bound( data().begin(), data().end(), varid,
                                  []( const VarBuf& a, const VarId& b ) { return a.metaData() < b; } );
      return ( it != data().end() && it->metaData() == varid ) ? it : nullptr;
    }

    inline bool CfgManip::hasSameValue( const CfgData& a, const CfgData& b, VarId varid )
    {
      auto buf_a = searchBuf( a, varid );
      auto buf_b = searchBuf( b, varid );
      if ( !buf_a && !buf_b )
        return true;//both unset
      if ( buf_a && buf_b ) {
        //both set:
        return 0 == varInfo( varid ).bufCmp(*buf_a,*buf_b);
      }
      return false;
    }

    inline bool CfgManip::hasValueSet( const CfgData& data, VarId varid )
    {
      return searchBuf(data,varid) != nullptr;
    }

    inline bool CfgManip::isLayeredCrystal(const CfgData& data)
    {
      return hasValueSet( data, Cfg::VarId::lcaxis );
    }

    template<class TVarDef, class TValType>
    inline void CfgManip::setValue( CfgData& data, const TValType& val )
    {
      constexpr auto varid = constexpr_varIdFromName(TVarDef::name);
      //NB: no need to capture &varid (even gives warning on clang/gcc), but
      //VSCode wants it:
#ifdef _MSC_VER
      auto doCreate = [&val,&varid](){ return TVarDef::set_val( varid, val ); };
#else
      auto doCreate = [&val](){ return TVarDef::set_val( varid, val ); };
#endif
      detail_setVar( data, varid, doCreate );
    }

    template<class TVarDef>
    inline decltype(TVarDef::get_val(VarBuf(NullOpt))) CfgManip::getValue(const CfgData& data)
    {
      constexpr auto varid = constexpr_varIdFromName(TVarDef::name);
      return getValueFromBufPtr<TVarDef>(searchBuf(data,varid));
    }

    template<class TVarDef>
    inline decltype(TVarDef::get_val(VarBuf(NullOpt))) CfgManip::getValueFromBufPtr(const VarBuf* bufptr)
    {
      if ( bufptr )
        return TVarDef::get_val( *bufptr );
      if ( !TVarDef::has_default_value() )
        NCRYSTAL_THROW2(MissingInfo,"Value for parameter "<<TVarDef::name<<" not available");
      //Copy default value to static var since we might be returning by reference:
      static typename TVarDef::value_type s_def_val = TVarDef::default_value_or_nonsense();
      return s_def_val;
    }

    inline void CfgManip::setVarByStr( CfgData& data, VarId varid, StrView str )
    {
      auto doCreate = [&varid,&str](){ return varInfo(varid).from_str( varid, str.trimmed() ); };
      detail_setVar( data, varid, doCreate );
    }

    template<class TBufCreateFct>
    inline void CfgManip::detail_setVar( CfgData& data, VarId varid, TBufCreateFct& doCreate )
    {
      auto it = std::lower_bound( data().begin(), data().end(), varid,
                                  []( const VarBuf& a, const VarId& b ) { return a.metaData() < b; } );
      if ( it != data().end() && it->metaData() == varid ) {
        //Override existing:
        *it = doCreate();
        return;
      }

      //Insert in properly sorted location...
      if ( it == data().end() ) {
        //... at the end
        data().emplace_back( doCreate() );
        return;
      }
      //... at the position "it". We must move existing entries to make space. To
      //avoid issues with iterator invalidation we first append an unset entry, then
      //revalidate the "it" iterator, and finally move all entries past "it" one
      //step towards the end:
      nc_assert(!data().empty());
      auto tgtpos = std::distance( data().begin(), it );//"it" might get invalidated during emplace_back
      data().emplace_back( NullOpt );
      it = std::next( data().begin(), tgtpos );//make sure "it" is valid again
      auto itEmptySlot = std::prev( data().end() );
      while ( itEmptySlot > it ) {
        auto it2 = std::prev(itEmptySlot);
        *itEmptySlot  = std::move( *it2 );
        itEmptySlot = it2;
      }
      nc_assert( itEmptySlot == it );
      *it = doCreate();
    }


    inline TopLvlVar::TopLvlVar( DensityState ds )
      : m_dbl( ds.value ), m_uint( enumAsInt(ds.type) )
    {
      ds.validate();
    }

    inline TopLvlVar::TopLvlVar( phasechoice_t, unsigned pc ) : m_dbl(-1.0), m_uint(pc)
    {
      if ( pc > phasechoice_max )
        NCRYSTAL_THROW2(BadInput,"Invalid phase choice index (too high): "<<pc);
    }
    inline bool TopLvlVar::isPhaseChoice() const noexcept { return m_dbl < 0.0; }
    inline bool TopLvlVar::isDensity() const noexcept { return !isPhaseChoice(); }
    inline DensityState TopLvlVar::getDensity() const
    {
      nc_assert( isDensity() );
      nc_assert( m_dbl > 0.0 && !std::isinf(m_dbl) );
      nc_assert( m_uint==0 || m_uint==1 || m_uint==2 );
      return DensityState{ static_cast<DensityState::Type>(m_uint), m_dbl };
    }

    inline unsigned TopLvlVar::getPhaseChoice() const
    {
      nc_assert( isPhaseChoice() );
      return m_uint;
    }

    template<class TSCOrientation>
    ncnodiscard17 inline TSCOrientation CfgManip::createSCOrientation( const CfgData& data )
    {
      //Validate the 4 relevant parameters and create SCOrientation object (if relevant).
      auto bufptr_mos = searchBuf( data, VarId::mos );
      auto bufptr_dir1 = searchBuf( data, VarId::dir1 );
      auto bufptr_dir2 = searchBuf( data, VarId::dir2 );

      int nOrient = (bufptr_dir1?1:0) + (bufptr_dir2?1:0) + (bufptr_mos?1:0);
      if (nOrient!=0 && nOrient<3)
        NCRYSTAL_THROW(BadInput,"Must set all or none of mos, dir1 and dir2 parameters");

      if ( !nOrient ) {
        if ( hasValueSet( data, VarId::dirtol ) )
          NCRYSTAL_THROW(BadInput,"mos, dir1 and dir2 parameters must all be set when dirtol is set");
        NCRYSTAL_THROW(BadInput,"Can only create SCOrientation object for oriented configurations");
      }

      //mos already validated:
      nc_assert( vardef_mos::get_val( *bufptr_mos ) > 0.0 && vardef_mos::get_val( *bufptr_mos ) <= kPiHalf );

      TSCOrientation orient;
      orient.setPrimaryDirection( vardef_dir1::get_val( *bufptr_dir1 ) );
      orient.setSecondaryDirection( vardef_dir2::get_val( *bufptr_dir2 ), get_dirtol(data) );
      nc_assert_always(orient.isComplete());
      return orient;
    }
  }
}

#endif
