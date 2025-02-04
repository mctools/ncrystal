#ifndef NCrystal_CfgTypes_hh
#define NCrystal_CfgTypes_hh

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

#include "NCrystal/core/NCImmutBuf.hh"
#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/core/NCVariant.hh"

namespace NCRYSTAL_NAMESPACE {

  //////////////////////////////////////////////////////////////
  // Infrastructure needed for the variables in NCCfgVars.hh  //
  //////////////////////////////////////////////////////////////

  namespace Cfg {

    using VarId = detail::VarId;
    using VarBuf = detail::VarBuf;
    using VarBufVector = detail::VarBufVector;

    static constexpr StrView forbidden_chars_non_multiphase = StrView::make("\"'|(){}[]*&><#$");
    static constexpr StrView forbidden_chars_value_strreps = StrView::make("\"'|(){}[]*&><;=#$");//+ "=;"
    static constexpr StrView forbidden_chars_multiphase = StrView::make("\"'|(){}[]#$");//*&<> are allowed only in certain locations of multiphase strings

    enum class VarGroupId { Info, ScatterBase, ScatterExtra, Absorption };

    //A few utilities for doubles:
    inline double sanitiseDblValue(double v,const char * varname) {
      if ( std::isnan(v) )
        NCRYSTAL_THROW2(BadInput,"NAN (not-a-number) value provided"
                        " for parameter \""<<varname<<"\"");
      return ( v == 0.0 ? 0.0 : v );
    }
    using ValDbl_ShortStrOrigRep = ShortStr<VarBuf::buffer_local_size - sizeof(double)>;
    struct units_notavail {
      static constexpr auto actual_unit = nullptr;
      static void listAvailableUnits(std::ostream&){}
      static Optional<std::pair<double,ValDbl_ShortStrOrigRep>> parse(StrView){ return NullOpt; }
    };

    void standardInputStrSanityCheck( const char * parname, StrView strrep );

    template <class Derived, class TValueType>
    class ValBase {
    public:
      using units = units_notavail;
      using value_type = TValueType;

      constexpr static bool has_default_value()
      {
        static_assert ( std::is_same<decltype(Derived::default_value()),value_type>::value
                        || std::is_same<decltype(Derived::default_value()),NullOptType>::value, "" );
        return std::is_same<decltype(Derived::default_value()),value_type>::value;
      }

      static VarBuf default_value_in_varbuf()
      {
        nc_assert(has_default_value());
        struct detail2 {
          static VarBuf default_varbuf( NullOptType ) noexcept { return VarBuf(NullOpt); }//not actually called
          static VarBuf default_varbuf( const value_type& val ) noexcept { return Derived::set_val( static_cast<VarId>(0), val ); }//dummy varid
        };
        return detail2::default_varbuf( Derived::default_value() );
      }

      static void stream_default_value( std::ostream& os )
      {
        //This does not need to be super effiecient, we only need it to generate
        //documentation. It is more important to be consistent, and for that
        //reason we first place the default value in a VarBuf and then stream it.
        if (!has_default_value())
          return;
        Derived::stream_val(os,default_value_in_varbuf());
      }
      struct detail_dvon {
        //Next fct won't actually get called
        static constexpr value_type dvon( NullOptType )
        {
          return nc_assert_always_rv(false),*(value_type*)(nullptr);// cppcheck-suppress nullPointer
        }
        static constexpr value_type dvon( const value_type& v ) { return v; }
      };
      constexpr static value_type default_value_or_nonsense()
      {
        return detail_dvon::dvon( Derived::default_value() );
      }

      static void stream_default_value_json( std::ostream& os )
      {
        //Like stream_default_value, but with proper type for json. Outputs json
        //"null" if there is no default value.
        if ( has_default_value() )
          Derived::asJSONObject( os, default_value_in_varbuf() );
        else
          streamJSON(os,json_null_t{});
      }

    };

    template <class Derived>
    class ValDbl : public ValBase<Derived,double> {
      // As ever, floating point numbers are a bit weird and annoying in
      // several ways, even for IEEE 754 doubles: Contain values that prevent
      // proper sorting (NaN), bit patterns that are not unique for a given
      // value (-0 and +0), string representations that are not unique for a
      // given value ("0.1" and "0.99999999999999999999999"). To prevent bugs
      // and weird user experiences, we adopt the following strategy. First of
      // all, we always store the actual floating point (double) value in the
      // first sizeof(double) chars of the buffer. The following chars
      // optionally provides a null-terminated c-string with the string
      // representation of the value. However, this might be absent
      // (i.e. start with a null char) in case the string representation does
      // not fit in the local buffer, in which case a string representation
      // will be created on-the-fly from the double value when needed.
      //
      // A few take-aways:
      //
      // 1) To avoid weird sorting and other issues: Forbid NAN, and map -0 to 0.
      // 2) Always store the double itself in first 8 bytes
      // 3) Store str-rep in last bytes if and only if it fits there (and when
      //    setting from string, prefer the original version only if it is shorter).
      // 4) When sorting, first sort the double value, then presence and
      //    contentsof str-reps.
      // 5) When streaming and strrep is missing, generate it with the same
      //    method call as during setting.

      static constexpr auto buf_maxstrlen = VarBuf::buffer_local_size - sizeof(double) - 1;//-1 for null char
      static_assert(VarBuf::buffer_local_size>=detail::varbuf_calc::ValDbl_buf_minsize,"");
      static_assert(VarBuf::buffer_local_size>=sizeof(double),"");
      static_assert(alignof(VarBuf)%alignof(double)==0,"");
      static VarBuf detail_set_val_with_str( VarId varid, double value, StrView valstr )
      {
        //Assumes value already sanitised, validated, and the shortest str rep was found.
        char tmp[VarBuf::buffer_local_size];
        std::memcpy(&tmp[0],(void*)&value,sizeof(double));
        std::size_t nused = sizeof(double) + 1;//at least the double and a null char
        if ( valstr.size() <= buf_maxstrlen ) {
          std::memcpy( &tmp[sizeof(double)], valstr.data(), valstr.size() );
          tmp[sizeof(double)+valstr.size()] = '\0';
          nused += valstr.size();//also the characters
        } else {
          //don't store strrep, indicate with a single null char
          tmp[sizeof(double)] = '\0';
        }
        return VarBuf( &tmp[0], nused, varid );
      }
    public:
      using Base = ValDbl;
      static_assert( ValDbl_ShortStrOrigRep::bufsize == buf_maxstrlen + 1,"");

      static constexpr auto value_type_descr = "floating point number";

      static double value_validate( double value ) { return value; }

      static int cmp( const VarBuf& a, const VarBuf& b )
      {
        //First just sort on values:
        double valA = get_val( a );
        double valB = get_val( b );
        //NB: sanitiseDblValue guarantees no NaNs:
        nc_assert( !std::isnan(valA) );
        nc_assert( !std::isnan(valB) );
        if ( valA != valB )
          return ( valA < valB ? -1 : 1 );
        //valA = valB. To consider the full state (and thus avoid weird issues),
        //we now also consider the str-rep's if any:
        const char * cstr_A = ( a.dataAssertLocal() + sizeof(double) );
        const char * cstr_B = ( b.dataAssertLocal() + sizeof(double) );
        if ( bool(cstr_A) != bool(cstr_B) )
          return bool(cstr_A) ? -1 : 1;//those with strreps before those without
        if (!cstr_A)
          return 0;//both strreps are absent
        //Both strreps are there:
        return std::strcmp(cstr_A,cstr_B);
      }

      static double get_val( const VarBuf& buf )
      {
        return *reinterpret_cast<const double*>(buf.dataAssertLocal());
      }

      static VarBuf set_val( VarId varid, double value ) {
        value = Derived::value_validate( sanitiseDblValue(value,Derived::name) );
        return detail_set_val_with_str( varid, value, dbl2shortstr(value).to_view() );
      }
      static void stream_val( std::ostream& os, const VarBuf& buf )
      {
        const char * buf_data = buf.dataAssertLocal();
        const char * cstr = ( buf_data + sizeof(double) );
        if ( *cstr ) {
          //str rep (including null termination) was short enough to keep cached:
          os << cstr;
        } else {
          //str rep was not short enough to cache, generate on demand (passing
          //through same method as the cached would have been, for
          //consistency):
          os << dbl2shortstr( *reinterpret_cast<const double*>(buf_data) );
        }
      }

      static void asJSONObject( std::ostream& os, const VarBuf& b )
      {
        streamJSON(os,get_val(b));
      }

      static VarBuf from_str( VarId varid, StrView orig_str ) {
        standardInputStrSanityCheck(Derived::name,orig_str);
        auto opt_val = Derived::units::parse(orig_str);//this allows derived class to handle units and ranges
        if (!opt_val.has_value())
          NCRYSTAL_THROW2(BadInput,"Syntax error - invalid value \""
                          <<orig_str<<"\" provided for parameter \""<<Derived::name<<"\"");
        double value = opt_val.value().first;
        const auto& orig_str_cleaned = opt_val.value().second;
        value = Derived::value_validate( sanitiseDblValue( value, Derived::name ) );
        auto valstr = dbl2shortstr( value );
        return detail_set_val_with_str( varid,
                                        value,
                                        ( ( orig_str_cleaned.empty() || valstr.size() <= orig_str_cleaned.size() )
                                          ? valstr.to_view()
                                          : orig_str_cleaned.to_view() ) );
      }
    };

    template <class Derived>
    class ValVector : public ValBase<Derived,Vector> {
    public:
      //Simply stores a Vector directly in the buffer.
      static_assert(VarBuf::buffer_local_size>=sizeof(Vector),"");
      static_assert(alignof(VarBuf)%alignof(Vector)==0,"");
      static_assert(std::is_trivially_copyable<Vector>::value,"");
      static_assert(std::is_trivially_destructible<Vector>::value,"");

      static constexpr bool auto_normalise = true;//will normalise all vectors to unit length unless overridden to false.
      using Base = ValVector;
      static constexpr auto value_type_descr = "vector (3D)";
      static const Vector& get_val( const VarBuf& buf )
      {
        return *reinterpret_cast<const Vector*>( buf.dataAssertLocal() );
      }
      static int cmp( const VarBuf& a, const VarBuf& b )
      {
        return get_val( a ).lexCmp( get_val( b ) );
      }
      static VarBuf set_val( VarId varid, const Vector& value ) {
        static_assert(std::is_trivially_destructible<Vector>::value,"");
        static_assert(std::is_trivially_copyable<Vector>::value,"");
        alignas(Vector) char tmp[sizeof(Vector)];
        Vector & v = *(new (&tmp[0]) Vector(no_init));
        v[0] = sanitiseDblValue( value[0], Derived::name );
        v[1] = sanitiseDblValue( value[1], Derived::name );
        v[2] = sanitiseDblValue( value[2], Derived::name );
        if ( Derived::auto_normalise ) {
          double m2 = v.mag2();
          if ( !m2 )
            NCRYSTAL_THROW2(BadInput,"Null vector provided for parameter \""<<Derived::name<<"\"");
          v.normalise();//slower but safer than "v /= std::sqrt(m2)"
        }
        //Extra sanity checks:
        for ( auto& e : v )
          e = sanitiseDblValue( e, Derived::name );
        Derived::extraChecks( v );
        return VarBuf( &tmp[0], sizeof(tmp), varid );
      }

      static void stream_val( std::ostream& os, const VarBuf& buf )
      {
        auto& v = get_val(buf);
        os << fmt( v[0] ) << ',' << fmt( v[1] ) << ',' << fmt( v[2] );
      }

      static void asJSONObject( std::ostream& os, const VarBuf& b )
      {
        streamJSON(os,get_val(b));
      }

      static VarBuf from_str( VarId varid, StrView sv_in ) {
        standardInputStrSanityCheck(Derived::name,sv_in);
        Vector v(no_init);
        auto sv = sv_in;
        for ( auto idx : ncrange(3) ) {
          auto i = sv.find(',');
          Optional<double> opt_dblval;
          if ( ( i == StrView::npos ) == ( idx == 2 ) ) {
            opt_dblval = sv.substr(0,i).trimmed().toDbl();
          }
          if ( !opt_dblval.has_value() )
            NCRYSTAL_THROW2(BadInput,"Syntax error - invalid value \""
                            <<sv_in<<"\" provided for parameter \""<<Derived::name<<"\"");
          v[idx] = opt_dblval.value();
          sv = sv.substr( i < StrView::npos ? i+1 : StrView::npos );//+1 to skip past comma
        }
        return set_val( varid, v );
      }

      static void extraChecks( const Vector& ) {}
    };

    template <class Derived>
    class ValOrientDir : public ValBase<Derived,OrientDir> {
    public:
      //Store OrientDir by converting to 2*3=6 doubles and a bool (1
      //char). Presumably this will always be too large for the local buffer,
      //but usage is rare.
      static_assert(alignof(VarBuf)%alignof(OrientDir)==0,"");
      using Base = ValOrientDir;
      static constexpr auto value_type_descr = "crystal axis orientation";

      struct DecodedBuffer{ const double* vals; bool crystal_is_hkl; };
      static DecodedBuffer detail_bufdecode( const VarBuf& buf )
      {
        const char* bufdata = buf.data();
        return { reinterpret_cast<const double*>( bufdata ),
                 bool( *std::next(bufdata,6*sizeof(double)) == 1 ) };
      }

      static OrientDir get_val( const VarBuf& buf )
      {
        auto db = detail_bufdecode(buf);
        Vector vcrys( db.vals[0], db.vals[1], db.vals[2] );
        LabAxis vlab( db.vals[3], db.vals[4], db.vals[5] );
        if ( db.crystal_is_hkl )
          return { vcrys.as<HKLPoint>(), vlab };
        else
          return { vcrys.as<CrystalAxis>(), vlab };
      }

      static int cmp( const VarBuf& a, const VarBuf& b )
      {
        auto db_A = detail_bufdecode(a);
        auto db_B = detail_bufdecode(b);
        if ( db_A.crystal_is_hkl != db_B.crystal_is_hkl )
          return db_A.crystal_is_hkl ? -1 : 1;
        for ( auto i : ncrange(6) ) {
          if ( db_A.vals[i] != db_B.vals[i] )
            return db_A.vals[i] < db_B.vals[i] ? -1 : 1;
        }
        return 0;
      }

      static VarBuf set_val( VarId varid, const OrientDir& value ) {
        //nc_assert_always(!value.crystal.empty());//fails if moved-from
        alignas(double) char tmp[6*sizeof(double)+1];
        const Vector * vcrys;
        if ( value.crystal.has_value<HKLPoint>() ) {
          tmp[6*sizeof(double)] = 1;
          vcrys = &value.crystal.get<HKLPoint>().as<Vector>();
        } else if ( value.crystal.has_value<CrystalAxis>() ) {
          tmp[6*sizeof(double)] = 0;
          vcrys = &value.crystal.get<CrystalAxis>().as<Vector>();
        } else {
          NCRYSTAL_THROW2(BadInput,"Moved-from crystal direction object"
                          " provided for parameter \""<<Derived::name<<"\"");
        }
        if ( std::min<double>( value.lab.as<Vector>().mag2(), vcrys->mag2() ) < 1e-100 )
          NCRYSTAL_THROW2(BadInput,"Null vector provided for parameter \""<<Derived::name<<"\"");
        double * tmp_vals = reinterpret_cast<double*>( &tmp[0] );
        tmp_vals[0] = sanitiseDblValue( (*vcrys)[0], Derived::name );
        tmp_vals[1] = sanitiseDblValue( (*vcrys)[1], Derived::name );
        tmp_vals[2] = sanitiseDblValue( (*vcrys)[2], Derived::name );
        tmp_vals[3] = sanitiseDblValue( value.lab[0], Derived::name );
        tmp_vals[4] = sanitiseDblValue( value.lab[1], Derived::name );
        tmp_vals[5] = sanitiseDblValue( value.lab[2], Derived::name );
        return VarBuf( &tmp[0], sizeof(tmp), varid );
      }

      static void stream_val( std::ostream& os, const VarBuf& buf )
      {
        auto db = detail_bufdecode(buf);
        os << ( db.crystal_is_hkl?"@crys_hkl:":"@crys:")
           << fmt(db.vals[0]) << "," << fmt(db.vals[1]) << ","
           << fmt(db.vals[2]) << "@lab:" << fmt(db.vals[3]) << ","
           << fmt(db.vals[4]) << "," << fmt(db.vals[5]);
      }

      static void asJSONObject( std::ostream& os, const VarBuf& buf )
      {
        auto db = detail_bufdecode(buf);
        os << "{\"crystal_is_hkl\":";
        os << (db.crystal_is_hkl?"true":"false");
        os<<",\"crystal\":[";
        streamJSON(os,db.vals[0]);
        os << ",";
        streamJSON(os,db.vals[1]);
        os << ",";
        streamJSON(os,db.vals[2]);
        os << "], \"lab\":[";
        streamJSON(os,db.vals[3]);
        os << ",";
        streamJSON(os,db.vals[4]);
        os << ",";
        streamJSON(os,db.vals[5]);
        os << "]}";
      }

      static VarBuf from_str( VarId varid, StrView sv_in ) {
        standardInputStrSanityCheck(Derived::name,sv_in);
        //"@crys:c1,c2,c3@lab:l1,l2,l3" or "@crys_hkl:c1,c2,c3@lab:l1,l2,l3",
        //and spaces allowed around all delimeters or at string ends.

        //Split into parts:
        const char delims[] = "@:,,@:,,";
        constexpr unsigned ndelims = sizeof(delims)-1;
        static_assert(ndelims==8,"");
        StrView parts[ndelims+1];
        StrView sv = sv_in;
        auto throw_error = [&sv_in]() {
          NCRYSTAL_THROW2(BadInput,"Syntax error - invalid value \""
                          <<sv_in<<"\" provided for parameter \""<<Derived::name<<"\"");
        };
        for ( auto ipart : ncrange(ndelims) ) {
          auto i = sv.find(delims[ipart]);
          if ( i == StrView::npos )
            throw_error();
          parts[ipart] = sv.substr(0,i).trimmed();
          sv = sv.substr(i+1);
        }
        parts[ndelims] = sv.trimmed();

        //Parse parts:
        if (!parts[0].empty())
          throw_error();//part before the initial '@'.
        if ( parts[5] != StrView::make("lab") )
          throw_error();

        if ( ! ( parts[1]==StrView::make("crys")
                 || parts[1] == StrView::make("crys_hkl")) )
          throw_error();

        const bool is_crystal_hkl = (parts[1]==StrView::make("crys_hkl"));

        LabAxis lab{no_init};
        Vector crystal{no_init};
        auto trySetVal = [&throw_error](double& val, const StrView& sv_val) {
          auto opt_val = sv_val.toDbl();
          if ( !opt_val.has_value() )
            throw_error();
          val = opt_val.value();
        };
        trySetVal(crystal[0],parts[2]);
        trySetVal(crystal[1],parts[3]);
        trySetVal(crystal[2],parts[4]);
        trySetVal(lab[0],parts[6]);
        trySetVal(lab[1],parts[7]);
        trySetVal(lab[2],parts[8]);
        if ( is_crystal_hkl )
          return set_val( varid, OrientDir{ crystal.as<HKLPoint>(), lab } );
        else
          return set_val( varid, OrientDir{ crystal.as<CrystalAxis>(), lab } );
      }
    };

    template <class Derived>
    class ValStr : public ValBase<Derived,StrView> {
    public:
      using Base = ValStr;
      static constexpr auto value_type_descr = "string";//NB: Don't change this without follow up
                                                        //changes elsewhere (used to recognise type
                                                        //needing quotation marks in printouts)

      //Encode value in buffer - including null char.
      static StrView get_val( const VarBuf& buf )
      {
        return StrView( buf.data() );
      }

      static int cmp( const VarBuf& a, const VarBuf& b )
      {
        StrView valA = get_val( a );
        StrView valB = get_val( b );
        return valA == valB ? 0 : ( valA < valB ? -1 : 1 );
      }

      static VarBuf set_val( VarId varid, StrView value ) {
        return actual_set_val( varid, value );
      }

      static VarBuf set_val( VarId varid, const std::string& value ) {
        return actual_set_val( varid, StrView(value) );
      }

      static VarBuf set_val( VarId varid, const char * value ) {
        return actual_set_val( varid, StrView(value) );
      }

      static VarBuf from_str( VarId varid, StrView sv_in ) {
        return actual_set_val( varid, sv_in );
      }

      static void stream_val( std::ostream& os, const VarBuf& buf )
      {
        os << buf.data();
      }

      static void asJSONObject( std::ostream& os, const VarBuf& b )
      {
        streamJSON(os,get_val(b));
      }

    protected:
      static Variant<StrView,std::string> str2val( StrView sv ) { return sv; }
    private:
      static VarBuf actual_set_val( VarId varid, StrView sv_in ) {
        standardInputStrSanityCheck(Derived::name,sv_in);
        //NB: Not sure why "auto opt_value = ..." in the next line does not work!!
        Variant<StrView,std::string> opt_value = Derived::str2val(sv_in);
        static_assert(std::is_same<Variant<StrView,std::string>,decltype(opt_value)>::value,"");
        if (opt_value.empty())
          NCRYSTAL_THROW2(BadInput,"Syntax error - invalid value \""
                          <<sv_in<<"\" provided for parameter \""<<Derived::name<<"\"");

        if ( opt_value.has_value<std::string>() ) {
          const auto& value = opt_value.get<std::string>();
          return VarBuf( &value[0], value.size()+1, varid );//+1 for null terminator!
        } else {
          nc_assert(opt_value.has_value<StrView>());
          const auto& value = opt_value.get<StrView>();
          //Must copy in order to null terminate (NB: we could possibly avoid
          //this by a bit of care and an additional constructor on VarBuf):
          SmallVector<char,256> tmp;
          tmp.setByCopy( value.begin(),value.end() );
          tmp.push_back('\0');
          return VarBuf( &tmp[0], tmp.size(), varid );
        }
      }
    };

    template <class Derived>
    class ValBool : public ValBase<Derived,bool> {
      //ValBuf is just holding a single char with the value 0 or 1.
    public:
      using Base = ValBool;
      static constexpr auto value_type_descr = "boolean";

      static int cmp( const VarBuf& a, const VarBuf& b )
      {
        bool valA = get_val( a );
        bool valB = get_val( b );
        return valA == valB ? 0 : ( valA ? -1 : 1 );
      }

      static bool get_val( const VarBuf& buf )
      {
        static_assert(VarBuf::buffer_local_size>=1,"");
        return 0 != buf.dataAssertLocal()[0];
      }
      static VarBuf set_val( VarId varid, bool value ) {
        char c(value?1:0);
        return VarBuf( (const char*)&c, 1, varid );
      }
      static void stream_val( std::ostream& os, const VarBuf& buf )
      {
        os << (get_val(buf)?"1":"0");
      }

      static void asJSONObject( std::ostream& os, const VarBuf& b )
      {
        streamJSON(os,get_val(b));
      }

      static VarBuf from_str( VarId varid, StrView sv ) {
        standardInputStrSanityCheck(Derived::name,sv);
        if ( sv == StrView::make("true") || sv == StrView::make("1") ) {
          return set_val( varid, true );
        } else if ( sv == StrView::make("false") || sv == StrView::make("0") ) {
          return set_val( varid, false );
        } else {
          NCRYSTAL_THROW2(BadInput,"Could not convert \""<<sv
                          <<"\" to boolean value (should be \"true\", \"1\", \"false\" or \"0\")");
        }
      }
    };

    template <class Derived>
    class ValInt : public ValBase<Derived,int64_t> {
      //ValBuf is just holding the bytes of a int64_t value directly (we don't
      //have the same value<->string issues as for a floating point number).
    public:
      using Base = ValInt;
      static constexpr auto value_type_descr = "integer";

      static int cmp( const VarBuf& a, const VarBuf& b )
      {
        int64_t valA = get_val( a );
        int64_t valB = get_val( b );
        return valA == valB ? 0 : ( valA < valB ? -1 : 1 );
      }

      static int64_t get_val( const VarBuf& buf )
      {
        static_assert(VarBuf::buffer_local_size>=sizeof(int64_t),"");
        static_assert(alignof(VarBuf)%alignof(int64_t)==0,"");
        return *reinterpret_cast<const int64_t*>(buf.dataAssertLocal());
      }
      static VarBuf set_val( VarId varid, int64_t value ) {
        value = Derived::value_validate(value);
        return VarBuf( (const char*)&value, sizeof(int64_t), varid );
      }

      static int64_t value_validate( int64_t value ) { return value; }

      static void stream_val( std::ostream& os, const VarBuf& buf )
      {
        os << get_val(buf);
      }

      static void asJSONObject( std::ostream& os, const VarBuf& b )
      {
        streamJSON(os,get_val(b));
      }

      static Optional<int64_t> str2val( StrView sv ) { return sv.toInt(); }

      static VarBuf from_str( VarId varid, StrView str ) {
        standardInputStrSanityCheck(Derived::name,str);
        Optional<int64_t> opt_val = Derived::str2val(str);
        if (!opt_val.has_value())
          NCRYSTAL_THROW2(BadInput,"Syntax error - invalid value \""
                          <<str<<"\" provided for parameter \""<<Derived::name<<"\"");
        return set_val( varid, opt_val.value() );
      }
    };

    using VarFromStrFct = VarBuf(*)(VarId, StrView);
    using StreamVarFct  = void(*)(std::ostream&, const VarBuf&);
    using BufCmpFct = int(*)(const VarBuf&,const VarBuf&);
    using BufToJSONFct = void(*)(std::ostream&, const VarBuf&);

    class VarInfo final {
    public:

      using StreamFct = void (*)( std::ostream& );

      constexpr VarInfo( VarFromStrFct fs,
                         StreamVarFct sv,
                         BufCmpFct bcf,
                         VarGroupId grId,
                         const char * thename,
                         const char * thedescr,
                         StreamFct streamDefaultFct,
                         StreamFct streamDefaultFctJSON,
                         StreamFct streamUnitsFct,
                         const char * actualUnitName,
                         const char * valueTypeDescr,
                         BufToJSONFct jsonfct )
        : m_fromStrFct(fs),
          m_streamFct(sv),
          m_bcf(bcf),
          m_groupId(grId),
          m_namecstr(thename),
          m_name(StrView::constexpr_t(),thename),
          m_description(thedescr),
          m_sdf(streamDefaultFct),
          m_sdfjson(streamDefaultFctJSON),
          m_suf(streamUnitsFct),
          m_aun(actualUnitName),
          m_vtd(valueTypeDescr),
          m_jsonfct(jsonfct)
      {
      }
      VarBuf from_str( VarId varid, StrView sv ) const { return m_fromStrFct(varid,sv); }
      void stream( std::ostream& os, const VarBuf& buf) const { m_streamFct(os,buf); }
      void streamAsJSON( std::ostream& os, const VarBuf& buf ) const { nc_assert(m_jsonfct); m_jsonfct(os,buf); }
      int bufCmp( const VarBuf& a, const VarBuf& b ) const { return m_bcf(a,b); }
      constexpr VarGroupId groupId() const noexcept { return m_groupId; }
      constexpr const char * name() const noexcept { return m_namecstr; }
      constexpr const StrView& nameSV() const noexcept { return m_name; }
      constexpr const char * description() const noexcept { return m_description; }
      constexpr bool hasDefaultValue() const noexcept { return m_sdf != nullptr; }
      void streamDefaultValue( std::ostream& os ) const { nc_assert(m_sdf); m_sdf(os); }
      void streamDefaultValueJSON( std::ostream& os ) const { nc_assert(m_sdfjson); m_sdfjson(os); }
      constexpr bool hasUnits() const noexcept { return m_suf != nullptr; }
      void streamUnitsDescr( std::ostream& os ) const { nc_assert(m_suf); m_suf(os); }
      const char * actualUnitName() const { return m_aun; }
      const char * valueTypeDescr() const { return m_vtd; }

    private:
      VarFromStrFct m_fromStrFct;
      StreamVarFct m_streamFct;
      BufCmpFct m_bcf;
      VarGroupId m_groupId;
      const char * m_namecstr;
      StrView m_name;//NB: Could potentially keep name in local buffer, for
                     //better search locality (ShortStr<32>, but requires
                     //constexpr ShortStr constructor)?
      const char * m_description;
      StreamFct m_sdf;
      StreamFct m_sdfjson;
      StreamFct m_suf;
      const char * m_aun;
      const char * m_vtd;
      BufToJSONFct m_jsonfct;
    };

    template<class TVarDef>
    inline constexpr VarInfo make_varinfo()
    {
      static_assert( constexpr_strlen(TVarDef::name) >= 3,"" );
      static_assert( constexpr_strlen(TVarDef::name) <= 17,"" );//"absorptionfactory" is 17
      return VarInfo( TVarDef::from_str,
                      TVarDef::stream_val,
                      TVarDef::cmp,
                      TVarDef::group,
                      TVarDef::name,
                      TVarDef::description,
                      (TVarDef::has_default_value() ? TVarDef::stream_default_value : nullptr),
                      (TVarDef::has_default_value() ? TVarDef::stream_default_value_json : nullptr),
                      ( TVarDef::units::actual_unit!=nullptr ? TVarDef::units::listAvailableUnits : VarInfo::StreamFct(nullptr) ),
                      TVarDef::units::actual_unit,
                      TVarDef::value_type_descr,
                      TVarDef::asJSONObject );
    }

    struct units_temperature {
      static constexpr auto actual_unit = "K";
      static void listAvailableUnits(std::ostream&);
      static Optional<std::pair<double,ValDbl_ShortStrOrigRep>> parse(StrView sv);
    };

    struct units_length {
      static constexpr auto actual_unit = "Aa";
      static void listAvailableUnits(std::ostream&);
      static Optional<std::pair<double,ValDbl_ShortStrOrigRep>> parse(StrView sv);
    };

    struct units_angle {
      static constexpr auto actual_unit = "rad";
      static void listAvailableUnits(std::ostream&);
      static Optional<std::pair<double,ValDbl_ShortStrOrigRep>> parse(StrView sv);
    };

    struct units_purenumberonly {
      static constexpr const char* actual_unit = nullptr;
      static void listAvailableUnits(std::ostream&){}
      static Optional<std::pair<double,ValDbl_ShortStrOrigRep>> parse(StrView sv);
    };

    namespace detail {
      template<class T>
      inline constexpr unsigned constexpr_name2Idx_helper( const T& list, const char * name, unsigned i)
      {
        return constexpr_cstrequal(name,list[i].name()) ? i : constexpr_name2Idx_helper(list, name,i+1);
      }
    }

    template<class T>
    inline constexpr unsigned constexpr_name2Idx( const T& list, const char * name )
    {
      return detail::constexpr_name2Idx_helper<T>(list,name,0);
    }

  }
}

#endif
