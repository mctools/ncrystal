#ifndef NCrystal_ExtinctionCfg_hh
#define NCrystal_ExtinctionCfg_hh

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

#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/internal/utils/NCStrView.hh"

namespace NCRYSTAL_NAMESPACE {

  //////////////////////////////////////////////////////////////
  // Infrastructure needed for the variables in NCCfgVars.hh  //
  //////////////////////////////////////////////////////////////

  namespace Cfg {

    namespace detail {
      struct TV_ph2 {};
      struct TV_ph3 {};
    }

    template<class T1, class T2 = detail::TV_ph2, class T3 = detail::TV_ph3>
    class TrivialVariant {
      //A Variant class for N trivially destructible and copyable types, easy to
      //serialise with a simple bit-wise copy. For now we support N=1, N=2 and
      //N=3. Easy to support others.
      //Fixme: Consider T4, T5, ...
      static_assert(std::is_trivially_destructible<T1>::value,"");
      static_assert(std::is_trivially_copyable<T1>::value,"");
      static_assert(std::is_trivially_destructible<T2>::value,"");
      static_assert(std::is_trivially_copyable<T2>::value,"");
      static_assert(std::is_trivially_destructible<T3>::value,"");
      static_assert(std::is_trivially_copyable<T3>::value,"");
      static_assert(!std::is_same<T1,T2>::value,"");
      static_assert(!std::is_same<T1,T3>::value,"");
      static_assert(!std::is_same<T2,T3>::value,"");

      constexpr static std::size_t maxsize = ncconstexpr_max( sizeof(T1),
                                                              sizeof(T2),
                                                              sizeof(T3) );
      constexpr static std::size_t maxalign = ncconstexpr_lcm( alignof(T1),
                                                               alignof(T2),
                                                               alignof(T3) );
      alignas(maxalign) char m_data[maxsize];
      int m_content = 0;//empty

    public:
      TrivialVariant() = default;
      ~TrivialVariant() = default;
      TrivialVariant( const TrivialVariant& ) = default;
      TrivialVariant& operator=( const TrivialVariant& ) = default;

      TrivialVariant( const T1& o ) { *this = o; }
      TrivialVariant( const T2& o ) { *this = o; }
      TrivialVariant( const T3& o ) { *this = o; }

      TrivialVariant( NullOptType ) {}
      TrivialVariant& operator=( const NullOptType& ) { m_content = 0; return *this; }

      void clear() { m_content = 0; }
      void reset() { clear(); }

      template<class Tn>
      bool has_value() const noexcept {
        static_assert( std::is_same<Tn,T1>::value
                       || std::is_same<Tn,T2>::value
                       || std::is_same<Tn,T3>::value, "" );
        constexpr int valid = ( std::is_same<Tn,T1>::value ? 1
                                : ( std::is_same<Tn,T2>::value ? 2 : 3 ) );
        return m_content == valid;
      }

      template<class Tn>
      const Tn& value() const
      {
        nc_assert_always(has_value<Tn>());//fixme _always
        return *reinterpret_cast<const Tn*>(m_data);
      }

      template<class Tn>
      Tn& value()
      {
        nc_assert_always(has_value<Tn>());//fixme _always
        return *reinterpret_cast<Tn*>(m_data);
      }

      bool empty() const { return m_content == 0; }

      template<class Tn>
      TrivialVariant& operator=( const Tn& o )
      {
        if ( m_data != (char*)(&o) )
          std::memcpy( m_data, (char*)&o, sizeof(Tn) );
        static_assert( std::is_same<Tn,T1>::value
                       || std::is_same<Tn,T2>::value
                       || std::is_same<Tn,T3>::value, "" );
        constexpr int valid = ( std::is_same<Tn,T1>::value ? 1
                                : ( std::is_same<Tn,T2>::value ? 2 : 3 ) );
        m_content = valid;
        return *this;
      }
    };

    static_assert(std::is_trivially_destructible<TrivialVariant<double,int,char>>::value,"");
    static_assert(std::is_trivially_copyable<TrivialVariant<double,int,char>>::value,"");

    template<class TData>
    class TrivialOptional {
      //An Optional class for trivially destructible and copyable types, easy to
      //serialise with a simple bit-wise copy.
      static_assert(std::is_trivially_destructible<TData>::value,"");
      static_assert(std::is_trivially_copyable<TData>::value,"");
      alignas(TData) char m_data[sizeof(TData)];
      bool m_hasData = false;
    public:

      using data_t = TData;

      TrivialOptional() = default;
      ~TrivialOptional() = default;
      TrivialOptional( const TrivialOptional& ) = default;
      TrivialOptional& operator=( const TrivialOptional& ) = default;

      TrivialOptional( const TData& o ) { *this = o; }
      TrivialOptional( NullOptType ) {}
      TrivialOptional& operator=( const NullOptType& ) { m_hasData = false; return *this; }

      void clear() { m_hasData = false; }
      void reset() { clear(); }
      bool has_value() const { return m_hasData; }
      const TData& value() const
      {
        nc_assert(has_value());
        return *reinterpret_cast<const TData*>(m_data);
      }
      TData& value()
      {
        nc_assert(has_value());
        return *reinterpret_cast<TData*>(m_data);
      }

      bool empty() const { return !m_hasData; }
      //char * rawData() const { return &m_data[0]; }

      TrivialOptional& operator=( const TData& o )
      {
        if ( m_data != (char*)(&o) )
          std::memcpy( m_data, (char*)&o, sizeof(TData) );
        m_hasData = true;
        return *this;
      }
    };
    static_assert(std::is_trivially_destructible<TrivialOptional<double>>::value,"");
    static_assert(std::is_trivially_copyable<TrivialOptional<double>>::value,"");

    //Fwd. declare type from NCCfgTypes.hh:
    using ValDbl_ShortStrOrigRep = ShortStr<detail::VarBuf::buffer_local_size
                                            - sizeof(double)>;
    //fixme:

    struct ExtnCfg_Sabine {
      Length blockSize; // l
      struct ScndExtn {
        enum class TiltType { Rectangular, Triangular };
        TiltType tilt_type = TiltType::Rectangular;
        double blockSpread; // G or g
        Length grainSize; // L
        bool use_correlated_model = true;
      };
      TrivialOptional<ScndExtn> secondary;
      static_assert(std::is_trivially_destructible<ScndExtn>::value,"");
      static_assert(std::is_trivially_copyable<ScndExtn>::value,"");
      static_assert(std::is_trivially_destructible<TrivialOptional<ScndExtn>>::value,"");
      static_assert(std::is_trivially_copyable<TrivialOptional<ScndExtn>>::value,"");

    };











    // struct ExtnCfg_CR {
    //   Length blockSize;// l
    //   double blockSpread;// G
    //   Length grainSize;// L
    // };

    // struct ExtnCfg_RED {
    //   Length blockSize; // L
    //   Length bendingRadius // R;
    //   Length mfpElasDeformedRegions; // l
    //   Optional<double> deformationGradient; // c
    // };
    //fixme: to NCTypes.hh in some place:
    // using ValDbl_ShortStrOrigRep = ShortStr<VarBuf::buffer_local_size
    //                                         - sizeof(double)>;/

    // class ExtnCfgData {
    // public:
    //   //Extinction cfg data, but with all the details hidden for ABI
    //   //stability reasons. Except that one can check if extinction is enabled or not.
    //   bool enabled() const { return ! m_data.isEmpty(); }

    //   //Implementation detail, do not use:
    //   const VarBuf& detail_accessRawData() { return m_data; }//fixme: better access control?
    // private:
    //   VarBuf m_data;//Holds an ExtinctionCfg object.
    // };

    class ExtinctionCfg final {
      struct DblVar {
        double value;
        ValDbl_ShortStrOrigRep strrep;//starts with null byte if not available
      };
      struct SabineData {
        SabineData( const ExtnCfg_Sabine& o ) : obj(o) {}
        SabineData() = default;
        ExtnCfg_Sabine obj;
        ValDbl_ShortStrOrigRep origstr_blockSize = NullOpt;
        ValDbl_ShortStrOrigRep origstr_grainSize = NullOpt;
        ValDbl_ShortStrOrigRep origstr_blockSpread = NullOpt;
      };
      using data_t = TrivialVariant<SabineData>;
      data_t m_data = NullOpt;
      static_assert(std::is_trivially_destructible<data_t>::value,"");
      static_assert(std::is_trivially_copyable<data_t>::value,"");

    public:

      //List of models, at most one can be set:
      bool enabled() const { return !m_data.empty(); }

      //Specific models:
      bool has_sabine() const {
        return m_data.has_value<SabineData>();
      }

      const ExtnCfg_Sabine& get_sabine() const {
        if ( !m_data.has_value<SabineData>() )
          NCRYSTAL_THROW(BadInput,"ExtinctionCfg::get_sabine called but Sabine"
                         " model not available. Check .has_sabine() first.");
        return m_data.value<SabineData>().obj;
      }

      //From cfg-string value:
      ExtinctionCfg( StrView );

      //Non-enabled object:
      ExtinctionCfg( no_init_t ) : m_data( NullOpt ) { nc_assert_always(!enabled()); }//fixme =default??

      //Re-encode in form for cfg-string or JSON:
      void stream( std::ostream& ) const;
      void streamJSON( std::ostream& ) const;

      //Disable extinction:
      void clear() { m_data.reset(); }
      void disable() { m_data.reset(); }

      //Assign to specific models:
      ExtinctionCfg( const ExtnCfg_Sabine& sb ) : m_data(SabineData{ sb }) {nc_assert_always(enabled()); }
      ExtinctionCfg& operator=( const ExtnCfg_Sabine& sb ) { m_data = SabineData{ sb }; nc_assert_always(enabled());return *this; }//fixme inline

      ExtinctionCfg( const ExtinctionCfg& ) = default;
      ExtinctionCfg& operator=( const ExtinctionCfg& ) = default;

      ExtinctionCfg( const ExtinctionCfgData&& );

      int cmp( const ExtinctionCfg& ) const;

      //Serialisation to/from VarBuf objects is supported, but is considered a
      //private implementation detail not to be generally used:
      struct detail_from_varbuf_t {};
      ExtinctionCfg( detail_from_varbuf_t, const detail::VarBuf& );
      detail::VarBuf detail_to_varbuf( detail::VarId ) const;

    };
    static_assert(detail::varbuf_calc::buf_align
                  % alignof(ExtinctionCfg) == 0,"");//fixme move this to where we reinterpret_cast


    // class ExtinctionCfg {
    // public:

    //   // Decoded ExtinctionCfg (Immutable, cheap to copy).

    //   //Decode cfg data string (throws BadInput in case of syntax issues):
    //   ExtinctionCfg( const ExtinctionCfgData& );

    //   //Empty, same as decoding an empty cfg data string:
    //   ExtinctionCfg();

    //   //Re-encode into data string:
    //   ExtinctionCfgData encode() const;

    //   //Check if extinction is enabled at all:
    //   bool enabled() const { return m_scale.has_value(); }


    //   //Actual extinction model parameters. At most one of these will have a
    //   //value.

    //   Optional<ExtnCfg_Sabine> model_Sabine() const;

    // private:
    //   struct Impl;
    //   std::shared_ptr<Impl> m_impl;
    // };


// struct ExtnCfg_Sabine {
//   Length blockSize; // l
//   struct ScndExtn {
//     enum class TiltType { Rectangular, Triangular };
//     TiltType tilt_type = TiltType::Rectangular;
//     double blockSpread; // G or g
//     Length grainSize; // L
//     bool use_correlated_model = true;
//   };
//   Optional<ScndExtn> secondary;
// };

// struct ExtnCfg_BC {
//   Length blockSize; // l
//   struct ScndExtn {
//     enum class TiltType { Gaussian, Lorentzian, Fresnel };
//     TiltType tilt_type = TiltType::Gaussian;
//     double blockSpread; // g
//     Length grainSize; // L
//   };
//   Optional<ScndExtn> secondary;
// };

// struct ExtnCfg_CR {
//   Length blockSize;// l
//   double blockSpread;// G
//   Length grainSize;// L
// };

// struct ExtnCfg_RED {
//   Length blockSize; // L
//   Length bendingRadius // R;
//   Length mfpElasDeformedRegions; // l
//   Optional<double> deformationGradient; // c
// };

//   }

// }

  }
}

#endif

// Sabine:

//       //l : block size, Aa
//       //G : integral breadth of the angular distribution of mosaic blocks, dimensionless
//       //L : grain size, Aa (A grain is formed by crystallites or blocks.)
//       //tilt_dist : distribution type for the tilts between mosaic blocks (rect/tri/...). Only needed for secondary.
//       //Alternatively G->g for corr model.

// BC:
//       //tilt_dist /opt : 0 for primary extinction, 1, 2, 3 for sencondary extinction following a
//       //      Gaussian, Lorentzian or Fresnel distribution, respectively
//       //l : "t", mean path length through a perfect crystal, equivalent to block size, Aa
//       //g : width parameter of the mosaic distribution, dimensionless
//       //L : "T-bar", mean path length through a mosaic crystal, Aa

// cooper-rouse:
// //lgL

//       //l : mean free path between elastically deformed regions
//       //R : Bending radius parameter
//       //L : crystal block size
//       //c : deformation gradient parameter
