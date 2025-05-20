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
#include "NCrystal/internal/utils/NCTrivOptVar.hh"

namespace NCRYSTAL_NAMESPACE {

  ////////////////////////////////////////////////
  // Types needed for extinction configuration  //
  ////////////////////////////////////////////////

  namespace Cfg {

    //Fwd. declare type from NCCfgTypes.hh (fixme):
    using ValDbl_ShortStrOrigRep = ShortStr<detail::VarBuf::buffer_local_size-sizeof(double)>;

    namespace detail {
      using StrOrigRep = ValDbl_ShortStrOrigRep;

      template<class TData, int NStrReps>
      struct ParsedDataHolder {
        TData obj;
      };
    }

    template<class TClass>
    struct FPValHolder {
      TClass value;
      detail::StrOrigRep detail_orig_str_rep;
    };

    struct ExtnCfg_Generic {
      FPValHolder<Length> domainSize;
      struct Grain {
        FPValHolder<Length> grainSize;
        FPValHolder<double> angularSpread;
      };
      Utils::TrivialOptional<Grain> grain;
    };

    struct ExtnCfg_Sabine {
      ExtnCfg_Generic generic;//the common parts (including str reps)
      //Extra parameters for secondary extinction:
      enum class TiltType { Rectangular, Triangular };
      TiltType tilt_type = TiltType::Rectangular;
      bool use_correlated_model = true;
    };

    class ExtinctionCfg final {
      // struct SabineData {
      //   SabineData( const ExtnCfg_Sabine& o ) : obj(o) {}
      //   SabineData() = default;
      //   ExtnCfg_Sabine obj;
      //   ValDbl_ShortStrOrigRep origstr_blockSize = NullOpt;
      //   ValDbl_ShortStrOrigRep origstr_grainSize = NullOpt;
      //   ValDbl_ShortStrOrigRep origstr_blockSpread = NullOpt;
      // };
      using data_t = Utils::TrivialVariant<ExtnCfg_Sabine>;

      //        SabineData>;
      //      ValDbl_ShortStrOrigRep m_common_strreps[3];
      data_t m_data = NullOpt;
      static_assert(std::is_trivially_destructible<data_t>::value,"");
      static_assert(std::is_trivially_copyable<data_t>::value,"");

    public:

      //List of models, at most one can be set:
      bool enabled() const { return !m_data.empty(); }

      //Specific models:
      bool has_sabine() const {
        return m_data.has_value<ExtnCfg_Sabine>();
      }

      const ExtnCfg_Sabine& get_sabine() const {
        if ( !m_data.has_value<ExtnCfg_Sabine>() )
          NCRYSTAL_THROW(BadInput,"ExtinctionCfg::get_sabine called but Sabine"
                         " model not available. Check .has_sabine() first.");
        return m_data.value<ExtnCfg_Sabine>();
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
      ExtinctionCfg( const ExtnCfg_Sabine& sb ) : m_data(sb) {nc_assert_always(enabled()); }
      ExtinctionCfg& operator=( const ExtnCfg_Sabine& sb ) { m_data = sb; nc_assert_always(enabled());return *this; }//fixme inline

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
  }
}

#endif
