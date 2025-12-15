#ifndef NCrystal_MMC_EngineOpts_hh
#define NCrystal_MMC_EngineOpts_hh

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
#include "NCrystal/internal/utils/NCHists.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    class TallyFlags final {
    public:

      //Flags for controlling the tally setup. All encoded in a single integer
      //for efficiency.

      struct Flags {
        //Flags enabling certain histograms:
        using value_type = uint32_t;
        constexpr static value_type mu = 0x1;
        constexpr static value_type cosmu = 0x2;
        constexpr static value_type nscat = 0x4;
        constexpr static value_type w = 0x8;
        constexpr static value_type e = 0x10;
        constexpr static value_type l = 0x20;
        constexpr static value_type nscat_uw = 0x40;
        constexpr static value_type q = 0x80;

        //Flags controlling behaviour of all histograms:
        constexpr static value_type nobreakdown = 0x10000;//no comp. breakdown
        constexpr static value_type lowres = 0x20000;//less bins & binrange
        constexpr static value_type highres = 0x40000;//more bins & binrange

        //Collections:
        constexpr static value_type ALLHISTS = 0xFF;//Enable all histograms
        constexpr static value_type ALL = 0x700FF;//all flags set at once
        constexpr static value_type NONE = 0x0;//No flags set
        constexpr static value_type DEFAULT = mu;//Default flag values
      };
      using value_type = Flags::value_type;
      using strlist_type = SmallVector<StrView,8>;
      TallyFlags( value_type value = Flags::DEFAULT );
      TallyFlags( const strlist_type& );
      strlist_type toStringList() const;
      value_type getValue() const;
      bool hasAny( value_type v ) const;
      bool has( value_type v ) const;
      void add( value_type v );
      void remove( value_type v );
      void add( StrView sv ) { add(lookup(sv)); }
      void remove( StrView sv ) { remove(lookup(sv)); }
      static value_type lookup( StrView );
      static const char * singleFlagToString( value_type );
      static bool isSingleFlag( value_type );
    private:
      value_type m_value = Flags::DEFAULT;
    };

    class TallyBinningOverrides final {
      //Keep DB of binnings overridden by user. Any TallyFlags arguments to add
      //or lookup must have a single histogram bit set only!
    public:
      using Binning = Hists::Binning;
      void add( TallyFlags, const Binning& );
      const Binning& lookup( TallyFlags, const Binning& fallback ) const;
      bool operator==( const TallyBinningOverrides& ) const;
    private:
      TallyFlags::value_type m_all = 0;
      SmallVector_IC<std::pair<TallyFlags::value_type,Binning>,4> m_db;
    };

    struct RouletteOptions {
      //Probability of surviving a roulette attempts (would be 5/6 in real
      //russian roulette):
      double survival_probability = 0.1;
      //Particles with weights above this won't get rouletted:
      double weight_threshold = 1e-2;
      //Particles with less than this many scatterings won't get rouletted:
      int nscat_threshold = 2;
      bool operator==( const RouletteOptions& ) const;
    };

    struct EngineOpts {
      //FIXME: Some docs here (or in nctool output / wiki).
      std::uint64_t seed = 0;//Simulation seed. Note that results are only
                             //reproducible in case of single-threaded
                             //simulations, since it is not guaranteed which
                             //particles will be handled by which threads.

      RouletteOptions roulette;
      enum class IgnoreMiss : uint32_t { NO=0, YES=1, Default=NO };
      enum class IncludeAbsorption : uint32_t { NO=0, YES=1, Default=YES };
      IgnoreMiss ignoreMiss = IgnoreMiss::Default;
      ThreadCount nthreads = ThreadCount::auto_detect();
      IncludeAbsorption includeAbsorption = IncludeAbsorption::Default;

      //Limit on number of scatterings. If set, it will contain a value in the
      //range 0..32000. This is the highest number of scatterings that will be
      //modelled for a particle. After that number is reached, further
      //scattering cross sections fpr that particle will artifically become
      //0. As an example, setting nscatlimit=1 will effectively disable multiple
      //scattering effects:
      Optional<unsigned> nScatLimit;

      TallyFlags tallyFlags;
      TallyBinningOverrides tallyBinnings;

      //base of exit angle plots (if not set it will be the source beam
      //direction or [0,0,1] if that is also unavailable):
      Optional<NeutronDirection> tallyBeamDir;

      //Nominal beam direction for tallies needing initial energy. If
      //unavailable, those tallies will be booked but not be filled (and a
      //warning emitted):
      Optional<NeutronEnergy> tallyBeamEnergy;
    };

    //Parse to/from a string representation (can also be used to normalise a
    //string representation). The string representation will be compact, and not
    //include entries at default values:
    EngineOpts parseEngineOpts( StrView );
    std::string engineOptsToString( const EngineOpts& eopts );

    //Default streaming adds "EngineOpts(...)" around the string:
    std::ostream& operator<<( std::ostream& os, const EngineOpts& eopts );

    //Output as a JSON dictionary. This will always include all options,
    //including those at default values.
    void engineOptsToJSON(std::ostream&, const EngineOpts&);

  }
}
////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {
    inline TallyFlags::TallyFlags( value_type value )
      : m_value(value)
    {
      nc_assert_always( (m_value & Flags::ALL)==m_value );
    }

    inline TallyFlags::value_type TallyFlags::getValue() const
    {
      return m_value;
    }

    inline bool TallyFlags::hasAny( value_type v ) const
    {
      nc_assert( (v & Flags::ALL)==v && v );
      return v & m_value;
    }

    inline bool TallyFlags::has( value_type v ) const
    {
      nc_assert( (v & Flags::ALL)==v );
      return ( v & m_value ) == v;
    }

    inline void TallyFlags::add( value_type v )
    {
      nc_assert_always( (v & Flags::ALL)==v );
      m_value |= v;
    }

    inline void TallyFlags::remove( value_type v )
    {
      nc_assert_always( (v & Flags::ALL)==v );
      m_value &= ~v;
    }

    inline void TallyBinningOverrides::add( TallyFlags f, const Binning& b )
    {
      auto vt = f.getValue();
      nc_assert( vt & TallyFlags::Flags::ALLHISTS );
      nc_assert( TallyFlags::singleFlagToString( vt ) );
      if ( vt & m_all ) {
        for ( auto& e : m_db ) {
          if ( e.first == vt ) {
            e.second = b;
            return;
          }
        }
      } else {
        m_db.emplace_back( vt, b );
        m_all |= vt;
      }
    }

    inline const TallyBinningOverrides::Binning&
    TallyBinningOverrides::lookup( TallyFlags f, const Binning& fb ) const
    {
      auto vt = f.getValue();
      nc_assert( vt & TallyFlags::Flags::ALLHISTS );
      //throws if not single flag:
      nc_assert( TallyFlags::singleFlagToString( vt ) );
      if ( vt & m_all ) {
        for ( auto& e : m_db )
          if ( e.first == vt )
            return e.second;
      }
      return fb;
    }

    inline bool
    TallyBinningOverrides::operator==( const TallyBinningOverrides& o ) const
    {
      return m_all == o.m_all && m_db == o.m_db;
    }

    inline bool RouletteOptions::operator==( const RouletteOptions& o ) const
    {
      return  ( survival_probability == o.survival_probability
                && weight_threshold == o.weight_threshold
                && nscat_threshold == o.nscat_threshold );
    }

  }
}

#endif
