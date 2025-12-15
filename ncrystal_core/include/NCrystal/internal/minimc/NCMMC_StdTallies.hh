#ifndef NCrystal_MMC_StdTallies_hh
#define NCrystal_MMC_StdTallies_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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
#include "NCrystal/internal/minimc/NCMMC_Tally.hh"
#include "NCrystal/internal/utils/NCHists.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    template<class TBasket>
    class TallyStdHists final : public Tally<TBasket> {
    public:
      using hist_t = Hists::Hist1D<Hists::AllowWeights::YES,
                                   Hists::OverflowHandling::Record,
                                   std::vector<double>>;//fixme: smallvector?

      enum class DetailedHistsID : std::size_t {
        NOSCAT = 0,
        SINGLESCAT_ELAS,
        SINGLESCAT_INELAS,
        MULTISCAT_PUREELAS,
        MULTISCAT_OTHER,//fixme: confusing name. MULTISCAT_MIX ?
        END
      };
      static constexpr unsigned ndetail
      = static_cast<unsigned>(DetailedHistsID::END);
      static_assert( ndetail == 5, "" );
    private:
      using this_class_t = TallyStdHists;
      using TFlags = TallyFlags::Flags;
      //Keep all configuration options on a single object (for robust cloning):
      struct Options {
        TallyFlags flags;
        Optional<NeutronDirection> beamDir;
        Optional<NeutronEnergy> beamEnergy;
        TallyBinningOverrides tallyBinnings;
        bool operator==( const Options& o ) const {
          return ( flags.getValue() == o.flags.getValue()
                   && beamDir == o.beamDir
                   && beamEnergy == o.beamEnergy
                   && tallyBinnings == o.tallyBinnings );
        }
      };
      Options m_opt;

      struct HistGroup final : public MoveOnly {
        //Each quantity tallied is kept in a group of 1 histogram (all
        //neutrons), and if enabled also 5 histograms for contribution
        //breakdowns. We equip it with a few utilities for merging and filling.
        hist_t main;
        SmallVector<hist_t,5> detailed;
        void merge( const HistGroup& o )
        {
          nc_assert( main.title() == o.main.title() );
          nc_assert( detailed.size() == o.detailed.size() );
          main.merge( o.main );
          for ( std::size_t i = 0; i < detailed.size(); ++i ) {
            nc_assert( detailed.at(i).title() == o.detailed.at(i).title() );
            detailed[i].merge( o.detailed[i] );
          }
        }
        HistGroup( const Hists::Binning& b, const char * title)
          : main( b, title )
        {
        }

        void bookDetailedHists()
        {
          nc_assert( detailed.empty() );
          auto binning = main.binning();
          detailed.emplace_back(binning,"NOSCAT");
          detailed.emplace_back(binning,"SINGLESCAT_ELAS");
          detailed.emplace_back(binning,"SINGLESCAT_INELAS");
          detailed.emplace_back(binning,"MULTISCAT_PUREELAS");
          detailed.emplace_back(binning,"MULTISCAT_OTHER");
        }

        template<class TVarType = double>
        void fill_main( const TVarType* x, const double* w, std::size_t n)
        {
          //no surprises:
          static_assert( std::is_same<TVarType,double>::value
                         || std::is_same<TVarType,int>::value, "");
          if ( w != nullptr ) {
            for ( std::size_t i = 0; i < n; ++i)
              main.fill( x[i], w[i] );
          } else {
            for ( std::size_t i = 0; i < n; ++i)
              main.fill( x[i] );
          }
        }

        template<class TVarType = double>
        void fill_detail( const DetailedHistsID* hid,
                          const TVarType* x, const double* w, std::size_t n )
        {
          //no surprises:
          static_assert( std::is_same<TVarType,double>::value
                         || std::is_same<TVarType,int>::value, "");
          if ( w != nullptr ) {
            for ( std::size_t i = 0; i < n; ++i)
              vectAt(detailed,static_cast<std::size_t>(hid[i]))
                .fill(static_cast<double>(x[i]),w[i]);
          } else {
            for ( std::size_t i = 0; i < n; ++i)
              vectAt(detailed,static_cast<std::size_t>(hid[i]))
                .fill(static_cast<double>(x[i]));
          }
        }
      };

      //The histograms. Including a cache of indices for efficient later access:
      std::vector<HistGroup> m_hists;
      std::size_t m_histidx_mu = std::numeric_limits<std::size_t>::max();
      std::size_t m_histidx_cosmu = std::numeric_limits<std::size_t>::max();
      std::size_t m_histidx_nscat = std::numeric_limits<std::size_t>::max();
      std::size_t m_histidx_w = std::numeric_limits<std::size_t>::max();
      std::size_t m_histidx_e = std::numeric_limits<std::size_t>::max();
      std::size_t m_histidx_l = std::numeric_limits<std::size_t>::max();
      std::size_t m_histidx_nscat_uw = std::numeric_limits<std::size_t>::max();
      std::size_t m_histidx_q = std::numeric_limits<std::size_t>::max();

      //Utility function for determining number of bins:
      static determineNbins( TallyFlags f, unsigned nbase = 150 )
      {
        nc_assert_always( nbase%2 == 0);
        if ( f.has(TFlags::lowres) )
          nbase /= 2;
        if ( f.has(TFlags::highres) )
          nbase *= 10;
        return nbase;
      }
    public:

      //Convenience constructor for initialising from engine and source options
      //(could in principle also depend on geometry options):
      TallyStdHists( const EngineOpts& eo, const Source& src )
        : TallyStdHists([&eo,&src]()
        {
          Options opt;
          //Tally options from EngineOpts:
          opt.flags = eo.tallyFlags;
          opt.tallyBinnings = eo.tallyBinnings;
          //Nominal beam dir + energy from source:
          Optional<NeutronDirection> optbd = eo.tallyBeamDir;
          if ( !optbd.has_value() )
            optbd = src.nominalBeamDirection();
          if ( optbd.has_value() ) {
            auto v = optbd.value().as<Vector>().unit();
            if ( v.z() != 1.0 )
              opt.beamDir = v.as<NeutronDirection>();
          }
          opt.beamEnergy = eo.tallyBeamEnergy;
          if ( !opt.beamEnergy.has_value() )
            opt.beamEnergy = src.nominalBeamEnergy();
          return opt;
        }())
      {
      }

      TallyStdHists( Options opt )
        : m_opt(opt)
      {

        if ( !m_opt.beamEnergy.has_value() ) {
          //fixme: actually test this branch (maybe allow
          //       enginecfg.beamenergy=none to remove actual src.beamEnergy?)!
          constexpr TallyFlags::value_type hists_needing_Ei = TFlags::q;
          auto bf = m_opt.flags.getValue() & hists_needing_Ei;
          if ( bf ) {
            std::string bfl(joinstr(TallyFlags(bf).toStringList(),"\", \""));
            NCRYSTAL_WARN("Disabling some tallies due to missing knowledge of "
                          "incoming energy: \""<<bfl<<"\" (consider provi"
                          "ding one with the enginecfg beamenergy parameter).");
            m_opt.flags.remove( bf );
          }
        }

        const bool breakdown = !m_opt.flags.has(TFlags::nobreakdown);
        auto addhist = [this]( TallyFlags::value_type f,
                               std::size_t& idxvar,
                               hist_t::size_type nbins,
                               double xmin, double xmax,
                               bool fixnbins = false )
        {
          //sanity check that this is not already booked:
          nc_assert( idxvar == std::numeric_limits<std::size_t>::max());
          //Determine binning automatically, but allow user overrides:
          if (!fixnbins)
            nbins = determineNbins( m_opt.flags, nbins );
          Hists::Binning b_auto( nbins, xmin, xmax);
          idxvar = m_hists.size();
          m_hists.emplace_back( m_opt.tallyBinnings.lookup( TallyFlags(f),
                                                            b_auto ),
                                TallyFlags::singleFlagToString(f) );
        };
        constexpr auto nbins_std = 200;
        constexpr auto nbins_angular = 180;
        //Fixme for the TFlags::w histogram (i.e. non-weighted hist), the
        //fractions shown in the breakdown plots are very misleading.
        if ( m_opt.flags.has(TFlags::mu) )
          addhist( TFlags::mu, m_histidx_mu, nbins_angular, 0.0, 180.0 );
        if ( m_opt.flags.has(TFlags::cosmu) )
          addhist( TFlags::cosmu, m_histidx_cosmu, nbins_std, -1.0, 1.0 );
        if ( m_opt.flags.hasAny(TFlags::nscat|TFlags::nscat_uw) ) {
          unsigned nscatmax = ( m_opt.flags.has(TFlags::lowres) ? 10 : 20 );
          if (m_opt.flags.has(TFlags::highres))
            nscatmax += 20;
          if ( m_opt.flags.has(TFlags::nscat ) )
            addhist( TFlags::nscat, m_histidx_nscat,
                     nscatmax+2,-1.5, nscatmax+0.5, true );
          if ( m_opt.flags.has(TFlags::nscat_uw ) )
            addhist( TFlags::nscat_uw, m_histidx_nscat_uw,
                     nscatmax+2,-1.5, nscatmax+0.5, true );
        }
        if ( m_opt.flags.has(TFlags::w) )
          addhist( TFlags::w, m_histidx_w, nbins_std, 0.0, 1.0 );//fixme: logw?
                                                                 //else upper
                                                                 //limit should
                                                                 //be related to
                                                                 //initial src
                                                                 //weight if
                                                                 //known (maybe
                                                                 //we should for
                                                                 //now simply
                                                                 //disallow
                                                                 //changing
                                                                 //initial
                                                                 //weights).
        //fixme: better max of e/l/q tallys? emax = 3*max(beamE,kT)?
        if ( m_opt.flags.has(TFlags::e) )
          addhist( TFlags::e, m_histidx_e, nbins_std, 0.0, 1.0 );
        if ( m_opt.flags.has(TFlags::l) )
          addhist( TFlags::l, m_histidx_l, nbins_std, 0.0, 10.0 );
        if ( m_opt.flags.has(TFlags::q) )
          addhist( TFlags::q, m_histidx_q, nbins_std, 0.0, 10.0 );
        if ( !m_opt.flags.has(TFlags::nobreakdown) ) {
          for ( auto& h : m_hists )
            h.bookDetailedHists();
        }
      }

      using basket_t = TBasket;
      void registerResults( const basket_t& b ) override
      {
        const std::size_t n = b.size();
        //fixme: only find cos_exit_angle if has tally: q, mu, or cosmu
        double cos_exit_angle[basket_N];
        double exit_angle[basket_N];
        if ( !m_opt.beamDir.has_value() ) {
          //use (0,0,1) as reference:
          for ( std::size_t i = 0; i < n; ++i )
            cos_exit_angle[i] = b.neutrons.uz[i];//fixme; just a memcopy? (or better yet, just reassign a pointer)
        } else {
          const double bx( m_opt.beamDir.value()[0] );
          const double by( m_opt.beamDir.value()[1] );
          const double bz( m_opt.beamDir.value()[2] );
          for ( std::size_t i = 0; i < n; ++i )
            cos_exit_angle[i] = bx * b.neutrons.ux[i];
          for ( std::size_t i = 0; i < n; ++i )
            cos_exit_angle[i] += by * b.neutrons.uy[i];
          for ( std::size_t i = 0; i < n; ++i )
            cos_exit_angle[i] += bz * b.neutrons.uz[i];
        }

#ifndef NDEBUG
        for ( std::size_t i = 0; i < n; ++i ) {
          nc_assert_always( cos_exit_angle[i] > -(1.0+1e-14) );
          nc_assert_always( cos_exit_angle[i] <  (1.0+1e-14) );
        }
#endif
        for ( std::size_t i = 0; i < n; ++i )
          cos_exit_angle[i] = ncclamp(cos_exit_angle[i],-1.0,1.0);

        auto&w = b.neutrons.w;

        using nscat_t_cv = typename std::remove_all_extents<decltype(b.cache.nscat)>::type;
        using nscat_t = typename std::remove_cv<nscat_t_cv>::type;
        static_assert( std::is_same<nscat_t,int>::value, "");

        if ( m_opt.flags.has(TFlags::nscat) )
          vectAt(m_hists,m_histidx_nscat).fill_main( b.cache.nscat, w, n );

        if ( m_opt.flags.has(TFlags::nscat_uw) )
          vectAt(m_hists,m_histidx_nscat_uw).fill_main( b.cache.nscat,
                                                        nullptr, n );

        if ( m_opt.flags.has(TFlags::cosmu) )
          vectAt(m_hists,m_histidx_cosmu).fill_main( cos_exit_angle, w, n );

        if ( m_opt.flags.has(TFlags::w) )
          vectAt(m_hists,m_histidx_w).fill_main( w, nullptr, n );

        if ( m_opt.flags.has(TFlags::e) )
          vectAt(m_hists,m_histidx_e).fill_main( b.neutrons.ekin, w, n );

        if ( m_opt.flags.has(TFlags::mu) ) {
          for ( std::size_t i = 0; i < n; ++i )
            exit_angle[i] = std::acos( cos_exit_angle[i] );
          for ( std::size_t i = 0; i < n; ++i )
            exit_angle[i] *= kToDeg;
          vectAt(m_hists,m_histidx_mu).fill_main( exit_angle, w, n );
        }


        double lambda[basket_N];
        if ( m_opt.flags.has(TFlags::l) ) {
          //Calculate to 1/sqrt(ekin):
          for ( std::size_t i = 0; i < n; ++i )
            lambda[i] = 1.0 / ncmax(std::numeric_limits<double>::denorm_min(),
                                    b.neutrons.ekin[i]);
          for ( std::size_t i = 0; i < n; ++i )
            lambda[i] = std::sqrt( lambda[i] );
          //Finally, convert 1/sqrt(eV) to wavelength in Aa:
          constexpr double f2l = constexpr_ekin2wl(1.0);
          for ( std::size_t i = 0; i < n; ++i )
            lambda[i] *= f2l;
          vectAt(m_hists,m_histidx_l).fill_main( lambda, w, n );
        }

        double q[basket_N];//fixme avoid both lambda and q and ... temp vectors
                          //by reordering the code.!
        if ( m_opt.flags.has(TFlags::q) ) {
          //fixme: put in a helper method?
          //|qbar|^2 = |kfbar|^2+|kibar|^2-2kfbarDOTkibar
          //         = |kf|^2+|ki|^2-2*cosmu*|kf||ki|
          //         = ( Ef+Ei-2*cosmu*sqrt(Ef*Ei) ) * fact_ekin2ksq
          nc_assert( m_opt.beamEnergy.has_value() );
          const double Ei = m_opt.beamEnergy.value().get();
          for ( std::size_t i = 0; i < n; ++i )
            q[i] = b.neutrons.ekin[i] * Ei;
#ifndef NDEBUG
          for ( std::size_t i = 0; i < n; ++i )
            nc_assert( q[i]>=0.0 );
#endif
          for ( std::size_t i = 0; i < n; ++i )
            q[i] = std::sqrt( q[i] );
          for ( std::size_t i = 0; i < n; ++i )
            q[i] *= cos_exit_angle[i];
          for ( std::size_t i = 0; i < n; ++i )
            q[i] *= (-2.0);
          for ( std::size_t i = 0; i < n; ++i )
            q[i] += Ei;
          for ( std::size_t i = 0; i < n; ++i )
            q[i] += b.neutrons.ekin[i];
          constexpr double fact_ekin2ksq = ekin2ksq(1.0);
          for ( std::size_t i = 0; i < n; ++i )
            q[i] *= fact_ekin2ksq;
          //fixme: take absolute value? And what about units?
          vectAt(m_hists,m_histidx_q).fill_main( q, w, n );
        }

        if ( m_opt.flags.has(TFlags::nobreakdown) )
          return;

        //Detailed histograms:
        DetailedHistsID histid[basket_N];
        for ( std::size_t i = 0; i < n; ++i ) {
          auto nscat = b.cache.nscat[i];
          auto saw_inelas = b.cache.sawinelas[i];
          if (nscat==1) {
            histid[i] = ( saw_inelas
                          ? DetailedHistsID::SINGLESCAT_INELAS
                          : DetailedHistsID::SINGLESCAT_ELAS );
          } else if ( nscat>1 ) {
            histid[i] = ( saw_inelas
                          ? DetailedHistsID::MULTISCAT_OTHER
                          : DetailedHistsID::MULTISCAT_PUREELAS );
          } else {
            histid[i] = DetailedHistsID::NOSCAT;
          }
        }

        if ( m_opt.flags.has(TFlags::nscat) )
          vectAt(m_hists,m_histidx_nscat).fill_detail( histid,b.cache.nscat,
                                                       w, n );

        if ( m_opt.flags.has(TFlags::nscat_uw) )
          vectAt(m_hists,m_histidx_nscat_uw).fill_detail( histid,b.cache.nscat,
                                                          nullptr, n );

        if ( m_opt.flags.has(TFlags::cosmu) )
          vectAt(m_hists,m_histidx_cosmu).fill_detail( histid,
                                                       cos_exit_angle, w, n );

        if ( m_opt.flags.has(TFlags::e) )
          vectAt(m_hists,m_histidx_e).fill_detail( histid,
                                                   b.neutrons.ekin, w, n );

        if ( m_opt.flags.has(TFlags::l) )
          vectAt(m_hists,m_histidx_l).fill_detail( histid, lambda, w, n );

        if ( m_opt.flags.has(TFlags::q) )
          vectAt(m_hists,m_histidx_q).fill_detail( histid, q, w, n );

        if ( m_opt.flags.has(TFlags::w) )
          vectAt(m_hists,m_histidx_w).fill_detail( histid, w, nullptr, n );

        if ( m_opt.flags.has(TFlags::mu) )
          vectAt(m_hists,m_histidx_mu).fill_detail( histid, exit_angle, w, n );
      }

      shared_obj<TallyBase> cloneSetup() const override
      {
        auto tmp = makeSO<this_class_t>( m_opt );
        //NB: Cloning without histogram contents!
        return tmp;
      }

      void merge(TallyBase&& o_base) override
      {
        auto optr = dynamic_cast<this_class_t*>(&o_base);
        nc_assert_always(optr!=nullptr);
        const this_class_t& o = *optr;
        nc_assert_always( m_opt == o.m_opt );
        nc_assert_always( m_hists.size() == o.m_hists.size() );
        for ( auto i : ncrange(m_hists.size()) )
          m_hists.at(i).merge( o.m_hists.at(i) );
      };

      //FIXME: Rename as "toJSONItems", since it will be adding key,value pairs!
      VectS tallyItemNames() const override
      {
        VectS res;
        res.reserve( m_hists.size() );
        for ( auto& h : m_hists )
          res.emplace_back(h.main.title());
        return res;
      }

      void tallyItemToJSON( std::ostream& os, StrView itemName ) const override
      {
        for ( auto& h : m_hists ) {
          if ( itemName != h.main.title() )
            continue;
          os << ("{\"datatype\":\"NCrystalMiniMCTallyHistBreakdown_v1\","
                 "\"tallyname\":");
          streamJSON(os,h.main.title());
          os << ",\"total\":";
          h.main.toJSON(os);
          os << ",\"breakdown\":{";
          if ( !m_opt.flags.has(TFlags::nobreakdown) ) {
            bool first(true);
            for ( auto& hdetail : h.detailed ) {
              if (first)
                first=false;
              else
                os<<',';
              streamJSON(os,hdetail.title());
              os << ':';
              hdetail.toJSON(os);
            }
          }
          os << "} }";
        }
      }
    };
  }
}

//fixme: rename name of this file? StdHists?
//fixme: "Aborting plotting of empty histogram" -> we should just show it!

#endif
