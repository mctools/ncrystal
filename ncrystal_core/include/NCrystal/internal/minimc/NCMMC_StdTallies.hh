#ifndef NCrystal_MMC_StdTallies_hh
#define NCrystal_MMC_StdTallies_hh

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
#include "NCrystal/internal/minimc/NCMMC_Tally.hh"
#include "NCrystal/internal/utils/NCHists.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    struct Tally_ExitAngle_Options
    {
      unsigned nbins = 1800;//10 per degree
      unsigned detail_level = 0;
      // 0: Just the 1D exit angle hist
      // 1: Also the running stats for that hist
      // 2: Also a bunch of other hists, depending on neutron history.
    };

    enum class TallyCollectRunningStats { YES, NO };

    template<class TBasket>
    class Tally_ExitAngle final : public Tally<TBasket> {
    public:
      using Options = Tally_ExitAngle_Options;
      using hist_exitangle_t = Hists::HistBinData1D<
        Hists::AllowWeights::YES,
        Hists::OverflowHandling::Clamp,
        SmallVector_IC<double,1800> >;

    private:
      using this_class_t = Tally_ExitAngle;
      hist_exitangle_t m_exitangle_binned;
      Hists::RunningStats1D m_exitangle_stats;
      Options m_opt;
      using extrahist_exitangle_t = Hists::Hist1D<hist_exitangle_t::opt_allow_weights,
                                                  hist_exitangle_t::opt_of_handling,
                                                  std::vector<double>>;
      std::vector<extrahist_exitangle_t> m_detailed_hists;
      enum class DetailedHistsID : std::size_t {
        NOSCAT = 0,
        SINGLESCAT_ELAS,
        SINGLESCAT_INELAS,
        MULTISCAT_PUREELAS,
        MULTISCAT_OTHER,
        END
      };

      static std::size_t histid2idx( DetailedHistsID val )
      {
        return static_cast<std::size_t>(val);
      }

    public:

      Tally_ExitAngle( Options opt = {} )
        : m_exitangle_binned( opt.nbins, 0.0, 180.0 ),
          m_opt( std::move(opt) )
      {

        if ( m_opt.detail_level >= 2 ) {
          auto adddh = [this]( DetailedHistsID histid, std::string title )
          {
            nc_assert_always( static_cast<std::size_t>(histid) == m_detailed_hists.size() );
            m_detailed_hists.emplace_back( m_opt.nbins, 0.0, 180.0, std::move(title) );
          };
          adddh( DetailedHistsID::NOSCAT, "NOSCAT" );
          adddh( DetailedHistsID::SINGLESCAT_ELAS, "SINGLESCAT_ELAS" );
          adddh( DetailedHistsID::SINGLESCAT_INELAS, "SINGLESCAT_INELAS" );
          adddh( DetailedHistsID::MULTISCAT_PUREELAS, "MULTISCAT_PUREELAS" );
          adddh( DetailedHistsID::MULTISCAT_OTHER, "MULTISCAT_OTHER" );
          nc_assert_always( histid2idx(DetailedHistsID::END) == m_detailed_hists.size() );
        }
      }

      using basket_t = TBasket;
      void registerResults( const basket_t& b ) override
      {
        const std::size_t n = b.size();
        NCRYSTAL_DEBUGMMCMSG("Got result basket with size "<<n);
        double exit_angle[basket_N];
        for ( std::size_t i = 0; i < n; ++i ) {
          nc_assert( b.neutrons.uz[i] > -(1.0+1e-14) );
          nc_assert( b.neutrons.uz[i] <  (1.0+1e-14) );
          exit_angle[i] = std::acos( ncclamp(b.neutrons.uz[i],-1.0,1.0) );
        }
        for ( std::size_t i = 0; i < n; ++i )
          exit_angle[i] *= kToDeg;
        for ( std::size_t i = 0; i < n; ++i )
          m_exitangle_binned.fill( exit_angle[i], b.neutrons.w[i] );

        if ( hasRunningStats() ) {
          for ( std::size_t i = 0; i < n; ++i )
            m_exitangle_stats.registerValue( exit_angle[i], b.neutrons.w[i] );
        }
        if ( m_opt.detail_level >= 2 ) {
          for ( std::size_t i = 0; i < n; ++i ) {
            auto nscat = b.cache.nscat[i];
            auto saw_inelas = b.cache.sawinelas[i];
            auto histid = DetailedHistsID::NOSCAT;
            if (nscat==1) {
              histid = ( saw_inelas
                         ? DetailedHistsID::SINGLESCAT_INELAS
                         : DetailedHistsID::SINGLESCAT_ELAS );
            } else if ( nscat>1 ) {
              histid = ( saw_inelas
                         ? DetailedHistsID::MULTISCAT_OTHER
                         : DetailedHistsID::MULTISCAT_PUREELAS );
            }
            vectAt(m_detailed_hists,histid2idx(histid)).fill(exit_angle[i],
                                                             b.neutrons.w[i]);
          }
        }
      }

      shared_obj<TallyBase> clone() const override
      {
        //NB: Cloning without histogram contents!
        return makeSO<this_class_t>( m_opt );
      }

      void merge(TallyBase&& o_base) override
      {
        auto optr = dynamic_cast<this_class_t*>(&o_base);
        nc_assert_always(optr!=nullptr);
        const this_class_t& o = *optr;
        m_exitangle_binned.merge( o.m_exitangle_binned );
        if ( hasRunningStats() )
          m_exitangle_stats.merge( o.m_exitangle_stats );
        nc_assert( m_detailed_hists.size() == o.m_detailed_hists.size() );
        for ( auto i : ncrange(m_detailed_hists.size()) )
          m_detailed_hists.at(i).merge( o.m_detailed_hists.at(i) );
      }

      const hist_exitangle_t& getExitAngleBinned() const
      {
        return m_exitangle_binned;
      }

      bool hasRunningStats() const noexcept { return m_opt.detail_level>=1; }
      bool hasMultiHists() const noexcept { return m_opt.detail_level>=2; }
      bool hasJSON() const noexcept { return hasRunningStats() || hasMultiHists(); }
      void toJSON( std::ostringstream& os ) const
      {
        if (!hasJSON()) {
          os << "{}";
          return;
        }
        os << '{';
        streamJSON(os,"main_stats");
        os << ':';
        this->getExitAngleStats().toJSON(os);
        os << ',';
        streamJSON(os,"breakdown_hists");
        if (!hasMultiHists()) {
          os << ": [] }";
          return;
        }
        os << ": [";
        bool first(true);
        for ( auto& h : m_detailed_hists ) {
          if (first) {
            first=false;
          } else {
            os<<',';
          }
          h.toJSON(os);
        }
        os << "] }";
      }

      const Hists::RunningStats1D& getExitAngleStats() const
      {
        nc_assert_always(hasRunningStats());
        return m_exitangle_stats;
      }

    };

  }
}

#endif
