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
#include "NCrystal/internal/minimc/NCMMC_Basket.hh"
#include "NCrystal/internal/minimc/NCMMC_EngineOpts.hh"
#include "NCrystal/internal/utils/NCHists.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    class Source;

    namespace TallyStdHistsImpl {

      using tally_hist_t =
        Hists::Hist1D<Hists::AllowWeights::YES,
                      Hists::OverflowHandling::Record,
                      std::vector<double>>;

      enum class DetailedHistsID : std::size_t {
        NOSCAT = 0,
        SINGLESCAT_ELAS,
        SINGLESCAT_INELAS,
        MULTISCAT_PUREELAS,
        MULTISCAT_OTHER,
        END
      };

      struct DetailedHistsIDVect final {
        DetailedHistsID data[basket_N];
      };

      struct HistGroup final : public MoveOnly {
        //Each quantity tallied is kept in a group of a main histogram (for all
        //neutrons), and (if enabled) also 5 histograms for contribution
        //breakdown by components.
        static constexpr std::size_t NDETAIL = 5;

        tally_hist_t main;
        SmallVector<tally_hist_t,NDETAIL> detailed;

        HistGroup( const Hists::Binning&, const char * title);
        void bookDetailedHists();
        void merge( const HistGroup& );
      };

      struct TallyStdHists_Options {
        TallyFlags flags;
        Optional<NeutronDirection> beamDir;
        Optional<NeutronEnergy> beamEnergy;
        TallyBinningOverrides tallyBinnings;
        bool operator==( const TallyStdHists_Options& o ) const {
          return ( flags.getValue() == o.flags.getValue()
                   && beamDir == o.beamDir
                   && beamEnergy == o.beamEnergy
                   && tallyBinnings == o.tallyBinnings );
        }
        static TallyStdHists_Options create( const EngineOpts& eo,
                                             const Source& src);
      };

      struct TallyStdHists_Data {
        //Histograms, including a cache of indices for efficient later access:
        std::vector<HistGroup> hists;
        std::size_t histidx_mu = std::numeric_limits<std::size_t>::max();
        std::size_t histidx_cosmu = std::numeric_limits<std::size_t>::max();
        std::size_t histidx_nscat = std::numeric_limits<std::size_t>::max();
        std::size_t histidx_w = std::numeric_limits<std::size_t>::max();
        std::size_t histidx_e = std::numeric_limits<std::size_t>::max();
        std::size_t histidx_l = std::numeric_limits<std::size_t>::max();
        std::size_t histidx_nscat_uw = std::numeric_limits<std::size_t>::max();
        std::size_t histidx_q = std::numeric_limits<std::size_t>::max();
        static TallyStdHists_Data create( const TallyStdHists_Options& );
        const tally_hist_t& accessHistogram( StrView histname,
                                                  Optional<DetailedHistsID>
                                                  detailid = NullOpt ) const;
        VectS titles() const;
        void histToJSONFindByTitle( std::ostream&, StrView title ) const;
        void merge(const TallyStdHists_Data&);
      };

      void tallyRecord( TallyStdHists_Data&,
                        const TallyStdHists_Options&,
                        const NeutronBasket&,
                        const BasketValBufInt&,
                        const BasketValBufBool& );

      template<class TBasket>
      class TallyStdHists final : public Tally<TBasket> {
      public:
        using Options = TallyStdHists_Options;
        using Data = TallyStdHists_Data;
      private:
        using this_class_t = TallyStdHists;
        using TFlags = TallyFlags::Flags;
        //Do not add other data members than these two:
        Options m_opt;
        Data m_data;
        //Private struct, to prevent accidental usage of a constructor which
        //should really be private, but can not be since we use makeSO<..>:
        struct private_constructor_t {};
      public:

        //Initialise from engine and source options (could in principle also
        //depend on geometry options):
        TallyStdHists( const EngineOpts& eo, const Source& src)
          : TallyStdHists( private_constructor_t{},
                           Options::create(eo,src) )
        {
        }

        using basket_t = TBasket;

        void registerResults( const basket_t& b ) override
        {
          tallyRecord( m_data, m_opt,
                       b.neutrons, b.cache.nscat, b.cache.sawinelas );
        }

        //Don't use this constructor from outside the class:
        TallyStdHists( private_constructor_t, Options opt )
          : m_opt(opt), m_data(Data::create(opt))
        {
        }

        shared_obj<TallyBase> cloneSetup() const override
        {
          //NB: Cloning without histogram contents!
          return makeSO<this_class_t>( private_constructor_t{}, m_opt );
        }

        void merge(TallyBase&& o_base) override
        {
          auto optr = dynamic_cast<this_class_t*>(&o_base);
          nc_assert_always(optr!=nullptr);
          const this_class_t& o = *optr;
          nc_assert_always( m_opt == o.m_opt );
          m_data.merge( o.m_data );
        };

        const tally_hist_t& accessHistogram( StrView histname,
                                             Optional<DetailedHistsID>
                                             detailid = NullOpt ) const
        {
          return m_data.accessHistogram(histname,detailid);
        }

        VectS tallyItemNames() const override
        {
          return m_data.titles();
        }

        void tallyItemToJSON( std::ostream& os, StrView itemName ) const override
        {
          m_data.histToJSONFindByTitle( os, itemName );
        }
      };
    }

    template<class TBasket>
    using TallyStdHists = TallyStdHistsImpl::TallyStdHists<TBasket>;
  }
}
//fixme: "Aborting plotting of empty histogram" -> we should just show it!

#endif
