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

#include "NCrystal/internal/minimc/NCMMC_Tally.hh"
#include "NCrystal/internal/minimc/NCMMC_EngineOpts.hh"
#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/internal/utils/NCHists.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    struct SourceMetaData;

    namespace Tallies {

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

      struct TallyNeutronInitialInfo {
        //Info related to quantities needing knowledge of incoming parameters of
        //neutron. Depends on which tallies are enabled and whether or not
        //source provides particles with well defined fixed or mean values.
        Optional<NeutronDirection> fixedDir;
        Optional<NeutronEnergy> fixedEnergy;
        bool needsExtendedBaskets = false;
        static TallyNeutronInitialInfo create( const EngineOpts&,
                                               const SourceMetaData& );
        bool operator==( const TallyNeutronInitialInfo& o ) const {
          return ( fixedDir == o.fixedDir
                   && fixedEnergy == o.fixedEnergy
                   && needsExtendedBaskets == o.needsExtendedBaskets );
        }
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
        TallyNeutronInitialInfo neutronInitialInfo;
        TallyBinningOverrides tallyBinnings;
        bool operator==( const TallyStdHists_Options& o ) const {
          return ( flags.getValue() == o.flags.getValue()
                   && neutronInitialInfo == o.neutronInitialInfo
                   && tallyBinnings == o.tallyBinnings );
        }
        static TallyStdHists_Options create( const EngineOpts& eo,
                                             const SourceMetaData& src);
      };

      struct TallyStdHists_Data {
        //Histograms, including a cache of indices for efficient later access:
        std::vector<HistGroup> hists;
        std::size_t histidx_theta = std::numeric_limits<std::size_t>::max();
        std::size_t histidx_mu = std::numeric_limits<std::size_t>::max();
        std::size_t histidx_nscat = std::numeric_limits<std::size_t>::max();
        std::size_t histidx_w = std::numeric_limits<std::size_t>::max();
        std::size_t histidx_e = std::numeric_limits<std::size_t>::max();
        std::size_t histidx_l = std::numeric_limits<std::size_t>::max();
        std::size_t histidx_nscat_uw = std::numeric_limits<std::size_t>::max();
        std::size_t histidx_q = std::numeric_limits<std::size_t>::max();
        std::size_t histidx_de = std::numeric_limits<std::size_t>::max();
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
                        const BasketValBufInt&,
                        const NeutronBasketFields* );

      class TallyStdHists final : public TallyBase {
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
        TallyStdHists( const EngineOpts&, const SourceMetaData&);
        bool needsExtendedBaskets() const override;
        void registerResultsUB( const UniversalBasket&) override;
        shared_obj<TallyBase> cloneSetup() const override;
        void merge(TallyBase&&) override;
        const tally_hist_t& accessHistogram( StrView histname,
                                             Optional<DetailedHistsID>
                                             detailid = NullOpt ) const;
        VectS tallyItemNames() const override;
        void tallyItemToJSON( std::ostream&, StrView itemName ) const override;

        //Don't use this constructor from outside the class:
        TallyStdHists( private_constructor_t, Options );
      };
    }
    using TallyStdHists = Tallies::TallyStdHists;
  }
}

#endif
