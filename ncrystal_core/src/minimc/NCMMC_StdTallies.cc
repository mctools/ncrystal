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

#include "NCrystal/internal/minimc/NCMMC_StdTallies.hh"
#include "NCrystal/internal/minimc/NCMMC_Source.hh"
#include "NCrystal/internal/minimc/NCMMC_Utils.hh"
#include "NCMMC_BasketUtils.hh"
#include "NCrystal/internal/utils/NCVector.hh"

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;
namespace NCMMCT = NCrystal::MiniMC::Tallies;

NCMMCT::HistGroup::HistGroup( const Hists::Binning& b, const char * title)
  : main( b, title )
{
}

void NCMMCT::HistGroup::merge( const HistGroup& o )
{
  nc_assert( main.title() == o.main.title() );
  nc_assert( detailed.size() == o.detailed.size() );
  main.merge( o.main );
  for ( std::size_t i = 0; i < detailed.size(); ++i ) {
    nc_assert( detailed.at(i).title() == o.detailed.at(i).title() );
    detailed[i].merge( o.detailed[i] );
  }
}

NCMMCT::TallyNeutronInitialInfo
NCMMCT::TallyNeutronInitialInfo::create( const EngineOpts& eopts,
                                         const SourceMetaData& srcmd )
{
  TallyNeutronInitialInfo res;
  using TF = TallyFlags::Flags;
  constexpr TF::value_type flags_dir0 = TF::q|TF::mu|TF::theta;
  static_assert( flags_dir0 == 0x83, "" );//sanity check, value doesn't matter.
  constexpr auto flags_e0 = TF::q|TF::de;
  const bool needs_dir0 = eopts.tallyFlags.hasAny(flags_dir0);
  const bool needs_e0 = eopts.tallyFlags.hasAny(flags_e0);

  if ( !needs_dir0 && !needs_e0 )
    return res;

  nc_assert( eopts.tallyRef == EngineOpts::TallyReference::Truth
             || eopts.tallyRef == EngineOpts::TallyReference::Source );
  const bool tallyref_truth = ( eopts.tallyRef
                                == EngineOpts::TallyReference::Truth );
  auto& srcvar_dir0 = ( tallyref_truth
                        ? srcmd.fixedDirection : srcmd.meanDirection );
  auto& srcvar_e0 = ( tallyref_truth ? srcmd.fixedEnergy : srcmd.meanEnergy );

  const bool src_missing_dir0 = needs_dir0 && !srcvar_dir0.has_value();
  const bool src_missing_e0 = needs_e0 && !srcvar_e0.has_value();
  const bool src_missing_anything = (src_missing_dir0||src_missing_e0);

  if (src_missing_anything) {
    if ( tallyref_truth ) {
      //not an error in this case, but it does mean we must enable extended
      //baskets.
      res.needsExtendedBaskets = true;
    } else {
      nc_assert( eopts.tallyRef == EngineOpts::TallyReference::Source );
      //"tallyref=src", but source is missing something => error.
      auto giveerror = [&eopts](const char * quantitydescr,
                                TallyFlags flags)
      {
        nc_assert_always( ( TF::ALLHISTS & flags.getValue() )
                          == flags.getValue() );
        auto badflags = TallyFlags( eopts.tallyFlags.getValue()
                                    & flags.getValue() ).toStringList();
        nc_assert_always(!badflags.empty());
        NCRYSTAL_THROW2(BadInput,
                        "tallyref=src but tally \""<<badflags.front()
                        <<"\" needs a "<<quantitydescr
                        <<" which is not well defined with the chosen source");
      };
      if ( src_missing_dir0 )
        giveerror("direction",flags_dir0);
      if ( src_missing_e0 )
        giveerror("energy",flags_e0);
      nc_assert_always(false);
    }
  }

  //Cache quantities if relevant and available:
  if ( needs_dir0 )
    res.fixedDir = srcvar_dir0;
  if ( needs_e0 )
    res.fixedEnergy = srcvar_e0;
  return res;
}

NCMMCT::TallyStdHists_Options
NCMMCT::TallyStdHists_Options::create( const EngineOpts& eo,
                                       const SourceMetaData& srcmd )
{
  TallyStdHists_Options opt;
  opt.neutronInitialInfo = TallyNeutronInitialInfo::create( eo, srcmd );
  opt.flags = eo.tallyFlags;
  opt.tallyBinnings = eo.tallyBinnings;
  return opt;
}

NCMMCT::TallyStdHists_Data
NCMMCT::TallyStdHists_Data::create( const TallyStdHists_Options& opt )
{
  using TFlags = TallyFlags::Flags;
  NCMMCT::TallyStdHists_Data data;

  auto addhist = [&opt,&data]( TallyFlags::value_type f,
                               std::size_t& idxvar,
                               tally_hist_t::size_type nbins,
                               double xmin, double xmax,
                               bool fixnbins )
  {
    //sanity check that this is not already booked:
    nc_assert( idxvar == std::numeric_limits<std::size_t>::max());
    //Determine binning automatically, but allow user overrides:
    if (!fixnbins) {
      nc_assert_always( nbins%2 == 0);
      if ( opt.flags.has(TFlags::lowres) )
        nbins /= 2;
      if ( opt.flags.has(TFlags::highres) )
        nbins *= 10;
    }
    Hists::Binning b_auto( nbins, xmin, xmax);
    idxvar = data.hists.size();
    data.hists.emplace_back( opt.tallyBinnings.lookup( TallyFlags(f),
                                                       b_auto ),
                             TallyFlags::singleFlagToString(f) );
  };
  constexpr auto nbins_std = 200;
  constexpr auto nbins_angular = 180;
  if ( opt.flags.has(TFlags::theta) )
    addhist( TFlags::theta, data.histidx_theta, nbins_angular,
             0.0, 180.0, false );
  if ( opt.flags.has(TFlags::mu) )
    addhist( TFlags::mu, data.histidx_mu, nbins_std,
             -1.0, 1.0, false );
  if ( opt.flags.hasAny(TFlags::nscat|TFlags::nscat_uw) ) {
    unsigned nscatmax = ( opt.flags.has(TFlags::lowres) ? 10 : 20 );
    if (opt.flags.has(TFlags::highres))
      nscatmax += 20;
    if ( opt.flags.has(TFlags::nscat ) )
      addhist( TFlags::nscat, data.histidx_nscat,
               nscatmax+2,-1.5, nscatmax+0.5, true );
    if ( opt.flags.has(TFlags::nscat_uw ) )
      addhist( TFlags::nscat_uw, data.histidx_nscat_uw,
               nscatmax+2,-1.5, nscatmax+0.5, true );
  }
  if ( opt.flags.has(TFlags::w) )
    addhist( TFlags::w, data.histidx_w, nbins_std,
             0.0, 1.0, false );

  //fixme: better max of e/l/de/q tallys? emax = 3*max(beamE,kT)? Also, we
  //       should make sure we have unit tests for all tallies.
  if ( opt.flags.has(TFlags::e) )
    addhist( TFlags::e, data.histidx_e, nbins_std, 0.0, 1.0, false );
  if ( opt.flags.has(TFlags::l) )
    addhist( TFlags::l, data.histidx_l, nbins_std, 0.0, 10.0, false );
  if ( opt.flags.has(TFlags::de) )
    addhist( TFlags::de, data.histidx_de, nbins_std, -1.0, 1.0, false );
  if ( opt.flags.has(TFlags::q) )
    addhist( TFlags::q, data.histidx_q, nbins_std, 0.0, 10.0, false );
  if ( !opt.flags.has(TFlags::nobreakdown) ) {
    for ( auto& h : data.hists )
      h.bookDetailedHists();
  }


  return data;
}

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {
    namespace Tallies {
      namespace {

        template<class TBasketValBuf>
        void hgfill_main( HistGroup& hg,
                          const TBasketValBuf& x,
                          const BasketValBufDbl& w, std::size_t n)
        {
          for ( std::size_t i = 0; i < n; ++i)
            hg.main.fill( static_cast<double>(x.data[i]), w.data[i] );
        }

        template<class TBasketValBuf>
        void hgfill_detail( HistGroup& hg,
                            const DetailedHistsIDVect& hid,
                            const TBasketValBuf& x,
                            const BasketValBufDbl& w,
                            std::size_t n )
        {
          for ( std::size_t i = 0; i < n; ++i)
            vectAt(hg.detailed,static_cast<std::size_t>(hid.data[i]))
              .fill(static_cast<double>(x.data[i]),w.data[i]);
        }

        template<class TBasketValBuf>
        void hgfill_main( HistGroup& hg, const TBasketValBuf& x, std::size_t n)
        {
          for ( std::size_t i = 0; i < n; ++i)
            hg.main.fill( static_cast<double>(x.data[i]) );
        }

        template<class TBasketValBuf>
        void hgfill_detail( HistGroup& hg,
                          const DetailedHistsIDVect& hid,
                          const TBasketValBuf& x,
                          std::size_t n )
        {
          for ( std::size_t i = 0; i < n; ++i)
            vectAt(hg.detailed,static_cast<std::size_t>(hid.data[i]))
              .fill(static_cast<double>(x.data[i]));
        }

        void init_dethistidvect( DetailedHistsIDVect& v,
                                 const BasketValBufInt& nscat,
                                 const BasketValBufInt& nscat_inelas,
                                 std::size_t n )
        {
          auto& d = v.data;
          for ( std::size_t i = 0; i < n; ++i ) {
            if (nscat[i]>1) nclikely {
                d[i] = ( nscat_inelas[i]
                         ? DetailedHistsID::MULTISCAT_OTHER
                         : DetailedHistsID::MULTISCAT_PUREELAS );
              } else {
              if ( nscat[i] == 1 ) {
                d[i] = ( nscat_inelas[i]
                         ? DetailedHistsID::SINGLESCAT_INELAS
                         : DetailedHistsID::SINGLESCAT_ELAS );
              } else {
                d[i] = DetailedHistsID::NOSCAT;
              }
            }
          }
        }

        void tallyRecord_lambda( TallyStdHists_Data& data,
                                 const NeutronBasket& neutrons,
                                 const DetailedHistsIDVect * dethistid )
        {
          //Calculate to 1/sqrt(ekin) without triggering zero division errors:
          BasketValBufDbl tmp;
          const std::size_t n = neutrons.size();
          auto& nf = neutrons.fields;
          for ( std::size_t i = 0; i < n; ++i )
            tmp.data[i]
              = 1.0 / ncmax( std::numeric_limits<double>::denorm_min(),
                             nf.ekin[i] );
          for ( std::size_t i = 0; i < n; ++i )
            tmp.data[i] = std::sqrt( tmp[i] );
          //Convert 1/sqrt(eV) to wavelength in Aa:
          constexpr double f2l = constexpr_ekin2wl(1.0);
          for ( std::size_t i = 0; i < n; ++i )
            tmp.data[i] *= f2l;
          //Fill:
          auto& h = vectAt(data.hists,data.histidx_l);
          hgfill_main( h, tmp, nf.w, n );
          if ( dethistid )
            hgfill_detail( h, *dethistid, tmp, nf.w, n );
        }

        void tallyRecord_de( TallyStdHists_Data& data,
                             const Optional<NeutronEnergy>& fixedE0,
                             const NeutronBasket& b,
                             const DetailedHistsIDVect * dethistid,
                             const NeutronBasketFields* neutrons0 )
        {
          BasketValBufDbl de;
          const std::size_t n = b.size();
          if ( fixedE0.has_value() ) {
            const double Ei = fixedE0.value().get();
            for ( std::size_t i = 0; i < n; ++i )
              de.data[i] = Ei - b.fields.ekin[i];
          } else {
            nc_assert( neutrons0 != nullptr );
            for ( std::size_t i = 0; i < n; ++i )
              de.data[i] = neutrons0->ekin[i] - b.fields.ekin[i];
          }
          auto& h = vectAt(data.hists,data.histidx_de);
          hgfill_main( h, de, b.fields.w, n );
          if ( dethistid )
            hgfill_detail( h, *dethistid, de, b.fields.w, n );
        }

        void tallyRecord_q(  TallyStdHists_Data& data,
                             const Optional<NeutronEnergy>& fixedE0,
                             const BasketValBufDbl& mu,
                             const NeutronBasket& neutrons,
                             const DetailedHistsIDVect * dethistid,
                             const NeutronBasketFields* neutrons0 )
        {
          BasketValBufDbl q;
          //|qbar|^2 = |kfbar|^2+|kibar|^2-2kfbarDOTkibar
          //         = |kf|^2+|ki|^2-2*mu*|kf||ki|
          //         = ( Ef+Ei-2*mu*sqrt(Ef*Ei) ) * fact_ekin2ksq

          const std::size_t n = neutrons.size();
          if ( fixedE0.has_value() ) {
            const double Ei = fixedE0.value().get();
            for ( std::size_t i = 0; i < n; ++i )
              q.data[i] = neutrons.fields.ekin[i] * Ei;
          } else {
            nc_assert( neutrons0 != nullptr );
            for ( std::size_t i = 0; i < n; ++i )
              q.data[i] = neutrons.fields.ekin[i] * neutrons0->ekin[i];
          }
#ifndef NDEBUG
          for ( std::size_t i = 0; i < n; ++i )
            nc_assert( q.data[i]>=0.0 );
#endif
          for ( std::size_t i = 0; i < n; ++i )
            q.data[i] = std::sqrt( q.data[i] );
          for ( std::size_t i = 0; i < n; ++i )
            q.data[i] *= mu.data[i];
          for ( std::size_t i = 0; i < n; ++i )
            q.data[i] *= (-2.0);
          if ( fixedE0.has_value() ) {
            const double Ei = fixedE0.value().get();
            for ( std::size_t i = 0; i < n; ++i )
              q.data[i] += Ei;
          } else {
            nc_assert( neutrons0 != nullptr );
            for ( std::size_t i = 0; i < n; ++i )
              q.data[i] += neutrons0->ekin[i];
          }
          for ( std::size_t i = 0; i < n; ++i )
            q.data[i] += neutrons.fields.ekin[i];
          constexpr double fact_ekin2ksq = ekin2ksq(1.0);
          for ( std::size_t i = 0; i < n; ++i )
            q.data[i] *= fact_ekin2ksq;
          //From q^2 to |q|:
          for ( std::size_t i = 0; i < n; ++i )
            q.data[i] = Utils::fast_sqrt_clippos(q.data[i]);

          auto& h = vectAt(data.hists,data.histidx_q);
          hgfill_main( h, q, neutrons.fields.w, n );
          if ( dethistid )
            hgfill_detail( h, *dethistid, q, neutrons.fields.w, n );
        }

        void tallyRecord_mu_theta_q( TallyStdHists_Data& data,
                                     const TallyStdHists_Options& opt,
                                     const NeutronBasket& b,
                                     const DetailedHistsIDVect * dethistid,
                                     const NeutronBasketFields* neutrons0 )
        {
          using TFlags = TallyFlags::Flags;
          nc_assert( opt.flags.has(TFlags::mu)
                     || opt.flags.has(TFlags::theta)
                     || opt.flags.has(TFlags::q) );

          //////////////////
          //First find mu
          const std::size_t n = b.size();
          BasketValBufDbl mu;
          auto& bf = b.fields;
          if ( opt.neutronInitialInfo.fixedDir.has_value() ) {
            const double bx = opt.neutronInitialInfo.fixedDir.value()[0];
            const double by = opt.neutronInitialInfo.fixedDir.value()[1];
            const double bz = opt.neutronInitialInfo.fixedDir.value()[2];
            if ( bz==1 ) {
              //using (0,0,1) as reference is such a normal case that it makes
              //sense to take advantage of the fact that mu=uz in this case:
              BasketUtils::memcpydata( mu.data, bf.uz.data, n );
            } else {
              for ( std::size_t i = 0; i < n; ++i )
                mu.data[i] = bx * bf.ux[i];
              for ( std::size_t i = 0; i < n; ++i )
                mu.data[i] += by * bf.uy[i];
              for ( std::size_t i = 0; i < n; ++i )
                mu.data[i] += bz * bf.uz[i];
            }
          } else {
            //We must take the individual initial value of neutrons (so
            //presumably we are using extended baskets).
            nc_assert( neutrons0 != nullptr );
            for ( std::size_t i = 0; i < n; ++i )
              mu.data[i] = neutrons0->ux[i] * bf.ux[i];
            for ( std::size_t i = 0; i < n; ++i )
              mu.data[i] += neutrons0->uy[i] * bf.uy[i];
            for ( std::size_t i = 0; i < n; ++i )
              mu.data[i] += neutrons0->uz[i] * bf.uz[i];
          }

#ifndef NDEBUG
          for ( std::size_t i = 0; i < n; ++i ) {
            nc_assert_always( mu.data[i] > -(1.0+1e-14) );
            nc_assert_always( mu.data[i] <  (1.0+1e-14) );
          }
#endif
          for ( std::size_t i = 0; i < n; ++i )
            mu.data[i] = ncclamp(mu.data[i],-1.0,1.0);

          /////////////
          //Tally mu
          if ( opt.flags.has(TFlags::mu) ) {
            auto& h = vectAt(data.hists,data.histidx_mu);
            hgfill_main( h, mu, bf.w, n );
            if ( dethistid )
              hgfill_detail( h, *dethistid, mu, bf.w, n );
          }

          //////////
          //Tally q
          if ( opt.flags.has(TFlags::q) )
            tallyRecord_q( data, opt.neutronInitialInfo.fixedEnergy,
                           mu, b, dethistid, neutrons0 );

          ////////////////////////////////////////////////
          //Tally theta (using mu array to hold theta values)
          if ( opt.flags.has(TFlags::theta) ) {
            for ( std::size_t i = 0; i < n; ++i )
              mu.data[i] = std::acos( mu.data[i] );
            for ( std::size_t i = 0; i < n; ++i )
              mu.data[i] *= kToDeg;
            auto& h = vectAt(data.hists,data.histidx_theta);
            hgfill_main( h, mu, bf.w, n );
            if ( dethistid )
              hgfill_detail( h, *dethistid, mu, bf.w, n );
          }
        }
      }
    }
  }
}

void NCMMCT::HistGroup::bookDetailedHists()
{
  static_assert( NDETAIL == 5, "" );
  static_assert( static_cast<std::size_t>(DetailedHistsID::END) == 5, "" );
  nc_assert( detailed.empty() );
  auto binning = main.binning();
  detailed.emplace_back(binning,"NOSCAT");
  detailed.emplace_back(binning,"SINGLESCAT_ELAS");
  detailed.emplace_back(binning,"SINGLESCAT_INELAS");
  detailed.emplace_back(binning,"MULTISCAT_PUREELAS");
  detailed.emplace_back(binning,"MULTISCAT_OTHER");
}

void NCMMCT::TallyStdHists_Data::merge(const TallyStdHists_Data& o)
{
  const auto n = hists.size();
  nc_assert_always( o.hists.size() == n
                    && histidx_theta == o.histidx_theta
                    && histidx_mu == o.histidx_mu
                    && histidx_nscat == o.histidx_nscat
                    && histidx_w == o.histidx_w
                    && histidx_e == o.histidx_e
                    && histidx_de == o.histidx_de
                    && histidx_l == o.histidx_l
                    && histidx_nscat_uw == o.histidx_nscat_uw
                    && histidx_q == o.histidx_q );
  for ( auto i : ncrange(n) )
    hists.at(i).merge( o.hists.at(i) );
}

NC::VectS NCMMCT::TallyStdHists_Data::titles() const
{
  VectS res;
  res.reserve( this->hists.size() );
  for ( auto& h : this->hists )
    res.emplace_back(h.main.title());
  return res;
}

void NCMMCT::TallyStdHists_Data::histToJSONFindByTitle( std::ostream& os,
                                                        StrView title ) const
{
  for ( auto& h : this->hists ) {
    if ( title != h.main.title() )
      continue;
    os << ("{\"datatype\":\"NCrystalMiniMCTallyHistBreakdown_v1\","
           "\"tallyname\":");
    streamJSON(os,h.main.title());
    os << ",\"total\":";
    h.main.toJSON(os);
    os << ",\"breakdown\":{";
    if ( !h.detailed.empty() ) {
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

void NCMMCT::tallyRecord( TallyStdHists_Data& data,
                          const TallyStdHists_Options& opt,
                          const NeutronBasket& neutrons,
                          const BasketValBufInt& nscat,
                          const BasketValBufInt& nscat_inelas,
                          const NeutronBasketFields* neutrons0 )
{
  using TFlags = TallyFlags::Flags;
  const std::size_t n = neutrons.size();
  nc_assert( n <= basket_N );
  if ( n == 0 )
    return;
  //////////////////////////////////
  //Find DetailedHistID if required:
  DetailedHistsIDVect histid_buf;
  const DetailedHistsIDVect * histid = nullptr;
  if ( !opt.flags.has(TFlags::nobreakdown) ) {
    init_dethistidvect( histid_buf, nscat, nscat_inelas, n );
    histid = &histid_buf;
  }

  //////////////////////////////
  //Fill cos(mu), theta, q tallies:
  constexpr auto flag_miscmu = (TFlags::mu|TFlags::theta|TFlags::q );
  if ( opt.flags.hasAny(flag_miscmu) )
    tallyRecord_mu_theta_q( data, opt, neutrons, histid, neutrons0 );

  constexpr auto flag_nonmiscmu = ( TFlags::ALLHISTS & (~flag_miscmu) );
  if ( !opt.flags.hasAny(flag_nonmiscmu) )
    return;

  ///////////////
  //Fill e tally:
  if ( opt.flags.has(TFlags::e) ) {
    auto& h = vectAt(data.hists,data.histidx_e);
    hgfill_main( h, neutrons.fields.ekin, neutrons.fields.w, n );
    if ( histid )
      hgfill_detail( h, *histid, neutrons.fields.ekin, neutrons.fields.w, n );
  }

  ////////////////
  //Fill de tally:
  if ( opt.flags.has(TFlags::de) )
    tallyRecord_de( data, opt.neutronInitialInfo.fixedEnergy,
                    neutrons, histid, neutrons0 );

  ////////////////////
  //Fill lambda tally:
  if ( opt.flags.has(TFlags::l) )
    tallyRecord_lambda( data, neutrons, histid );

  constexpr auto flag_other=(TFlags::nscat|TFlags::nscat_uw|TFlags::w);
  if ( !opt.flags.hasAny(flag_other) )
    return;

  using nscat_t_cv = typename std::remove_all_extents<decltype(nscat.data)>::type;
  using nscat_t = typename std::remove_cv<nscat_t_cv>::type;
  static_assert( std::is_same<nscat_t,int>::value, "");

  if ( opt.flags.has(TFlags::nscat) ) {
    auto& h = vectAt(data.hists,data.histidx_nscat);
    hgfill_main( h, nscat, neutrons.fields.w, n );
    if ( histid )
      hgfill_detail( h, *histid, nscat, neutrons.fields.w, n );
  }

  if ( opt.flags.has(TFlags::nscat_uw) ) {
    auto& h = vectAt(data.hists,data.histidx_nscat_uw);
    hgfill_main( h, nscat, n );
    if ( histid )
      hgfill_detail( h, *histid, nscat, n );
  }

  if ( opt.flags.has(TFlags::w) ) {
    auto& h = vectAt(data.hists,data.histidx_w);
    hgfill_main( h, neutrons.fields.w, n );
    if ( histid )
      hgfill_detail( h, *histid, neutrons.fields.w, n );
  }
}

void NCMMCT::TallyStdHists::merge(Tally&& o_base)
{
  auto optr = dynamic_cast<this_class_t*>(&o_base);
  nc_assert_always(optr!=nullptr);
  const this_class_t& o = *optr;
  nc_assert_always( m_opt == o.m_opt );
  m_data.merge( o.m_data );
}

NCMMCT::TallyStdHists::TallyStdHists( const EngineOpts& eo,
                                      const SourceMetaData& srcmd )
  : TallyStdHists( private_constructor_t{},
                   Options::create(eo,srcmd) )
{
}

void NCMMCT::TallyStdHists::registerResults( const Basket& b)
{
  nc_assert( b.neutrons );
  nc_assert( b.nscat );
  nc_assert( b.nscat_inelas );
  nc_assert( !( m_opt.neutronInitialInfo.needsExtendedBaskets
                && b.neutrons_initial==nullptr) );

  tallyRecord( m_data, m_opt,
               *b.neutrons, *b.nscat, *b.nscat_inelas, b.neutrons_initial );
}

NCMMCT::TallyStdHists::TallyStdHists( private_constructor_t, Options opt )
  : m_opt(opt), m_data(Data::create(opt))
{
}

bool NCMMCT::TallyStdHists::needsExtendedBaskets() const
{
  return m_opt.neutronInitialInfo.needsExtendedBaskets;
}

NC::shared_obj<NCMMC::Tally> NCMMCT::TallyStdHists::cloneSetup() const
{
  //NB: Cloning without histogram contents!
  return makeSO<this_class_t>( private_constructor_t{}, m_opt );
}

NC::VectS NCMMCT::TallyStdHists::tallyItemNames() const
{
  return m_data.titles();
}

void NCMMCT::TallyStdHists::tallyItemToJSON( std::ostream& os,
                                             StrView itemName ) const
{
  m_data.histToJSONFindByTitle( os, itemName );
}
