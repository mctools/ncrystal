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
#include "NCrystal/internal/utils/NCVector.hh"

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;
namespace NCMMCT = NCrystal::MiniMC::TallyStdHistsImpl;

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

NCMMCT::TallyStdHists_Options
NCMMCT::TallyStdHists_Options::create( const EngineOpts& eo,
                                       const Source& src)
{
  using TFlags = TallyFlags::Flags;
  TallyStdHists_Options opt;
  //Tally options from EngineOpts:
  opt.flags = eo.tallyFlags;
  opt.tallyBinnings = eo.tallyBinnings;
  //Nominal beam dir + energy from source:
  Optional<NeutronDirection> optbd = eo.tallyBeamDir;
  if ( !optbd.has_value() )
    optbd = src.nominalBeamDirection();
  if ( optbd.has_value() ) {
    //Note: We only store beamdir if it is not (0,0,1), since it makes
    //it cheaper to calculate cos(theta) for this common case.
    auto v = optbd.value().as<Vector>().unit();
    if ( v.z() != 1.0 )
      opt.beamDir = v.as<NeutronDirection>();
  }
  opt.beamEnergy = eo.tallyBeamEnergy;
  if ( !opt.beamEnergy.has_value() )
    opt.beamEnergy = src.nominalBeamEnergy();

  if ( !opt.beamEnergy.has_value() ) {
    //fixme: actually test this branch (maybe allow
    //       enginecfg.beamenergy=none to remove actual src.beamEnergy?)!
    constexpr TallyFlags::value_type hists_needing_Ei = TFlags::q;
    auto bf = opt.flags.getValue() & hists_needing_Ei;
    if ( bf ) {
      std::string bfl(joinstr(TallyFlags(bf).toStringList(),"\", \""));
      NCRYSTAL_WARN("Disabling some tallies due to missing knowledge of "
                    "incoming energy: \""<<bfl<<"\" (consider provi"
                    "ding one with the enginecfg beamenergy parameter).");
      opt.flags.remove( bf );
    }
  }
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
  //Fixme for the TFlags::w histogram (i.e. non-weighted hist), the
  //fractions shown in the breakdown plots are very misleading.
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
             0.0, 1.0, false );//fixme: logw? else upper limit should be
  //related to initial src weight if known
  //(maybe we should for now simply disallow
  //changing initial weights).

  //fixme: better max of e/l/q tallys? emax = 3*max(beamE,kT)?
  if ( opt.flags.has(TFlags::e) )
    addhist( TFlags::e, data.histidx_e, nbins_std, 0.0, 1.0, false );
  if ( opt.flags.has(TFlags::l) )
    addhist( TFlags::l, data.histidx_l, nbins_std, 0.0, 10.0, false );
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
    namespace TallyStdHistsImpl {
      namespace {

        template<class TBasketValBuf>
        void hgfill_main( HistGroup& hg,
                          const TBasketValBuf& x,
                          const BasketValBufDbl& w, std::size_t n)
        {
          for ( std::size_t i = 0; i < n; ++i)
            hg.main.fill( static_cast<double>(x[i]), w[i] );
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
              .fill(static_cast<double>(x[i]),w[i]);
        }

        template<class TBasketValBuf>
        void hgfill_main( HistGroup& hg, const TBasketValBuf& x, std::size_t n)
        {
          for ( std::size_t i = 0; i < n; ++i)
            hg.main.fill( static_cast<double>(x[i]) );
        }

        template<class TBasketValBuf>
        void hgfill_detail( HistGroup& hg,
                          const DetailedHistsIDVect& hid,
                          const TBasketValBuf& x,
                          std::size_t n )
        {
          for ( std::size_t i = 0; i < n; ++i)
            vectAt(hg.detailed,static_cast<std::size_t>(hid.data[i]))
              .fill(static_cast<double>(x[i]));
        }

        void init_dethistidvect( DetailedHistsIDVect& v,
                                 const BasketValBufInt& nscat,
                                 const BasketValBufBool& sawinelas,
                                 std::size_t n )
        {
          auto& d = v.data;
          for ( std::size_t i = 0; i < n; ++i ) {
            if (nscat[i]>1) nclikely {
                d[i] = ( sawinelas[i]
                         ? DetailedHistsID::MULTISCAT_OTHER
                         : DetailedHistsID::MULTISCAT_PUREELAS );
              } else {
              if ( nscat[i] == 1 ) {
                d[i] = ( sawinelas[i]
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
          for ( std::size_t i = 0; i < n; ++i )
            tmp[i]
              = 1.0 / ncmax( std::numeric_limits<double>::denorm_min(),
                             neutrons.ekin[i] );
          for ( std::size_t i = 0; i < n; ++i )
            tmp[i] = std::sqrt( tmp[i] );
          //Convert 1/sqrt(eV) to wavelength in Aa:
          constexpr double f2l = constexpr_ekin2wl(1.0);
          for ( std::size_t i = 0; i < n; ++i )
            tmp[i] *= f2l;
          //Fill:
          auto& h = vectAt(data.hists,data.histidx_l);
          hgfill_main( h, tmp, neutrons.w, n );
          if ( dethistid )
            hgfill_detail( h, *dethistid, tmp, neutrons.w, n );
        }

        void tallyRecord_q(  TallyStdHists_Data& data,
                             const NeutronEnergy& beamEnergy,
                             const BasketValBufDbl& mu,
                             const NeutronBasket& neutrons,
                             const DetailedHistsIDVect * dethistid )
        {
          BasketValBufDbl q;
          //|qbar|^2 = |kfbar|^2+|kibar|^2-2kfbarDOTkibar
          //         = |kf|^2+|ki|^2-2*mu*|kf||ki|
          //         = ( Ef+Ei-2*mu*sqrt(Ef*Ei) ) * fact_ekin2ksq
          const double Ei = beamEnergy.get();
          const std::size_t n = neutrons.size();
          for ( std::size_t i = 0; i < n; ++i )
            q.data[i] = neutrons.ekin[i] * Ei;
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
          for ( std::size_t i = 0; i < n; ++i )
            q.data[i] += Ei;
          for ( std::size_t i = 0; i < n; ++i )
            q.data[i] += neutrons.ekin[i];
          constexpr double fact_ekin2ksq = ekin2ksq(1.0);
          for ( std::size_t i = 0; i < n; ++i )
            q.data[i] *= fact_ekin2ksq;
          //From q^2 to |q|:
          for ( std::size_t i = 0; i < n; ++i )
            q.data[i] = Utils::fast_sqrt_clippos(q.data[i]);

          auto& h = vectAt(data.hists,data.histidx_q);
          hgfill_main( h, q, neutrons.w, n );
          if ( dethistid )
            hgfill_detail( h, *dethistid, q, neutrons.w, n );
        }

        void tallyRecord_mu_theta_q( TallyStdHists_Data& data,
                                     const TallyStdHists_Options& opt,
                                     const NeutronBasket& b,
                                     const DetailedHistsIDVect * dethistid )
        {
          using TFlags = TallyFlags::Flags;
          nc_assert( opt.flags.has(TFlags::mu)
                     || opt.flags.has(TFlags::theta)
                     || opt.flags.has(TFlags::q) );

          //////////////////
          //First find mu
          const std::size_t n = b.size();
          BasketValBufDbl mu;
          //fixme: if beamDIR has a value, we could look for (0,0,1) as a
          //special case. This can be done in the constructor.
          if ( !opt.beamDir.has_value() ) {
            //use (0,0,1) as reference, so just copy over uz values.
            detail::memcpydata( mu.data, b.uz.data, n );
          } else {
            const double bx( opt.beamDir.value()[0] );
            const double by( opt.beamDir.value()[1] );
            const double bz( opt.beamDir.value()[2] );
            for ( std::size_t i = 0; i < n; ++i )
              mu.data[i] = bx * b.ux[i];
            for ( std::size_t i = 0; i < n; ++i )
              mu.data[i] += by * b.uy[i];
            for ( std::size_t i = 0; i < n; ++i )
              mu.data[i] += bz * b.uz[i];
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
            hgfill_main( h, mu, b.w, n );
            if ( dethistid )
              hgfill_detail( h, *dethistid, mu, b.w, n );
          }

          //////////
          //Tally q
          if ( opt.flags.has(TFlags::q) ) {
            nc_assert( opt.beamEnergy.has_value() );
            tallyRecord_q( data, opt.beamEnergy.value(), mu, b, dethistid );
          }

          ////////////////////////////////////////////////
          //Tally theta (using mu array to hold theta values)
          if ( opt.flags.has(TFlags::theta) ) {
            for ( std::size_t i = 0; i < n; ++i )
              mu.data[i] = std::acos( mu.data[i] );
            for ( std::size_t i = 0; i < n; ++i )
              mu.data[i] *= kToDeg;
            auto& h = vectAt(data.hists,data.histidx_theta);
            hgfill_main( h, mu, b.w, n );
            if ( dethistid )
              hgfill_detail( h, *dethistid, mu, b.w, n );
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

const NCMMCT::tally_hist_t&
NCMMCT::TallyStdHists_Data::accessHistogram( StrView histname,
                                             Optional<DetailedHistsID>
                                             detailid ) const
{
  //fixme: unit test
  const tally_hist_t * histptr = nullptr;
  for ( auto& h : this->hists ) {
    if ( histname != h.main.title() )
      continue;
    if ( !detailid.has_value() ) {
      histptr = &h.main;
    } else {
      auto histidx = static_cast<std::size_t>( detailid.value() );
      if ( h.detailed.empty() || !(histidx < h.detailed.size()) )
        NCRYSTAL_THROW( BadInput,
                        "Detailed breakdown histograms not available" );
      histptr = & vectAt(h.detailed,histidx);
    }
    break;
  }
  if ( !histptr )
    NCRYSTAL_THROW2( BadInput,
                     "Tally histogram not available: \""<<histname<<"\".");
  return *histptr;
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
                          const BasketValBufBool& sawinelas )
{
  using TFlags = TallyFlags::Flags;
  const std::size_t n = neutrons.size();

  //////////////////////////////////
  //Find DetailedHistID if required:
  DetailedHistsIDVect histid_buf;
  const DetailedHistsIDVect * histid = nullptr;
  if ( !opt.flags.has(TFlags::nobreakdown) ) {
    init_dethistidvect( histid_buf, nscat, sawinelas, n );
    histid = &histid_buf;
  }

  //////////////////////////////
  //Fill cos(mu), theta, q tallies:
  constexpr auto flag_miscmu = (TFlags::mu|TFlags::theta|TFlags::q );
  if ( opt.flags.hasAny(flag_miscmu) )
    tallyRecord_mu_theta_q( data, opt, neutrons, histid );

  constexpr auto flag_nonmiscmu = ( TFlags::ALLHISTS & (~flag_miscmu) );
  if ( !opt.flags.hasAny(flag_nonmiscmu) )
    return;

  ///////////////
  //Fill e tally:
  if ( opt.flags.has(TFlags::e) ) {
    auto& h = vectAt(data.hists,data.histidx_e);
    hgfill_main( h, neutrons.ekin, neutrons.w, n );
    if ( histid )
      hgfill_detail( h, *histid, neutrons.ekin, neutrons.w, n );
  }

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
    hgfill_main( h, nscat, neutrons.w, n );
    if ( histid )
      hgfill_detail( h, *histid, nscat, neutrons.w, n );
  }

  if ( opt.flags.has(TFlags::nscat_uw) ) {
    auto& h = vectAt(data.hists,data.histidx_nscat_uw);
    hgfill_main( h, nscat, n );
    if ( histid )
      hgfill_detail( h, *histid, nscat, n );
  }

  if ( opt.flags.has(TFlags::w) ) {
    auto& h = vectAt(data.hists,data.histidx_w);
    hgfill_main( h, neutrons.w, n );
    if ( histid )
      hgfill_detail( h, *histid, neutrons.w, n );
  }

}
