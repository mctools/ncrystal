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

#include "NCrystal/internal/minimc/NCMMC_Utils.hh"
#include "NCrystal/internal/extd_utils/NCABIUtils.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCrystal/internal/cfgutils/NCCfgTypes.hh"
#include "NCrystal/factories/NCFactImpl.hh"

namespace NC = NCrystal;
namespace NCMMCU = NCrystal::MiniMC::Utils;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    inline void safe_mult_strongzero( double * __restrict a,
                                      const double * __restrict b,
                                      std::size_t n )
    {
      //Assuming all a[i], b[i]>=0 and that a[i] is infinite, implements a[i] *=
      // b[i] in a way so the result is always 0 if a[i] is zero. Specifically,
      // this is intended to allow the case of (a[i],b[i])=(0,inf) to yield 0
      // without triggering Nan or FPE.
      //
      //It was tested that this seems to vectorize nicely.
      for ( std::size_t i = 0; i < n; ++i ) {
        nc_assert( std::isfinite(a[i]) );
        nc_assert( b[i] >= 0.0 );
        a[i] *= ( a[i] ? b[i] : 0.0 );
      }
    }
  }
}

void NCMMCU::calcProbTransm( NumberDensity nd, std::size_t n,
                             bool geom_is_unbounded,
                             const double * ncrestrict xs_or_nullptr,
                             const double * ncrestrict dists,
                             double * ncrestrict out )
{
  //Calculate probability of transmission. Note that if the geometry is
  //unbounded, we need to consider that dist=inf is a possibility. In that case
  //we want the transmission probability to be 0 (to not tally something not
  //leaving the geometry).

  //Fixme: Ensure we unit test the behaviour here, especially regarding
  //geom_is_unbounded, dists[i]=inf, and xs=0.0. Optionally xs=inf would be nice
  //to get to work (should give ptransm=0)
  if ( !xs_or_nullptr ) {
    if ( geom_is_unbounded ) {
      for ( auto i : ncrange( n ) ) {
        out[i] = static_cast<double>(!std::isinf(dists[i]));
      }
#ifndef NDEBUG
      for ( auto i : ncrange( n ) ) {
        nc_assert_always(out[i]==(std::isinf(dists[i])?0.0:1.0));
      }
#endif
    } else {
      for ( auto i : ncrange( n ) )
        out[i] = 1.0;
    }
    return;
  }

  for ( auto i : ncrange( n ) )
    out[i] = macroXS( nd, CrossSect{xs_or_nullptr[i]} );
  if ( geom_is_unbounded ) {
    safe_mult_strongzero( out, dists, n );//fixme would fail if xs=inf
  } else {
    for ( auto i : ncrange( n ) )
      out[i] *= dists[i];
  }
  for ( auto i : ncrange( n ) )
    out[i] = -out[i];
  for ( auto i : ncrange( n ) )
    out[i] = std::exp( out[i] );
  if ( geom_is_unbounded ) {
    for ( auto i : ncrange( n ) )
      out[i] *= (!std::isinf(dists[i]));
  }
#ifndef NDEBUG
  for ( auto i : ncrange( n ) ) {
    nc_assert_always( !std::isnan(out[i]) );
    nc_assert_always( ( !std::isinf(dists[i]) ) || out[i]==0.0 );
  }
#endif
}

void NCMMCU::propagate( NeutronBasket& b,
                        bool geom_is_unbounded,
                        const double* ncrestrict dists )
{
  if ( geom_is_unbounded ) {
    //dists might be inf, so we can't just multiply direction vectors with
    //dists, or we might end up with 0*inf.
    double tmp[basket_N];
    std::copy( &b.ux[0], &b.ux[0]+b.size(), &tmp[0] );
    safe_mult_strongzero( tmp, dists, b.size() );
    for ( auto i : ncrange( b.size() ) )
      b.x[i] += tmp[i];
    std::copy( &b.uy[0], &b.uy[0]+b.size(), &tmp[0] );
    safe_mult_strongzero( tmp, dists, b.size() );
    for ( auto i : ncrange( b.size() ) )
      b.y[i] += tmp[i];
    std::copy( &b.uz[0], &b.uz[0]+b.size(), &tmp[0] );
    safe_mult_strongzero( tmp, dists, b.size() );
    for ( auto i : ncrange( b.size() ) )
      b.z[i] += tmp[i];
  } else {
    for ( auto i : ncrange( b.size() ) )
      b.x[i] += dists[i] * b.ux[i];
    for ( auto i : ncrange( b.size() ) )
      b.y[i] += dists[i] * b.uy[i];
    for ( auto i : ncrange( b.size() ) )
      b.z[i] += dists[i] * b.uz[i];
  }
}

void NCMMCU::propagateAndAttenuate( NeutronBasket& b,
                                    NumberDensity nd,
                                    bool geom_is_unbounded,
                                    const double* ncrestrict dists,
                                    const double* ncrestrict xsvals )
{
  //FIXME: if !xsvals but we have infinities, we need to check for those!!!!

  propagate( b, geom_is_unbounded, dists );
  if ( xsvals ) {
    double tmp[basket_N];
    for ( auto i : ncrange( b.size() ) )
      tmp[i] = macroXS( nd, CrossSect{ xsvals[i] } );
    if ( geom_is_unbounded ) {
      safe_mult_strongzero( tmp, dists, b.size() );
    } else {
      for ( auto i : ncrange( b.size() ) )
        tmp[i] *= dists[i];
    }
    for ( auto i : ncrange( b.size() ) )
      tmp[i] = -tmp[i];
    for ( auto i : ncrange( b.size() ) )
      tmp[i] = std::exp( tmp[i] );
    for ( auto i : ncrange( b.size() ) ) {
      b.w[i] *= tmp[i];
    }
    //fixme: for the special case of macroxs=0 and dist=inf we need to set w=0.
  }

  if ( geom_is_unbounded ) {
    //Make sure that particles propagating to infinity without scattering are
    //always "lost" to the tallies by setting w=0. This is needed for when cross
    //sections are 0 but distances are infinite.
    for ( auto i : ncrange( b.size() ) )
      b.w[i] *= (1.0-static_cast<double>(std::isinf(dists[i])));
  }
}

void NCMMCU::sampleRandDists( RNG& rng, NumberDensity nd,
                              const double * ncrestrict dists,
                              const double * ncrestrict xsvals,
                              std::size_t N,
                              double * ncrestrict tgt )
{
  NewABI::generateMany( rng, N, tgt );
  for ( auto i : ncrange(N) ) {
    double macroxs = macroXS( nd, CrossSect{ xsvals[i] } );
    if ( !macroxs ) ncunlikely {
      tgt[i] = kInfinity;
      continue;
    }
    if ( std::isinf(dists[i]) ) ncunlikely {
      tgt[i] = std::log( tgt[i] )/(-macroxs);
    } else {
      RandExpIntervalSampler rs( 0.0, dists[i], macroxs );
      tgt[i] = rs.sample( tgt[i] );
    }
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {
    namespace {
      std::string detail_format_count( StrView sv_count )
      {
        nc_assert_always( sv_count.has_value() );
        auto count_dbl = sv_count.toDbl();
        std::size_t count_i = 0;
        static_assert( sizeof(std::size_t) >= 8, "" );

        //value must be able to fit in std::uint64_t, i.e. be less than
        //~1.8446..e19. We pick 1e19 as maximum.
        if ( count_dbl.has_value()
             && count_dbl.value()>0.0
             && count_dbl.value()<=1e19 ) {
          count_i = static_cast<std::size_t>(count_dbl.value()+0.5);
          if ( count_i != count_dbl.value() )
            count_i = 0;
        }
        if (!count_i)
          NCRYSTAL_THROW2(BadInput,"Invalid count specification \""
                          <<sv_count<<"\". Count must be a positive"
                          " integral value (and at most 1e19).");
        return fmtUInt64AsNiceDbl( static_cast<std::uint64_t>(count_i) );
      }

      Optional<NeutronWavelength>
      detail_extract_braggThreshold( const Info& info ) {
        Optional<NeutronWavelength> res;
        if ( info.isSinglePhase() ) {
          res = info.getBraggThreshold();
        } else {
          //multiphase: take the longest!
          for ( auto& phase : info.getPhases() ) {
            auto bt = phase.second->getBraggThreshold();
            if ( !bt.has_value() )
              continue;
            if ( !res.has_value() || res.value() < bt.value() )
              res = bt;
          }
        }
        return res;
      }
    }
  }
}

NCMMCU::ScenarioDecoded NCMMCU::decodeScenario( const MatCfg& matcfg,
                                                const char* scenario )
{
  ScenarioDecoded res;

  StrView sv_input((scenario?scenario:""));
  auto badchar_strrep = findForbiddenChar( sv_input,
                                           Cfg::forbidden_chars_value_strreps,
                                           ExtraForbidOpt::RequireSimpleASCII );
  if ( badchar_strrep.has_value() )
    NCRYSTAL_THROW2(BadInput,"Forbidden character "<<badchar_strrep.value()
                    <<" in MiniMC scenario string: \""<<sv_input<<'"');

  //Split on whitespace + ":_":
  auto parts = sv_input.split_any<
    8,StrView::SplitKeepEmpty::No,StrView::SplitTrimParts::Yes>(" \t\n\r:_");

  OptionalInfoPtr info;
  ProcImpl::OptionalProcPtr scatter;

  if ( parts.empty() ) {
    info = FactImpl::createInfo( matcfg );
    const bool has_BT = detail_extract_braggThreshold(*info).has_value();
    const bool has_T = info->hasTemperature();
    const char * default_scenario = ( has_BT
                                      ? "0.8BT"
                                      : ( has_T ? "1kT" : "1.8Aa" ) );
    res = decodeScenario( matcfg, default_scenario );
    return res;
  }

  auto it = parts.begin();
  auto itE = parts.end();

  //Parse ENERGY and pencil:
  StrView sv_energy = *(it++);
  bool is_pencil = false;
  if ( it!=itE && *it == "pencil" ) {
    is_pencil = true;
    ++it;
  }

  //Parse COUNT:
  Optional<StrView> sv_count;
  if ( it!=itE && *std::prev(itE) == "times" ) {
    --itE;
    if ( it == itE )
      NCRYSTAL_THROW2(BadInput,"Invalid MiniMC scenario string: \"" << sv_input
                      <<"\". Expected integral value like \"10000\" or"
                      " \"1e6\" in front of keyword \"times\".");
    sv_count = *std::prev(itE);
    --itE;
  }
  std::string count_formatted;
  if ( sv_count.has_value()) {
    count_formatted = detail_format_count( sv_count.value() );
  } else {
    if ( scatter == nullptr )
      scatter = FactImpl::createScatter( matcfg );
    count_formatted = ( scatter->isOriented() ? "1e5" : "1e6" );
  }

  bool is_sphere(true);
  StrView sv_thickness = "1mfp";
  if ( it != itE ) {
    if ( *it != "on" )
      NCRYSTAL_THROW2(BadInput,"Invalid MiniMC scenario string: \""<<sv_input
                      <<"\". Expected keyword \"on\" but got \""<<*it<<"\".");
    ++it;
    bool had_params_after_on = false;

    if ( it != itE ) {
      //more parts, scrabe off final "sphere" or "slab" if found.
      if ( *std::prev(itE) == "sphere" ) {
        //selected sphere (already the default)
        --itE;
        had_params_after_on = true;
      } else if ( *std::prev(itE) == "slab" ) {
        //selected slab (already the default)
        had_params_after_on = true;
        is_sphere = false;
        --itE;
      }
    }
    //Now collect a thickness:
    if ( it == itE ) {
      if (!had_params_after_on)
        NCRYSTAL_THROW2(BadInput,"Invalid MiniMC scenario string: \""<<sv_input
                        <<"\". Missing parameters after keyword \"on\".");
    } else {
      sv_thickness = *(it++);
      if ( it != itE )
        NCRYSTAL_THROW2(BadInput,"Invalid MiniMC scenario string: \""<<sv_input
                        <<"\". Unexpected parameter: \""<<*it<<"\"");
    }
  }

  const bool is_slab = !is_sphere;
  if ( is_slab && !is_pencil )
    is_pencil=true;//always pencil beam in case of infinite slab.

  auto round6 = []( double val ) -> double
  {
    std::ostringstream ss;
    ss<<fmt(val,"%.6g");
    std::string sval = ss.str();
    return str2dbl(sval);
  };

  //Finally decode energy:
  auto parse_e = Cfg::unitSplit(sv_energy);
  if (!parse_e.has_value())
    NCRYSTAL_THROW2(BadInput,
                    "Invalid energy specification in \"" <<sv_energy<<"\".");
  StrView e_unit = parse_e.value().unit;
  NeutronEnergy neutron_energy;
  NeutronWavelength neutron_wavelength;
  bool wavelength_mode = false;
  const double e_val = parse_e.value().value;
  if ( e_unit == "Aa" ) {
    wavelength_mode = true;
    neutron_wavelength = NeutronWavelength{ e_val };
    neutron_energy = neutron_wavelength.energy();
  } else {
    if ( e_unit == "eV" ) {
      neutron_energy = NeutronEnergy{ e_val };
    } else if ( e_unit == "meV" ) {
      neutron_energy = NeutronEnergy{ e_val * 0.001 };
    } else if ( e_unit == "keV" ) {
      neutron_energy = NeutronEnergy{ e_val * 1e3 };
    } else if ( e_unit == "MeV" ) {
      neutron_energy = NeutronEnergy{ e_val * 1e6 };
    } else if ( e_unit == "GeV" ) {
      neutron_energy = NeutronEnergy{ e_val * 1e9 };
    } else if ( e_unit == "kT" ) {
      if (!info)
        info = FactImpl::createInfo( matcfg );
      if (!info->hasTemperature())
        NCRYSTAL_THROW(BadInput,"Can not use kT as a unit for a multi-phase "
                       "material without a single well-defined temperature.");
      const double kT = info->getTemperature().kT();
      nc_assert_always( kT > 0.0 );
      wavelength_mode = false;
      neutron_energy = NeutronEnergy{ round6(e_val * kT) };
      neutron_wavelength = neutron_energy.wavelength();
    } else if ( e_unit == "BT" ) {
      //unit is bragg threshold (or 4Aa)
      if (!info)
        info = FactImpl::createInfo( matcfg );
      wavelength_mode = true;
      auto braggthr = detail_extract_braggThreshold( *info );
      double bt_Aa = braggthr.value_or( NeutronWavelength{ 4.0 } ).get();
      neutron_wavelength = NeutronWavelength( round6(e_val * bt_Aa) );
      neutron_energy = neutron_wavelength.energy();
    } else {
      NCRYSTAL_THROW2(BadInput,
                      (e_unit.empty()?"Missing":"Invalid")
                      <<" energy unit on \"" <<sv_energy<<"\"."
                      << " Possible units include Aa, meV, eV,"
                      " kT (=kB*temperature),"
                      " BT (=Bragg threshold in Aa, or 4Aa if not available).");
    }
    neutron_wavelength = neutron_energy.wavelength();
  }

  //Finally decode thickness:
  double thickness_meter = -1.0;
  auto parse_t = Cfg::unitSplit(sv_thickness);
  if (!parse_t.has_value()) {
    NCRYSTAL_THROW2(BadInput,"Invalid thickness specification in \""
                    <<sv_thickness<<"\".");
  }

  if (parse_t.value().unit.empty() ) {
    std::ostringstream ss;
    Cfg::units_length::listAvailableUnitsNoDefault(ss);
    NCRYSTAL_THROW2(BadInput,"Missing length unit on: \""<<sv_thickness<<"\". "
                    <<"Accepts either the special unit"
                    " \"mfp\" (mean-free-path between scatterings) or one of"
                    " the standard length units: "<<ss.str());
  }
  if ( parse_t.value().unit == "mfp" ) {
    if ( info == nullptr )
      info = FactImpl::createInfo( matcfg );
    CachePtr cache;
    if ( scatter == nullptr )
      scatter = FactImpl::createScatter( matcfg );
    auto xs = scatter->crossSection( cache, neutron_energy,
                                     NeutronDirection{ 0.0, 0.0, 1.0 } );
    auto numdens = info->getNumberDensity();
    double macroxs = macroXS( numdens, xs );// [1/m]
    if ( macroxs < 1e-99 || macroxs > 1e99 ) {
      thickness_meter = 0.01;//fall back value of 1cm
    } else {
      thickness_meter = round6( parse_t.value().value / macroxs );
      //Limit precision slightly, to avoid FP instabilities:

    }
  } else {
    auto pv = Cfg::units_length::parse(sv_thickness);
    if ( !pv.has_value() ) {
      std::ostringstream ss;
      Cfg::units_length::listAvailableUnitsNoDefault(ss);
      NCRYSTAL_THROW2(BadInput,"invalid length: \""<<sv_thickness<<"\". "
                      <<"Must be a value followed by either the special unit"
                      " \"mfp\" (mean-free-path between scatterings) or one of"
                      " the standard length units: "<<ss.str());
    }
    thickness_meter = pv.value().first * 1e-10;//units_length::parse returns Aa
  }
  if ( thickness_meter < 1e-120 || thickness_meter > 1e120 )
    NCRYSTAL_THROW2(BadInput,"Could not determine suitable thickness from \""
                    <<sv_thickness<<"\" (thickness is invalid or out of range)");

  ///////////////////////////////////////////////////////////////////
  //Ok, we now have count, thickness, energy, is_sphere/is_slab and
  //is_pencil. Time to translate to actual cfg strings:

  const char * g15 = "%.15g";

  std::ostringstream ss_src;
  std::ostringstream ss_geom;
  const double src_z = - thickness_meter*0.5;
  //Slightly more efficient if we start pencil beams inside the geometry
  //(1e-14 is small but will not be obscured by "%.15g".
  const double src_z_plus_epsilon = src_z * (1.0-1e-14);
  if ( is_sphere ) {
    nc_assert_always( !is_slab );
    const double sphere_r = thickness_meter*0.5;
    ss_geom << "sphere;r="<<fmt(sphere_r);
    if ( is_pencil )
      ss_src << "constant;z=" << fmt(src_z_plus_epsilon,g15);//fixme: do away
                                                             //with the epsilon?
                                                             //Also needs unit
                                                             //test that
                                                             //proptovolentry
                                                             //always gives a
                                                             //proptovolexit of
                                                             //the correct
                                                             //magnitude in this
                                                             //case.
    else
      ss_src << "circular;z=" << fmt(src_z);
    ss_src << ";r="<<fmt(sphere_r);
  } else {
    nc_assert_always( is_slab );
    const double slab_dz = thickness_meter*0.5;
    ss_geom << "slab;dz="<<fmt(slab_dz);
    nc_assert_always( is_pencil );
    ss_src << "constant;z=" << fmt(src_z);
  }

  //Add energy to src:
  if ( wavelength_mode ) {
    ss_src << ";wl="<<fmt(neutron_wavelength.get());
  } else {
    ss_src << ";ekin="<<fmt(neutron_energy.get());
  }

  //Add count:
  ss_src << ";n=" << count_formatted;

  res.srccfg = ss_src.str();
  res.geomcfg = ss_geom.str();
  res.enginecfg = "";

  return res;
}

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {
    namespace {
#ifndef NDEBUG
      void check_x_ux_thickness_sane( const double * ncrestrict x,
                                      const double * ncrestrict ux,
                                      std::size_t n,
                                      double thickness )
      {
        nc_assert( n <= basket_N );
        nc_assert( thickness > 0.0);
        nc_assert( std::isfinite(thickness) );
        for( std::size_t i = 0; i < n; ++i ) {
          nc_assert( std::isfinite(x[i]) );
          nc_assert( ux[i] >= -1.0 );
          nc_assert( ux[i] <= 1.0 );
        }
      }
#endif
    }
  }
}

void NCMMCU::distToSlabExit( const double * ncrestrict x,
                             const double * ncrestrict ux,
                             double * ncrestrict out_dist,
                             std::size_t n,
                             double slab_halfthickness )
{
  const double dx = slab_halfthickness;
#ifndef NDEBUG
  check_x_ux_thickness_sane( x, ux, n, dx );
#endif
  for( std::size_t i = 0; i < n; ++i ) {
    if ( ux[i] > 0.0 ) {
      out_dist[i] = ( dx-x[i] )/ux[i];
    } else {
      if ( ux[i] < 0.0 ) nclikely {
        out_dist[i] = -( dx+x[i] )/ux[i];
      } else ncunlikely {
        out_dist[i] = kInfinity;
      }
    }
  }
}

void NCMMCU::distToSlabEntry( const double * ncrestrict x,
                              const double * ncrestrict ux,
                              double * ncrestrict out_dist,
                              std::size_t n,
                              double slab_halfthickness )
{
  const double dx = slab_halfthickness;
#ifndef NDEBUG
  check_x_ux_thickness_sane( x, ux, n, dx );
#endif

  // a = |x|-dx
  double a[basket_N];
  for( std::size_t i = 0; i < n; ++i )
    a[i] = ncabs( x[i] );
  for( std::size_t i = 0; i < n; ++i )
    a[i] -= dx;

  double x_ux[basket_N];
  for( std::size_t i = 0; i < n; ++i )
    x_ux[i] = x[i] * ux[i];

  for( std::size_t i = 0; i < n; ++i ) {
    if ( a[i] <= 0 ) {
      // |x|<=dx inside or at edge.
      if ( a[i] ) nclikely {
        out_dist[i] = 0.0;//inside
      } else {
        //edge (unlikely)
        out_dist[i] = ( x_ux[i] > 0.0
                        ? -1.0 //points out
                        : 0.0//points in or skirts along edge
                        );
      }
    } else {
      //|x|>dx, outside.
      if ( x_ux[i] >= 0.0 ) {
        //outside, misses.
        out_dist[i] = -1.0;
      } else {
        //outside, hits. We know that ux is non-zero and has opposite sign of x.
        out_dist[i] = a[i] / ncabs(ux[i]);
      }
    }
  }
}
