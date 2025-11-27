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

void NCMMCU::calcProbTransm( NumberDensity nd, std::size_t N,
                             const double * ncrestrict xs_or_nullptr,
                             const double * ncrestrict dist,
                             double * ncrestrict out )
{
  if ( !xs_or_nullptr ) {
    for ( auto i : ncrange( N ) )
      out[i] = 1.0;
    return;
  }
  for ( auto i : ncrange( N ) )
    out[i] = macroXS( nd, CrossSect{xs_or_nullptr[i]} );
  for ( auto i : ncrange( N ) )
    out[i] *= dist[i];
  for ( auto i : ncrange( N ) )
    out[i] = -out[i];
  for ( auto i : ncrange( N ) )
    out[i] = std::exp( out[i] );
}

void NCMMCU::propagate( NeutronBasket& b,
                        const double* ncrestrict dists )
{
  for ( auto i : ncrange( b.size() ) )
    b.x[i] += dists[i] * b.ux[i];
  for ( auto i : ncrange( b.size() ) )
    b.y[i] += dists[i] * b.uy[i];
  for ( auto i : ncrange( b.size() ) )
    b.z[i] += dists[i] * b.uz[i];
}

void NCMMCU::propagateAndAttenuate( NeutronBasket& b,
                                    NumberDensity nd,
                                    const double* ncrestrict dists,
                                    const double* ncrestrict xsvals )
{
  propagate( b, dists );
  if ( xsvals ) {
    for ( auto i : ncrange( b.size() ) )
      b.w[i] *= std::exp( -macroXS( nd, CrossSect{ xsvals[i] } ) * dists[i] );
  }
}

void NCMMCU::sampleRandDists( RNG& rng, NumberDensity nd,
                              const double * ncrestrict dists,
                              const double * ncrestrict xsvals,
                              std::size_t N,
                              double * ncrestrict tgt )
{
  NewABI::generateMany( rng, N, tgt );
  //TODO: Better vectorisation possible?
  for ( auto i : ncrange(N) ) {
    RandExpIntervalSampler rs( 0.0, dists[i],
                               macroXS( nd, CrossSect{ xsvals[i] } ) );
    tgt[i] = rs.sample( tgt[i] );
  }
}

void NCMMCU::scatterGivenMu( RNG& rng,
                             NeutronBasket& b,
                             double * ncrestrict mu_vals )
{
  for ( auto i : ncrange(b.size()) ) {
    //TODO: better vectorisation possible?
    auto newdir
      = randNeutronDirectionGivenScatterMu( rng,
                                            CosineScatAngle{ mu_vals[i] },
                                            NeutronDirection{ b.ux[i],
                                                              b.uy[i],
                                                              b.uz[i] } );
    b.ux[i] = newdir[0];
    b.uy[i] = newdir[1];
    b.uz[i] = newdir[2];
  }
}

NCMMCU::ScenarioDecoded NCMMCU::decodeScenario( const MatCfg& matcfg,
                                                const char* scenario )
{
  //FIXME: More specific and helpful exception messages.

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

  if ( parts.empty() )
    return decodeScenario( matcfg, "0.8BT" );

  auto it = parts.begin();
  auto itE = parts.end();
  StrView sv_energy = *(it++);
  bool is_pencil = false;
  if ( it!=itE && *it == "pencil" ) {
    is_pencil = true;
    ++it;
  }

  auto bad = [sv_input]()
  {
    //fixme: remove this general thing
    NCRYSTAL_THROW2(BadInput,
                    "Invalid MiniMC scenario string: \""<<sv_input<<'"');
  };

  bool is_sphere(true);
  StrView sv_thickness = "1mfp";
  if ( it != itE ) {
    if ( *it != "on" )
      bad();
    ++it;
    if ( it == itE )
      bad();
    if ( *std::prev(itE) == "sphere" ) {
      //selected sphere (already the default)
      --itE;
    } else if ( *std::prev(itE) == "slab" ) {
      //selected slab (already the default)
      is_sphere = false;
      --itE;
    }
    if ( it != itE )
      sv_thickness = *(it++);
    if ( it != itE )
      bad();
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
    bad();
  StrView e_unit = parse_e.value().unit;
  NeutronEnergy neutron_energy;
  NeutronWavelength neutron_wavelength;
  bool wavelength_mode = false;
  OptionalInfoPtr info;
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
    } else if ( e_unit == "BT" ) {
      //unit is bragg threshold (if any)
        if (!info)
          info = FactImpl::createInfo( matcfg );
        Optional<NeutronWavelength> braggthr;
        if ( info->isSinglePhase() )
          braggthr = info->getBraggThreshold();
        if ( !braggthr.has_value() ) {
          //fixme: in case of multiphase, look though all the phases and take
          //the longest bragg threshold found?
          braggthr = NeutronWavelength{ 4.0 };
        }
        wavelength_mode = true;
        double bt_Aa = braggthr.value_or( NeutronWavelength{ 4.0 } ).get();
        neutron_wavelength = NeutronWavelength( round6(e_val * bt_Aa) );
        neutron_energy = neutron_wavelength.energy();
    } else {
      NCRYSTAL_THROW2(BadInput,
                      (e_unit.empty()?"Missing":"Invalid")
                      <<" energy unit on \"" <<sv_energy<<"\"."
                      << " Possible units include Aa, meV, eV, or"
                      " BT (=Bragg threshold in Aa, or 4Aa if not available).");
    }
    neutron_wavelength = neutron_energy.wavelength();
  }

  //Finally decode thickness:
  double thickness_meter = -1.0;
  auto parse_t = Cfg::unitSplit(sv_thickness);
  if (!parse_t.has_value() || parse_t.value().unit.empty() )
    bad();
  if ( parse_t.value().unit == "mfp" ) {
    //fixme: round to nearest 3 digits? Or perhaps just return descriptive
    //string we can use in titles etc.?
    //We need Info+Scatter for this.
    if ( info == nullptr )
      info = FactImpl::createInfo( matcfg );
    CachePtr cache;
    auto scatter = FactImpl::createScatter( matcfg );
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
  //ok, we now have thickness, energy, is_sphere/is_slab and is_pencil. Time to
  //translate to actual cfg strings:

  const char * g15 = "%.15g";

  ScenarioDecoded res;
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
      ss_src << "constant;z=" << fmt(src_z_plus_epsilon,g15);
    else
      ss_src << "circular;z=" << fmt(src_z);
    ss_src << ";r="<<fmt(sphere_r);
  } else {
    nc_assert_always( is_slab );
    double box_dz = thickness_meter;
    double box_dxy = std::max(1.0,thickness_meter*1e5);
    ss_geom << "box;dz="<<fmt(box_dz)<<";dx="
            <<fmt(box_dxy)<<";dy="<<fmt(box_dxy);
    nc_assert_always( is_pencil );
    ss_src << "constant;z=" << fmt(src_z_plus_epsilon,g15);
  }

  //Add energy to src:
  if ( wavelength_mode ) {
    ss_src << ";wl="<<fmt(neutron_wavelength.get());
  } else {
    ss_src << ";ekin="<<fmt(neutron_energy.get());
  }

  res.srccfg = ss_src.str();
  res.geomcfg = ss_geom.str();
  res.enginecfg = "";
  return res;
}
