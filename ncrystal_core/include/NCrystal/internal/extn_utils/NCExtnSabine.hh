#ifndef NCrystal_ExtnSabineMdl_hh
#define NCrystal_ExtnSabineMdl_hh

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

//#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/internal/utils/NCExtraTypes.hh"
#include "NCrystal/internal/utils/NCMsg.hh"//fixme
#include "NCrystal/internal/extn_utils/NCExtnUtils.hh"
//#include "NCrystal/internal/extn_utils/NCExtnEval.hh"

namespace NCRYSTAL_NAMESPACE {

  //////////////////////////////////////////////////////////////////////////////
  //
  // Implementation of models described in:
  //
  // Sabine, T. M. "The flow of radiation in a real crystal." (2006): 609-616.
  // doi: 10.1107/97809553602060000603
  //
  //////////////////////////////////////////////////////////////////////////////

  namespace Extn {

    struct SabineMdlPurePrimary final {

      struct ModelData {
        double k1;//unit is [ 1/(barn * Aa^2) ]
      };

      static ModelData initModelData( const PowderBraggInput::CellData& cell,
                                      Length blockSize )
      {
        ModelData res;
        double blksz_Aa = blockSize.get() / Length::angstrom;

        //To get k1 in units of 1/(barn*Aa^2) we must multiply with
        //barn/Aa^2 = 1e-8:
        res.k1 = ncsquare( blksz_Aa / cell.volume ) * 1e-8;
        return res;
      }

      struct NeutronData {
        double wlsq;  // [Aa^2]
        double factor_fsq2x; // [fixme]
      };

      static NeutronData initNeutronData( const ModelData& md,
                                          NeutronEnergy ekin )
      {
        NeutronData res;
        res.wlsq = ekin2wlsq( ekin.get() );
        res.factor_fsq2x = res.wlsq * md.k1;
        return res;
      }

      struct PlaneData {
        //Fields required by all models:
        double dsp;
        double fdm;
        //Other fields, as needed by this particular model:
        double fsq;//structure factor squared [barn]
        double inv2dsp_squared; // = (2*dspacing)^(-2)
      };

      static PlaneData initPlaneData( const PowderBraggInput::Plane& p )
      {
        PlaneData res;
        res.dsp = p.dsp;
        res.fdm = p.fsq * p.dsp * p.mult;
        res.fsq = p.fsq;
        res.inv2dsp_squared = ncsquare( 1.0 / ( 2.0 * p.dsp ) );
        return res;
      }

      static double extinctionFactor( const ModelData& /*used through neutron data*/,
                                      const NeutronData& n,
                                      const PlaneData& p )
      {
        const double sabine_x = n.factor_fsq2x * p.fsq;
        const double El = calcSabineEl_y0( sabine_x );
        const double Eb = calcSabineEb_y0( sabine_x );
        const double sinth_sq = n.wlsq * p.inv2dsp_squared;
        return El*(1.0-sinth_sq) + Eb*sinth_sq;
      }
    };

    struct SabineMdlCorrelatedScnd final {

      static double calc_g( MosaicityFWHM m )
      {
        //fixme: to util hdr
        constexpr double k = 1.0 / (2* kSqrtPi);// ~= 0.8862269254527579
        //Note: to get 2*sigma (2*stddev), we would instead have had to multiply
        //with ~0.85 (4.3% less), so a bit annoying that there is yet another
        //value, almost the same.
        return k / m.sigma().get();//fixme: not fwhm? Seems like eq. 6.4.8.2
                                   //indicates sigma.
      }

      struct ModelData {
        double k1;//unit is [ 1/(sqrt(barn)*Aa) ]
        double k2;//unit is [ 1/(barn*Aa^3) ]
      };

      // struct ScndCfg {
      //   Length grainSize
      //   MosaicityFWHM angularSpread;
      // };

      static ModelData initModelData( const PowderBraggInput::CellData& cell,
                                      Length blockSize,
                                      Length grainSize,
                                      MosaicityFWHM angularSpread )
      {
        nc_assert_always( grainSize > blockSize );
        nc_assert_always( blockSize.get() > 0.0 );//fixme: allow 0 in cfg strings!
        ModelData res;
        double g = calc_g( angularSpread );
        //g = 500.0;//fixme
        NCRYSTAL_MSG("TKTEST g="<<g);
        const double blksz_Aa = blockSize.get() / Length::angstrom;
        const double grnsz_Aa = grainSize.get() / Length::angstrom;
        NCRYSTAL_MSG("TKTEST blksz_Aa="<<blksz_Aa);
        NCRYSTAL_MSG("TKTEST grnsz_Aa="<<grnsz_Aa);
        //To get k1 in units of 1/sqrt(barn*Aa^2) we must multiply with
        //sqrt(barn)/angstrom = 1e-4:
        res.k1 = 1e-4 * blksz_Aa / cell.volume;
        //To get k2 in units of 1/(barn*Aa^3) we must multiply with
        //barn/angstrom^2 = 1e-8:
        NCRYSTAL_MSG("TKTEST cell.volume = "<<cell.volume);
        res.k2 = g * 1e-8 * ( grnsz_Aa - blksz_Aa ) / ncsquare(cell.volume);//fixme square????
        return res;
      }

      struct NeutronData {
        double wlsq;
        double k1_mult_wl;
        double k2_mult_wlsq;
      };

      static NeutronData initNeutronData( const ModelData& md,
                                          NeutronEnergy ekin )
      {
        const double wlsq = ekin2wlsq( ekin.get() );
        const double wl = std::sqrt(wlsq);
        NeutronData res;
        res.wlsq = wlsq;
        res.k1_mult_wl = md.k1 * wl;
        res.k2_mult_wlsq = md.k2 * wlsq;
        //Normal:
        // x = (Nc*F*wl*l)^2 = fsq*wl^2*(Nc*l)^2
        // Here: Nc*l is (blksz_Aa / cell.volume)*1e-4
        //
        //Now we get instead (ASSUMING SABINE MEANT N_c after eq. 6.4.8.6):
        //
        // x = (Nc*F*wl*l + g*Q0*(L-l))^2

        // Q0 = Nc^2*wl^3*F^2/sin20 = (Nc^2*wl^2*F^2) * wl/sin2th
        //    == (Nc*wl*F)^2 * wl/(2*sinth*costh)
        //    == (Nc*wl*F)^2 * dhkl/costh
        //
        // Thus:
        //
        // x = ( Nc*l*wl*F + g*Nc^2*(L-l)*F^2*dhkl*wl^2/costh )^2
        //
        //  ModelData: k1 = Nc*l
        //             k2 = g*Nc^2*(L-l)
        //
        //  NeutronData:
        //               wl, wl^2
        //
        //  With plane (cached: f, fsq, dhk) we must find: costh and then x.
        //
        //  We should make sure to reuse costh that we anyway need for combining later.

        //shuqi:
        //double sin_theta = 0.5 * wl / d_hkl; //2*d_hkl*sin(theta_hkl)=wl
        //and wl/sin_theta = wl/(2*sin_theta*cos_theta)
        //                 = dhkl / costheta


        // wl^3/sin(2theta) = 2 * wl^2 * d_hkl


        //lambda = 2*dhkl * sintheta (theta in 0..pi/2)

        //T=theta (max pi/2, this is not scattering angle, D=2*dhkl
        //Bragg condition: wl = D*sin(T)
        //sin(2T) = 2*sinT*sqrt(1-sinT^2)

        //Bragg condition: wl = D*sinT, sinT=wl/D

        //So wl^3/sin(2T) = wl^3 / ( 2*(wl/D)*sqrt(1-(wl/D)^2) )
        //
        //If u=wl/D = sinTh:

        //wl/sin2T = D / sqrt(1-u^2) = D/cosT

        //
        //What if wl ~= D, i.e. cosT ~= 0? ( backscat)
        //
        //
        //wl^3/sin(2T) = (0.5*D^3) * u^2 / sqrt(1-u^2)
        //
        //With y = u^2 = (sinT)^2 we get:
        //
        //wl^3/sin(2T) = (0.5*D^3) * y / sqrt(1-y)
        //
        //Which tends to infinity as y->1 (backscattering). So backscattering has infinite sabine_x (meaning El->0,Eb->

        return res;
      }

      struct PlaneData {
        //Fields required by all models:
        double dsp;
        double fdm;
        //Other fields, as needed by this particular model:
        double fsq;//structure factor squared [barn]
        double f;//sqrt(fsq)
        double inv2dsp_squared; // = (2*dspacing)^(-2)
      };

      static PlaneData initPlaneData( const PowderBraggInput::Plane& p )
      {
        PlaneData res;
        res.dsp = p.dsp;
        res.fdm = p.fsq * p.dsp * p.mult;
        res.fsq = p.fsq;
        res.f = std::sqrt( p.fsq );
        res.inv2dsp_squared = ncsquare( 1.0 / ( 2.0 * p.dsp ) );
        return res;
      }

      static double extinctionFactor( const ModelData& /*used through neutron data*/,
                                      const NeutronData& n,
                                      const PlaneData& p )
      {
        const double sinth_sq = n.wlsq * p.inv2dsp_squared;
        if ( sinth_sq > (1.0-1e-14) ) {
          //Backscattering limit (sabine x->0, E->0):
          //(fixme: better to make the check on what we divide with, i.e. costh )?
          return 0.0;
        };
        const double costh_sq = std::max<double>(0.0,1.0 - sinth_sq);
        const double costh = std::sqrt(costh_sq);
        const double sabine_x = ncsquare( n.k1_mult_wl * p.f
                                          + n.k2_mult_wlsq * p.fsq * p.dsp / costh);//fixme absorp 2 on k2, and cache fsq*dsp?

        //wl/sin2th = wl/2*sinth*costh = wl/(2*(wl/2d)*costh) = d/costh

        //const double sabine_x = n.factor_fsq2x * p.fsq;
        const double El = calcSabineEl_y0( sabine_x );
        const double Eb = calcSabineEb_y0( sabine_x );
        return El*costh_sq + Eb*sinth_sq;
      }
    };


    template<bool IS_RECTANGULAR>
    struct SabineMdlUncorrelatedScnd final {

      static double calc_G_Rectangular( MosaicityFWHM m )
      {
        //fixme: to utils

        // RMS of Sabine's rectangular function is 1/(2*sqrt(3)*G), so we
        // determine the value of G which would match up that RMS with that of a
        // Gaussian with the given mosaicity.

        const double rms = m.sigma().get();
        nc_assert_always( rms > 0.0 );
        constexpr double k = 1.0 / ( constexpr_sqrt(12.0) ); // ~= 0.2887..
        return k / rms;
      }

      static double calc_G_Triangular( MosaicityFWHM m )
      {
        //fixme: to utils

        // RMS of Sabine's triangular function is 1/(sqrt(6)*G), so we determine
        // the value of G which would match up that RMS with that of a Gaussian
        // with the given mosaicity.

        const double rms = m.sigma().get();
        nc_assert_always( rms > 0.0 );
        constexpr double k = 1.0 / constexpr_sqrt(6.0); // ~= 0.4082..
        return k / rms;
      }

      struct ModelData {
        double k1;//unit is [ 1/(barn*Aa^2) ]
        double k2;//unit is [ 1/(barn*Aa^3) ]
      };

      // struct ScndCfg {
      //   Length grainSize
      //   MosaicityFWHM angularSpread;
      // };

      static ModelData initModelData( const PowderBraggInput::CellData& cell,
                                      Length blockSize,
                                      Length grainSize,
                                      MosaicityFWHM angularSpread )
      {

        nc_assert_always( grainSize > blockSize );
        nc_assert_always( blockSize.get() > 0.0 );//fixme: allow 0 in cfg strings!
        const double blksz_Aa = blockSize.get() / Length::angstrom;
        const double grnsz_Aa = grainSize.get() / Length::angstrom;

        //To get k1 in units of 1/(barn*Aa^2) we must multiply with
        //barn/Aa^2 = 1e-8:
        ModelData res;
        res.k1 = ncsquare( blksz_Aa / cell.volume ) * 1e-8;

        //Secondary extinction:
        const double G = ( IS_RECTANGULAR
                           ? calc_G_Rectangular( angularSpread )
                           : calc_G_Triangular( angularSpread ) );
        nc_assert_always( G > 0.0 );
        //G = 500.0;//fixme
        NCRYSTAL_MSG("TKTEST G="<<G);

        //Factor of 1e-8 is to get k1 in units of 1/(barn*Aa^3):
        res.k2 = 1e-8 * G * ncsquare( 1.0 / cell.volume ) * grnsz_Aa; //FIXME: consider (grnsz_Aa-blksz_Aa) ??
        return res;
      }

      struct NeutronData {
        double wlsq;  // [Aa^2]
        double factor_fsq2x; // [fixme]
        double k2_mult_wlsq; // [fixme]
      };

      static NeutronData initNeutronData( const ModelData& md,
                                          NeutronEnergy ekin )
      {
        NeutronData res;
        res.wlsq = ekin2wlsq( ekin.get() );
        res.factor_fsq2x = res.wlsq * md.k1;
        res.k2_mult_wlsq = res.wlsq * md.k2;
        return res;
      }

      struct PlaneData {
        //Fields required by all models:
        double dsp;
        double fdm;
        //Other fields, as needed by this particular model:
        double fsq;//structure factor squared [barn]
        double inv2dsp_squared; // = (2*dspacing)^(-2)
      };

      static PlaneData initPlaneData( const PowderBraggInput::Plane& p )
      {
        PlaneData res;
        res.dsp = p.dsp;
        res.fdm = p.fsq * p.dsp * p.mult;
        res.fsq = p.fsq;
        res.inv2dsp_squared = ncsquare( 1.0 / ( 2.0 * p.dsp ) );
        return res;
      }

      static double extinctionFactor( const ModelData& /*used through neutron data*/,
                                      const NeutronData& n,
                                      const PlaneData& p )
      {
        //First primary extinction:
        const double sabine_x = n.factor_fsq2x * p.fsq;
        const double El = calcSabineEl_y0( sabine_x );
        const double Eb = calcSabineEb_y0( sabine_x );
        const double sinth_sq = n.wlsq * p.inv2dsp_squared;
        if ( sinth_sq > (1.0-1e-14) ) {
          //Backscattering limit (sabine x->0, E->0): (fixme: better to make the
          //check on what we divide with later, i.e. costh )?
          return 0.0;
        };
        const double costh_sq = std::max<double>(0.0,1.0 - sinth_sq);
        const double costh = std::sqrt(costh_sq);
        const double Eprimary = El*costh_sq + Eb*sinth_sq;

        //Then secondary extinction:
        const double sabine_x_scnd
          = Eprimary * n.k2_mult_wlsq * p.fsq * p.dsp / costh;//fixme cache fsq*dsp?

        const double El_scnd = ( IS_RECTANGULAR
                                 ? calcSabineEl_ScndRect_y0( sabine_x_scnd )
                                 : calcSabineEl_ScndTriang_y0( sabine_x_scnd ) );
        const double Eb_scnd = ( IS_RECTANGULAR
                                 ? calcSabineEb_ScndRect_y0( sabine_x_scnd )
                                 : calcSabineEb_ScndTriang_y0( sabine_x_scnd ) );
        const double Escnd = El_scnd*costh_sq + Eb_scnd*sinth_sq;
        return Eprimary * Escnd;
      }
    };

    using SabineMdlUncorrelatedScnd_Rec = SabineMdlUncorrelatedScnd<true>;
    using SabineMdlUncorrelatedScnd_Tri = SabineMdlUncorrelatedScnd<false>;

    //fixme: adding non-Sabine models here for now, needs cleanup or file
    //renaming later.

    enum class BC_YpParameterisation { Classic1974, ClassicUpdated2025, Lux2025 };
    template<BC_YpParameterisation TModel = BC_YpParameterisation::Lux2025>
    struct BCMdlPurePrimary final {

      struct ModelData {
        double k1;//unit is [ 1/(barn * Aa^2) ]
      };

      static ModelData initModelData( const PowderBraggInput::CellData& cell,
                                      Length blockSize )
      {
        ModelData res;
        double blksz_Aa = blockSize.get() / Length::angstrom;
        //To get k1 in units of 1/(barn*Aa^2) we must multiply with
        //barn/Aa^2 = 1e-8:
        res.k1 = ( 2.0 / 3.0 ) * ncsquare( blksz_Aa / cell.volume ) * 1e-8;
        return res;
      }

      struct NeutronData {
        double wlsq;  // [Aa^2]
        double factor_fsq2x; // [fixme]
      };

      static NeutronData initNeutronData( const ModelData& md,
                                          NeutronEnergy ekin )
      {
        NeutronData res;
        res.wlsq = ekin2wlsq( ekin.get() );
        res.factor_fsq2x = res.wlsq * md.k1;
        return res;
      }

      struct PlaneData {
        //Fields required by all models:
        double dsp;
        double fdm;
        //Other fields, as needed by this particular model:
        double fsq;//structure factor squared [barn]
        double inv2dsp_squared; // = (2*dspacing)^(-2)
      };

      static PlaneData initPlaneData( const PowderBraggInput::Plane& p )
      {
        PlaneData res;
        res.dsp = p.dsp;
        res.fdm = p.fsq * p.dsp * p.mult;
        res.fsq = p.fsq;
        res.inv2dsp_squared = ncsquare( 1.0 / ( 2.0 * p.dsp ) );
        return res;
      }

      static double extinctionFactor( const ModelData& /*used through neutron data*/,
                                      const NeutronData& n,
                                      const PlaneData& p )
      {
        const double bc_x = n.factor_fsq2x * p.fsq;
        const double sinth_sq = n.wlsq * p.inv2dsp_squared;

        constexpr bool is_classic = (TModel == BC_YpParameterisation::Classic1974);
        constexpr bool is_updatedclassic = (TModel == BC_YpParameterisation::ClassicUpdated2025);
        constexpr bool is_lux2025 = (TModel == BC_YpParameterisation::Lux2025);
        constexpr bool is_any_classic = ( is_classic || is_updatedclassic );

        if ( is_any_classic ) {
          const double cos2th = 1.0 - 2.0 * sinth_sq;
          double A, B;
          if ( is_classic ) {
            A = 0.20 + 0.45 * cos2th;
            B = 0.22 - 0.12 * ncsquare(0.5 - cos2th);
          } else {
            nc_assert(is_updatedclassic);
            A = 0.559 + 0.537 * cos2th;
            B = 0.604 - 0.222 * ncsquare(0.5 - cos2th);
          }
          double t = 1. + B * bc_x;
          if ( std::fabs( t ) < 1e-20 )
            t = ( t > 0 ? 1e-20 : -1e-20 );
          const double u = 1. + 2. * bc_x + ncsquare( bc_x ) * A / t;
          return 1. / std::sqrt(std::max<double>(1.0,u));//fixme: clamp u to
                                                         //[1,..] to avoid
                                                         //extinction factors
                                                         //larger than 1 (and to
                                                         //avoid FPEs in this
                                                         //line).
        } else {
          //fixme: to utils header?:
          nc_assert(is_lux2025);
          const double inv1px = 1.0 / ( 1.0 + bc_x );
          const double s = std::sqrt(sinth_sq);
          const double A = 0.8968+s*(-0.2928+s*(2.528+s*(-11.21+s*(21.07+s*(-20.17+7.32*s)))));
          const double B = 0.4697+s*(-0.1688+s*(5.37+s*(-20.68+s*(40.52+s*(-38.09+12.89*s)))));
          const double C = -0.7486+s*(-0.5671+s*(-1.248+s*(17.56+s*(-47.8+s*(59.14-26.61*s)))));
          const double t = 1. + bc_x * ( B + C * inv1px );
          nc_assert_always( t > 0.0 );//fixme _always
          const double u = 1. + 2. * bc_x + ( ncsquare( bc_x ) * A - 0.1 * bc_x ) / t;
          nc_assert_always(u>=0.0);//fixme _always
          return 1. / std::sqrt(u);
        }
      }
    };

  }
}

#endif
