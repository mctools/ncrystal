#ifndef NCrystal_ExtnSabine_hh
#define NCrystal_ExtnSabine_hh

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

//#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/internal/utils/NCExtraTypes.hh"
#include "NCrystal/internal/utils/NCMsg.hh"//fixme
#include "NCrystal/internal/extn_utils/NCExtnUtils.hh"
#include "NCrystal/internal/extn_utils/NCExtnBC2025.hh"//fixme: untangle from sabine hdr
//#include "NCrystal/internal/extn_utils/NCExtnEval.hh"
#include "NCrystal/internal/extn_utils/NCExtnMdlBC.hh"//fixme: untangle from sabine hdr

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

    //In Sabine's 1988 comparison paper, he finds in section 10.b that to yield
    //consistent results, his definition of x must be exactly twice of that of
    //the Becker-Coppens x-value. Since the Becker-Coppens x-definition contains
    //a geometric factor of 2/3 compared to the standard Sabine x definition, we
    //modify the Sabine x value by a factor of 4/3. Since x \propto blocksize^2,
    //this means block-sizes fitted with the Sabine model will change by a
    //factor of sqrt(3)/2!=93%.
    //Fixme: consider this some more. Especially how it should be used in
    //secondary extinction models.
    constexpr static double sabine_x_extra_factor = 4./3.;

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
        res.k1 = sabine_x_extra_factor* ncsquare( blksz_Aa / cell.volume ) * 1e-8;
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
        double sabine_x = n.factor_fsq2x * p.fsq;
        sabine_x /= sabine_x_extra_factor;
        sabine_x *= 1.0/std::sqrt(1e-14+n.wlsq * p.inv2dsp_squared);//fixme
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
        //fixme: not sure if sabine_x_extra_factor should be used here????
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
        //fixme: sabine_x_extra_factor for secondary models???
        const double sabine_x = sabine_x_extra_factor * ncsquare( n.k1_mult_wl * p.f
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
        res.k1 = sabine_x_extra_factor * ncsquare( blksz_Aa / cell.volume ) * 1e-8;

        //Secondary extinction:
        const double G = ( IS_RECTANGULAR
                           ? calc_G_Rectangular( angularSpread )
                           : calc_G_Triangular( angularSpread ) );
        nc_assert_always( G > 0.0 );
        //G = 500.0;//fixme
        NCRYSTAL_MSG("TKTEST G="<<G);

        //NB: Not using sabine_x_extra_factor here?!?
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

  }
}

#endif
