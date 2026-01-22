#ifndef NCrystal_ExtnBC_hh
#define NCrystal_ExtnBC_hh

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

#include "NCrystal/internal/extn_utils/NCExtnBC1974.hh"
#include "NCrystal/internal/extn_utils/NCExtnBC2025.hh"

namespace NCRYSTAL_NAMESPACE {

  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // Utilities for selecting Becker-Coppens recipes based on enums selecting  //
  // recipe version and components.                                           //
  //                                                                          //
  // BCEval (any: P=Primary, G=ScndGauss, L=ScndLorentz, F=ScndFresnel).      //
  // BCScndEval (just the secondary models, G/L/F).                           //
  // BCPrimaryEval (just primary, P).                                         //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////

  namespace Extn {

    struct BCXCalc_P final {
      struct ModelData final {
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
      class NeutronXCalc final {
        double m_factor_fsq2x; // [1/barn]
      public:
        NeutronXCalc( const ModelData& md, NeutronWavelength wl )
          : m_factor_fsq2x( ncsquare(wl.get() ) * md.k1 )
        {
        }
        template<class TPlaneData>
        double calcX( const TPlaneData& p, const double /*sintheta*/ ) const
        {
          return m_factor_fsq2x * p.fsq;
        }
      };
    };

    struct BCXCalc_G final {
      struct ModelData final {
        double k1; // [ 1/(barn * Aa^2) ]
        double k2; // [ Aa^2 ]
      };
      static ModelData initModelData( const PowderBraggInput::CellData& cell,
                                      Length grainSize,
                                      MosaicityFWHM angularSpread )
      {
        //2*pi*g^2 = 1/(2*sigma^2) <=>g^2 = 1/(4*pi*sigma^2)
        const double sigma = angularSpread.sigma().get();
        nc_assert_always(sigma>0.0);
        const double inv_g_squared = k4Pi * ncsquare(sigma);
        NCRYSTAL_MSG("TKTEST g="<<std::sqrt(1.0/inv_g_squared));

        //X = (2/3)*Q*alphaG*Tbar
        //alphaG = alphabar/sqrt( 1+alphabar^2/(2*g^2))
        //alphabar = (3/4)*(4/3)*tbar*(sin2theta/lambda)
        //Q = |F/V|^2*lambda^3/sin2theta
        // => x = (2/3)*|F/V|^2*lambda^2*tbar/sqrt( 1+alphabar^2/(2*g^2)))*Tbar
        //      = (2/3)*|F/V|^2*lambda^2*tbar*TBar * 1/sqrt( 1+alphabar^2/(2*g^2))
        //
        //So same as in primary case except for final factor. To get k1 in units
        //of 1/(barn*Aa^2) we must multiply with barn/Aa^2 = 1e-8:

        const double tbar_Aa = grainSize.get() / Length::angstrom;
        ModelData res;
        res.k1 = ( 2.0 / 3.0 ) * ncsquare( tbar_Aa / cell.volume ) * 1e-8;

        //Now:
        //
        //1/sqrt( 1+alphabar^2/(2*g^2))
        //   = 1/sqrt( 1+(tbar*(sin2theta/lambda))^2/(2*g^2)) =
        //   = 1/sqrt( 1+0.5*((tbar*sin2theta/(lambda*g))^2))
        //sin2theta = 2*sinth*sqrt(1-sinth^2), lambda/2d = sinth,

        //Now we want to define k2 as everything in the last term except
        //sin2theta/lambda:

        res.k2 = 0.5 * ncsquare(tbar_Aa)*inv_g_squared;//[Aa^2]

        return res;
      }
      class NeutronXCalc final {
        double m_k1wlsq = 0.0; // [1/barn]
        double m_k2 = 0.0; // [1/barn]
      public:
        NeutronXCalc( const ModelData& md, NeutronWavelength ekin )
          : m_k1wlsq(ekin2wlsq( ekin.get() ) * md.k1),
            m_k2(md.k2)
        {
        }
        template<class TPlaneData>
        double calcX( const TPlaneData& p, const double sintheta ) const
        {
          //const double inv2dsq = ncsquare(p.inv2d);
          //ksl = sin(2*theta) / lambda:
          //sin2theta = 2*sinth*sqrt(1-sinth^2), lambda/2d = sinth,
          //(sin2theta/lambda)^2 = 4*sinth^2*(1-sinth^2)/ ( (2d)^2 * sinth^2)
          //                     = 4*(1-sinth^2)/(2d)^2
          const double ksl = 4.0*( 1.0 - ncsquare(sintheta) )*ncsquare(p.inv2dsp);
          return m_k1wlsq * p.fsq / std::sqrt( std::max(1e-199,1.0 + m_k2 * ksl));

          return ( ( 2.0 / 3.0 ) * ncsquare( tbar_Aa / cell.volume ) * 1e-8 * wl^2 * p.fsq )

            / std::sqrt( 1.0 + 2 * ncsquare(tbar_Aa)*( 1.0 - ncsquare(sintheta) ) / ((2d)^2*g^2)


          //sin2th^2 = 4*sinth^2*(1-sinth^2)
                         //sinth = wl/D2




          //shuqi:
         xs = 2. / 3. * Q_theta * L / [(wl/l)*std::sqrt(1 + NC::ncsquare(l / wl)*NC::ncsquare(sin_2theta) / (2. * g * g))];




        }
      };
    };


    enum class BC_RecipeVersion { Classic1974 = 0, Std2025 = 1, Lux2025 = 2,
                                  Default = Std2025 };
    enum class BC_M{ P, G, L, F };//Primary, ScndGauss, ScndLorentz, ScndFresnel
    using BC_Component = BC_M;

    template<BC_M, BC_RecipeVersion = BC_RecipeVersion::Default>
    struct BCEval {};

    template<> struct BCEval<BC_M::P,BC_RecipeVersion::Classic1974> { static double eval( double x, double sinth ) { return BC1974::y_primary(x,sinth); } };
    template<> struct BCEval<BC_M::P,BC_RecipeVersion::Std2025> { static double eval( double x, double sinth ) { return BC2025::y_primary(x,sinth); } };
    template<> struct BCEval<BC_M::P,BC_RecipeVersion::Lux2025> { static double eval( double x, double sinth ) { return BC2025::y_primary_lux(x,sinth); } };
    template<> struct BCEval<BC_M::G,BC_RecipeVersion::Classic1974> { static double eval( double x, double sinth ) { return BC1974::y_scndgauss(x,sinth); } };
    template<> struct BCEval<BC_M::G,BC_RecipeVersion::Std2025> { static double eval( double x, double sinth ) { return BC2025::y_scndgauss(x,sinth); } };
    template<> struct BCEval<BC_M::G,BC_RecipeVersion::Lux2025> { static double eval( double x, double sinth ) { return BC2025::y_scndgauss_lux(x,sinth); } };
    template<> struct BCEval<BC_M::L,BC_RecipeVersion::Classic1974> { static double eval( double x, double sinth ) { return BC1974::y_scndlorentz(x,sinth); } };
    template<> struct BCEval<BC_M::L,BC_RecipeVersion::Std2025> { static double eval( double x, double sinth ) { return BC2025::y_scndlorentz(x,sinth); } };
    template<> struct BCEval<BC_M::L,BC_RecipeVersion::Lux2025> { static double eval( double x, double sinth ) { return BC2025::y_scndlorentz_lux(x,sinth); } };
    template<> struct BCEval<BC_M::F,BC_RecipeVersion::Classic1974> { static double eval( double x, double sinth ) { return BC1974::y_scndfresnel(x,sinth); } };
    template<> struct BCEval<BC_M::F,BC_RecipeVersion::Std2025> { static double eval( double x, double sinth ) { return BC2025::y_scndfresnel(x,sinth); } };
    template<> struct BCEval<BC_M::F,BC_RecipeVersion::Lux2025> { static double eval( double x, double sinth ) { return BC2025::y_scndfresnel_lux(x,sinth); } };

    //Just primary:

    template<BC_RecipeVersion recV = BC_RecipeVersion::Default>
    struct BCPrimaryEval {
      static double eval( double x, double sinth )
      {
        return BCEval<BC_M::P,recV>::eval(x,sinth);
      }
    };

    //Just secondary:

    enum class BC_scndM{ G, L, F };//ScndGauss, ScndLorentz, ScndFresnel
    constexpr BC_scndM BC_M_to_scndM( BC_M m )
    {
      //      static_assert( m == BC_M::G || m == BC_M::L || m == BC_M::F, "" );
      return nc_assert_rv(m == BC_M::G || m == BC_M::L || m == BC_M::F),
        ( m == BC_M::G ? BC_scndM::G
          : ( m == BC_M::L ? BC_scndM::L : BC_scndM::F ) );
    }

    constexpr BC_M BC_scndM_to_M( BC_scndM m )
    {
      return ( m == BC_scndM::G ? BC_M::G
               : ( m == BC_scndM::L ? BC_M::L : BC_M::F ) );
    }

    template< BC_scndM scndM,
              BC_RecipeVersion recV = BC_RecipeVersion::Default >
    struct BCScndEval {
      static double eval( double x, double sinth )
      {
        constexpr BC_M m = BC_scndM_to_M(scndM);
        return BCEval<m,recV>::eval(x,sinth);
      }
    };

  }
}

#endif
