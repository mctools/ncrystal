#ifndef NCrystal_ExtnDummy_hh
#define NCrystal_ExtnDummy_hh

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

#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/internal/extn_utils/NCExtnEval.hh"

namespace NCRYSTAL_NAMESPACE {

  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // Dummy extinction model for testing.                                      //
  //                                                                          //
  // It supports no secondary extinction, and extinction factors are given as //
  // E=max(0,1-<blocksize>/100um). In other words, going linearly from 1.0 to //
  // 0.0 as blocksize is increased from 0um to 100um.                         //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////


  namespace Extn {

    struct DummyMdlPurePrimary final {

      struct ModelData {
        double fixed_E;
      };

      static ModelData initModelData( const PowderBraggInput::CellData& cell,
                                      Length blockSize )
      {
        ModelData res;
        double x = blockSize.get() / ( 100.0 * Length::micrometer );
        res.fixed_E = std::min( 1.0, std::max<double>( 0.0, 1.0 - x ) );
        return res;
      }

      struct NeutronData {
        double wlsq;
      };

      static NeutronData initNeutronData( const ModelData& md,
                                          NeutronEnergy ekin )
      {
        NeutronData res;
        res.wlsq = ekin2wlsq( ekin.get() );
        return res;
      }

      static double extinctionFactorScalar( const ModelData& md,
                                            const NeutronData& nd,
                                            const double fsq,
                                            const double inv2dsp )
      {
        return md.

        //DblBasket sabine_x, El, Eb, sinth_sq, costh_sq;
        const double sabine_x = nd.factor_fsq2x * fsq;
        const double El = calcSabineEl_y0( sabine_x );
        const double Eb = calcSabineEb_y0( sabine_x );
        double sinth_sq = nd.wlsq * ncsquare( inv2dsp );//fixme: cache inv2dspsq ?
        //Fixme: just do max(0.0,...) on final E value just before returning?
        //Or have the framework clip the values to [0,1] for absolute certainty!
        sinth_sq = std::min<double>(1.0,sinth_sq);
        double costh_sq = 1.0 - sinth_sq;
        costh_sq = std::max<double>(costh_sq,0.0);
        return El*costh_sq + Eb*sinth_sq;
      }

      template<class TDblBasket, class TFDMBasket>
      static TDblBasket extinctionFactor( const ModelData& md,
                                          const NeutronData& nd,
                                          const FDMBasket& fdm )
      {
        constexpr auto N = DblBasket::count;
        DblBasket sabine_x, El, Eb, sinth_sq, costh_sq;
        for ( auto i : ncrange( N ) )
          sabine_x[i] = nd.factor_fsq2x * fdm.fsq[i];

        for ( auto i : ncrange( N ) )
          El[i] = calcSabineEl_y0( sabine_x[i] );

        for ( auto i : ncrange( N ) )
          Eb[i] = calcSabineEb_y0( sabine_x[i] );

        for ( auto i : ncrange( N ) )
          sinth_sq[i] = nd.wlsq * ncsquare( fdm.inv2dsp[i] );//fixme: cache inv2dspsq ?

        //Fixme: just do max(0.0,...) on final E value just before returning?
        //Or have the framework clip the values to [0,1] for absolute certainty!
        for ( auto i : ncrange( N ) )
          sinth_sq[i] = std::min<double>(1.0,sinth_sq[i]);

        for ( auto i : ncrange( N ) )
          costh_sq[i] = 1.0 - sinth_sq[i];

        for ( auto i : ncrange( N ) )
          costh_sq[i] = std::max<double>(costh_sq[i],0.0);

        //E = El*cos^2(th) + Eb*sin^2(th):
        for ( auto i : ncrange( N ) )
          El[i] *= costh_sq[i];
        for ( auto i : ncrange( N ) )
          Eb[i] *= sinth_sq[i];
        for ( auto i : ncrange( N ) )
          El[i] += Eb[i];
        return El;
      }


      template<TAccumulator>
      void accumulatePlaneContribs( NeutronEnergy, TAccumulator& acc )
      {
        if ( ekin < m_threshold )
          return;
        const auto wl = ekin.wavelength();
        const double v0 = m_data.cell.volume;
        const double kkk = ncsquare(wl.get()*m_domainSizeAa/v0)*1e-8;



        const double wlhalf = wl.get()*0.5;
        const double factor
          = ncsquare(wl.get())/(2.0*m_data.cell.volume*m_data.cell.n_atoms);


        StableSum contrib;
        for ( auto& e : m_data.planes ) {
          if ( e.dsp < wlhalf )
            break;



          double sabine_x = kkk * e.fsq;
          double El = calcSabineEl_y0( sabine_x );
          double Eb = calcSabineEb_y0( sabine_x );
          double sinth_sq = std::min<double>(1.0,ncsquare(0.5 * wl.get() / e.dsp));
          double costh_sq = std::max<double>(1.0 - sinth_sq,0.0);
          double extinction_correction = El * costh_sq + Eb * sinth_sq;
          contrib.add( e.dsp * e.fsq * e.mult * extinction_correction  );
        }
        return CrossSect{ factor * contrib.sum()  };


      }

std::vector<double>&

Length blockSize );
    };

    double calcSabineEb( double x, double y );//Sabine 6.4.5.5
    double calcSabineEl( double x, double y );//Improved version of Sabine
                                              //eqs. 6.4.5.3-4 as discussed
                                              //above. See "Original" functions
                                              //below for reference.

    //Same but y=0 or Eb with pre-calculated A(y) and B(y) factors:
    double calcSabineEb_y0( double x );
    double calcSabineEl_y0( double x );
    double calcSabineA( double y );//Sabine 6.4.5.6
    double calcSabineB( double y );//Sabine 6.4.5.7
    double calcSabineEb_cachedAB( double x, double A, double B );

    //Original recipe for El as it appeared in Sabine's paper:
    double calcSabineElOriginal( double x, double y );
    double calcSabineElOriginal_y0( double x );

    //We also provide evaluations of the secondary extinction factors in the
    //uncorrelated model, for both rectangular and triangular tilt functions:
    double calcSabineEl_ScndRect( double x, double y );//Sabine 6.4.9.2
    double calcSabineEb_ScndRect( double x, double y );//Sabine 6.4.9.3
    double calcSabineEl_ScndTriang( double x, double y );//Sabine 6.4.9.4
    double calcSabineEb_ScndTriang( double x, double y );//Sabine 6.4.9.5

    //Same but y=0 or Eb with pre-calculated A(y) and B(y) factors:
    double calcSabineEl_ScndRect_y0( double x );
    double calcSabineEb_ScndRect_y0( double x );
    double calcSabineEl_ScndTriang_y0( double x );
    double calcSabineEb_ScndTriang_y0( double x );
    double calcSabineEb_ScndRect_cachedAB( double x, double A, double B );
    double calcSabineEb_ScndTriang_cachedAB( double x, double A, double B );
  }

}

////////////////////////////
// Inline implementations //
////////////////////////////

inline double NCrystal::Extn::calcSabineEl( double x, double y )
{
  nc_assert( y>=0.0 );
  nc_assert( std::isfinite(y) );
  return std::exp(-y)*calcSabineEl_y0(x);
}

inline double NCrystal::Extn::calcSabineElOriginal( double x, double y )
{
  nc_assert( y>=0.0 );
  nc_assert( std::isfinite(y) );
  return std::exp(-y)*calcSabineElOriginal_y0(x);
}

inline double NCrystal::Extn::calcSabineEb_y0( double x )
{
  nc_assert( x>=0.0 );
  nc_assert( std::isfinite(x) );
  return 1.0 / std::sqrt( 1.0 + x );
}

inline double NCrystal::Extn::calcSabineEb_cachedAB( double x,
                                                     double A,
                                                     double B )
{
  nc_assert( x>=0.0 );
  nc_assert( std::isfinite(x) );
  nc_assert( A>=0.0 );
  nc_assert( std::isfinite(A) );
  nc_assert( B>=0.0 );
  nc_assert( std::isfinite(B) );
  return A / std::sqrt( 1.0 + B * x );
}

inline double NCrystal::Extn::calcSabineEb( double x, double y )
{
  nc_assert( y>=0.0 );
  nc_assert( std::isfinite(y) );
  nc_assert( x>=0.0 );
  nc_assert( std::isfinite(x) );
  return calcSabineA( y ) / std::sqrt( 1.0 + calcSabineB( y ) * x );
}

inline double NCrystal::Extn::calcSabineEb_ScndRect_y0( double x )
{
  nc_assert( x>=0.0 );
  nc_assert( std::isfinite(x) );
  return 1.0 / ( 1.0 + x );
}

inline double NCrystal::Extn::calcSabineEl_ScndRect( double x, double y )
{
  nc_assert( y>=0.0 );
  nc_assert( std::isfinite(y) );
  return std::exp(-y)*calcSabineEl_ScndRect_y0(x);
}

inline double NCrystal::Extn::calcSabineEb_ScndRect_cachedAB( double x,
                                                              double A,
                                                              double B )
{
  nc_assert( x>=0.0 );
  nc_assert( std::isfinite(x) );
  nc_assert( A>=0.0 );
  nc_assert( std::isfinite(A) );
  nc_assert( B>=0.0 );
  nc_assert( std::isfinite(B) );
  return A / ( 1.0 + B*x );
}

inline double NCrystal::Extn::calcSabineEb_ScndRect( double x, double y )
{
  nc_assert( x>=0.0 );
  nc_assert( std::isfinite(x) );
  const double A = calcSabineA( y );
  const double B = calcSabineB( y );
  return A / ( 1.0 + B*x );
}

inline double NCrystal::Extn::calcSabineEl_ScndTriang( double x, double y )
{
  nc_assert( y>=0.0 );
  nc_assert( std::isfinite(y) );
  return std::exp(-y)*calcSabineEl_ScndTriang_y0(x);
}

inline double NCrystal::Extn::calcSabineEb_ScndTriang_cachedAB( double x,
                                                                double A,
                                                                double B )
{
  nc_assert( x>=0.0 );
  nc_assert( std::isfinite(x) );
  nc_assert( A>=0.0 );
  nc_assert( std::isfinite(A) );
  nc_assert( B>=0.0 );
  nc_assert( std::isfinite(B) );
  return A * calcSabineEb_ScndTriang_y0( B * x );
}

inline double NCrystal::Extn::calcSabineEb_ScndTriang( double x, double y )
{
  nc_assert( x>=0.0 );
  nc_assert( std::isfinite(x) );
  const double A = calcSabineA( y );
  const double B = calcSabineB( y );
  return A * calcSabineEb_ScndTriang_y0( B * x );
}

#endif
