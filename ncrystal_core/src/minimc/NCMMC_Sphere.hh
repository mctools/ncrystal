#ifndef NCrystal_MMC_Sphere_hh
#define NCrystal_MMC_Sphere_hh

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

#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/minimc/NCMMC_Basket.hh"
#include "NCrystal/internal/minimc/NCMMC_Defs.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace {

    inline double fast_sqrt_clippos( double x ) noexcept {
      //Using this in a loop does not actually allow for vectorization,
      //https://gcc.gnu.org/bugzilla/show_bug.cgi?id=91645 . However, perhaps
      //this will be solved in GCC 13/14 or is already solved in clang/msvc?
      //Further investigations needed.
#if defined(__clang__) || defined(__GNUC__)
      return __builtin_sqrt(__builtin_fmax(0.0,x));
#else
      return std::sqrt(std::max<double>(0.0,x));
#endif
    }
  }

  namespace MiniMC {

    constexpr double geom_tolerance_factor = 1e-12;

    class Sphere {
      //double m_radius;
      double m_radiusSq;
    public:
      Sphere( Length radius )
        : m_radiusSq( ncsquare( radius.dbl() ) )
      {
        nc_assert_always( radius.dbl() > 0.0 );
        nc_assert_always( m_radiusSq < 1e199 );
        nc_assert_always( m_radiusSq > 0.0 );
      }

      bool pointIsInside( const Vector& v ) const
      {
        return v.mag2() <= m_radiusSq;
      }

      void distToVolumeEntry( const NeutronBasket& nb,
                              Span<double> tgt ) const
      {
        nc_assert( tgt.size() >= nb.nused);
        distToVolumeEntryImpl( nb.x, nb.y, nb.z,
                               nb.ux, nb.uy, nb.uz,
                               tgt.data(), nb.nused );
      }

      void distToVolumeExit( const NeutronBasket& nb,
                             Span<double> tgt ) const
      {
        nc_assert( tgt.size() >= nb.nused);
        distToVolumeExitImpl( nb.x, nb.y, nb.z,
                              nb.ux, nb.uy, nb.uz,
                              tgt.data(), nb.nused );
      }

      static void unit_test()
      {
        {
          const double x[]  = { -30.0, -30.0, 30.0, 0.0, 0.0, };
          const double y[]  = {   0.0,   0.0,  0.0, 0.0, 0.0, };
          const double z[]  = {   0.0,   0.0,  0.0, 0.0, 10.0*(1.0-1.0e-14), };
          const double ux[] = {   1.0,   0.0,  1.0, 0.0, 0.0, };
          const double uy[] = {   0.0,  -1.0,  0.0, 1.0, 0.0, };
          const double uz[] = {   0.0,   0.0,  0.0, 0.0, 1.0 };
          const double dist_to_entry[] = { 20.0, -1.0, -1.0, 0.0, 0.0 };
          constexpr std::size_t n = sizeof(x) / sizeof(*x);
          double buf[n];
          Sphere(Length{10.0}).distToVolumeEntryImpl( x,y,z,ux,uy,uz,buf,n);
          for ( std::size_t i = 0; i < n; ++i ) {
            nc_assert_always(floateq(buf[i],dist_to_entry[i]));
          }
        }
        {
          const double x[]  = { -9.999,  0.0,  5.0, 9.999, 0.0, -10.0, -10.0, -10.0, -10.0 };
          const double y[]  = {   0.0,   0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
          const double z[]  = {   0.0,   0.0,  0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0 };
          const double ux[] = {   1.0,   0.0,  1.0, 1.0, 0.0, 0.0, -1.0, 1.0, kInvSqrt2 };
          const double uy[] = {   0.0,  -1.0,  0.0, 0.0, 0.0, 1.0, 0.0, 0.0, kInvSqrt2 };
          const double uz[] = {   0.0,   0.0,  0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
          const double dist_to_exit[] = { 19.999, 10.0, 5.0, 0.001, 0.0, 0.0, 0.0, 20.0, 10.0*kSqrt2 };
          constexpr std::size_t n = sizeof(x)/sizeof(*x);
          double buf[n];
          Sphere(Length{10.0}).distToVolumeExitImpl( x,y,z,ux,uy,uz,buf,n);
          for ( std::size_t i = 0; i < n; ++i ) {
            nc_assert_always(floateq(buf[i],dist_to_exit[i]));
          }
        }
      }

    private:
      void distToVolumeEntryImpl( const double * ncrestrict x,
                                  const double * ncrestrict y,
                                  const double * ncrestrict z,
                                  const double * ncrestrict ux,
                                  const double * ncrestrict uy,
                                  const double * ncrestrict uz,
                                  double * ncrestrict tgt,
                                  std::size_t n ) const ncnoexceptndebug
      {
        //TODO: Properly test this for auto-vectorisation (at least check with
        //-fopt-info-vec-missed).
#ifndef NDEBUG
        for ( std::size_t i = 0; i < n; ++i ) {
          nc_assert( std::isfinite(x[i]) && std::isfinite(y[i]) && std::isfinite(z[i]) );
          nc_assert( std::isfinite(ux[i]) && std::isfinite(uy[i]) && std::isfinite(uz[i]) );
          nc_assert( ncabs(ncsquare(ux[i])+ncsquare(uy[i])+ncsquare(uz[i]) - 1.0) < 1e-13 );
        }
#endif
        for ( std::size_t i = 0; i < n; ++i ) {
          const double psq_mr2 = ncsquare(x[i])+ncsquare(y[i])+ncsquare(z[i]) - m_radiusSq;
          if ( psq_mr2 <= 0.0 ) {
            tgt[i] = 0.0;
          } else {
            const double pd = (x[i]) * (ux[i]) + (y[i]) * (uy[i]) +(z[i]) * (uz[i]);
            const double D = ncsquare( pd ) - psq_mr2;
            if ( D < 0 ) {
              tgt[i] = -1.0;
            } else {
              const double sqrtD = std::sqrt( D );
              const double tmp = -(sqrtD+pd);
              tgt[i] = tmp > 0.0 ? tmp : -1.0;
            }
          }
        }
      }

      void distToVolumeExitImpl( const double * ncrestrict x,
                                 const double * ncrestrict y,
                                 const double * ncrestrict z,
                                 const double * ncrestrict ux,
                                 const double * ncrestrict uy,
                                 const double * ncrestrict uz,
                                 double * ncrestrict tgt,
                                 std::size_t n ) const ncnoexceptndebug
      {
        //TODO: Properly test this for auto-vectorisation (at least check with
        //-fopt-info-vec-missed).
#ifndef NDEBUG
        for ( std::size_t i = 0; i < n; ++i ) {
          nc_assert( std::isfinite(x[i]) && std::isfinite(y[i]) && std::isfinite(z[i]) );
          nc_assert( std::isfinite(ux[i]) && std::isfinite(uy[i]) && std::isfinite(uz[i]) );
          nc_assert( ncabs(ncsquare(ux[i])+ncsquare(uy[i])+ncsquare(uz[i]) - 1.0) < 1e-13 );
        }
#endif
        //Split into a pre-loop, using tgt as temporary cache area, for
        //efficient loop auto-vectorisation:
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] = m_radiusSq - (x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);//-psq_mr2 = r^2-p^2, since inside tgt[i]>=0

        double buf_pd[NeutronBasket::N];
        for ( std::size_t i = 0; i < n; ++i )
          buf_pd[i] = x[i]*ux[i]+y[i]*uy[i]+z[i]*uz[i];

        //discriminant (should be non-negative since we can't miss the sphere from inside it.
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] += buf_pd[i]*buf_pd[i];//tgt[i] = pd^2-psq_mr2 = D

        for ( std::size_t i = 0; i < n; ++i ) {
          //We should be inside when this function is called! So any negative
          //distances should arise from FP instability only, and we must be on
          //the edge.
          const double sqrtD = fast_sqrt_clippos(tgt[i]);
          tgt[i] = ncmax(0.0,sqrtD-buf_pd[i]);
        }
      }

    };

  }
}

#endif
