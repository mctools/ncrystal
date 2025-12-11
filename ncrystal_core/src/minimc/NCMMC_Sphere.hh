#ifndef NCrystal_MMC_Sphere_hh
#define NCrystal_MMC_Sphere_hh

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

#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
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
      return std::sqrt(ncmax(0.0,x));
#endif
    }
  }

  namespace MiniMC {

    class Sphere {
      //double m_radius;
      double m_radiusSq;
      double m_input_radius_m;
    public:
      Sphere( Length radius )
        : m_radiusSq( ncsquare( radius.dbl() ) ),
          m_input_radius_m( radius.dbl() )
      {
        static_assert( Length::meter == 1.0, "" );
        nc_assert_always( radius.dbl() > 0.0 );
        nc_assert_always( m_radiusSq < 1e199 );
        nc_assert_always( m_radiusSq > 0.0 );
      }

      void toJSONDecodedItems(std::ostream& os) const
      {
        os << "\"name\":\"sphere\",\"r\":";
        streamJSON(os,m_input_radius_m);
      }

      void toCfgString(std::ostream& os) const
      {
        os << "sphere;r="<<fmt(m_input_radius_m);
      }

      bool pointIsInside( const Vector& v ) const
      {
        return v.mag2() <= m_radiusSq;
      }

      void distToVolumeEntry( const NeutronBasket& nb,
                              Span<double> tgt,
                              std::size_t offset ) const
      {
        nc_assert( tgt.size() >= nb.nused);
        nc_assert( offset < nb.nused);
        nc_assert( nb.nused - offset > 0 );
        distToVolumeEntryImpl( nb.x + offset, nb.y + offset, nb.z + offset,
                               nb.ux + offset, nb.uy + offset, nb.uz + offset,
                               tgt.data() + offset, nb.nused - offset );
      }

      void distToVolumeExit( const NeutronBasket& nb,
                             Span<double> tgt ) const
      {
        nc_assert( tgt.size() >= nb.nused);
        distToVolumeExitImpl( nb.x, nb.y, nb.z,
                              nb.ux, nb.uy, nb.uz,
                              tgt.data(), nb.nused );
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
#ifndef NDEBUG
        for ( std::size_t i = 0; i < n; ++i ) {
          nc_assert( std::isfinite(x[i]) && std::isfinite(y[i]) && std::isfinite(z[i]) );
          nc_assert( std::isfinite(ux[i]) && std::isfinite(uy[i]) && std::isfinite(uz[i]) );
          nc_assert( ncabs(ncsquare(ux[i])+ncsquare(uy[i])+ncsquare(uz[i]) - 1.0) < 1e-13 );
        }
#endif
        double pdotu[basket_N];
        for ( std::size_t i = 0; i < n; ++i )
          pdotu[i] = x[i] * ux[i];
        for ( std::size_t i = 0; i < n; ++i )
          pdotu[i] += y[i] * uy[i];
        for ( std::size_t i = 0; i < n; ++i )
          pdotu[i] += z[i] * uz[i];

        double psq_mr2[basket_N];
        for ( std::size_t i = 0; i < n; ++i )
          psq_mr2[i] = x[i] * x[i];
        for ( std::size_t i = 0; i < n; ++i )
          psq_mr2[i] += y[i] * y[i];
        for ( std::size_t i = 0; i < n; ++i )
          psq_mr2[i] += z[i] * z[i];
        for ( std::size_t i = 0; i < n; ++i )
          psq_mr2[i] -= m_radiusSq;

        for ( std::size_t i = 0; i < n; ++i ) {
          if ( psq_mr2[i] <= 0.0 ) {
            //Inside or on edge. If on edge and normal is not inwards, we should return
            //-1. Otherwise 0.:
            tgt[i] = ( ncmin(psq_mr2[i],pdotu[i])<0.0 ? 0.0 : -1.0 );
          } else {
            const double D = ncsquare( pdotu[i] ) - psq_mr2[i];
            if ( D < 0 ) {
              tgt[i] = -1.0;
            } else {
              const double sqrtD = std::sqrt( D );
              const double tmp = -(sqrtD+pdotu[i]);
              tgt[i] = tmp >= 0.0 ? tmp : -1.0;
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
#ifndef NDEBUG
        for ( std::size_t i = 0; i < n; ++i ) {
          nc_assert( std::isfinite(x[i]) && std::isfinite(y[i]) && std::isfinite(z[i]) );
          nc_assert( std::isfinite(ux[i]) && std::isfinite(uy[i]) && std::isfinite(uz[i]) );
          nc_assert( ncabs(ncsquare(ux[i])+ncsquare(uy[i])+ncsquare(uz[i]) - 1.0) < 1e-13 );
        }
#endif
        //Split into a pre-loop, using tgt as temporary cache area, for
        //efficient loop auto-vectorisation:
        //( -psq_mr2 = r^2-p^2, since inside tgt[i]>=0 )
        //tgt[i] = m_radiusSq - (x[i]*x[i]+y[i]*y[i]+z[i]*z[i]) :
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] = -x[i]*x[i];
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] -= y[i]*y[i];
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] -= z[i]*z[i];
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] += m_radiusSq;

        double buf_pd[NeutronBasket::N];
        //buf_pd[i] = x[i]*ux[i]+y[i]*uy[i]+z[i]*uz[i];
        for ( std::size_t i = 0; i < n; ++i )
          buf_pd[i] = x[i]*ux[i];
        for ( std::size_t i = 0; i < n; ++i )
          buf_pd[i] += y[i]*uy[i];
        for ( std::size_t i = 0; i < n; ++i )
          buf_pd[i] += z[i]*uz[i];

        //discriminant (should be non-negative since we can't miss the sphere
        //from inside it).
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] += buf_pd[i]*buf_pd[i];//tgt[i] = pd^2-psq_mr2 = D

        //We should be inside when this function is called! So any negative
        //distances should arise from FP instability only, and we must be on
        //the edge.
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] = fast_sqrt_clippos(tgt[i]);//sqrtD
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] -= buf_pd[i];
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] = ncmax(0.0,tgt[i]);
      }

    };

  }
}

#endif
