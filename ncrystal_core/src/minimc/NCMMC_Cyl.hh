#ifndef NCrystal_MMC_Cyl_hh
#define NCrystal_MMC_Cyl_hh

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

#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/minimc/NCMMC_Basket.hh"
#include "NCrystal/internal/minimc/NCMMC_Utils.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    //Unbounded or bounded cylinder around axis (0,1,0) and with radius r,
    //containing all points with x^2+z^2<=r*2. If bounded, an additional
    //constraint is added: |y|<dy. The special value dy=0 is used to indicate an
    //unbounded cylinder.

    class Cyl {
      double m_radiusSq;
      double m_input_radius_m;
      double m_dy;
    public:
      Cyl( Length radius, Length dy )
        : m_radiusSq( ncsquare( radius.dbl() ) ),
          m_input_radius_m( radius.dbl() ),
          m_dy( dy.dbl() )
      {
        static_assert( Length::meter == 1.0, "" );
        nc_assert_always( radius.dbl() > 0.0 );
        nc_assert_always( m_radiusSq < 1e199 );
        nc_assert_always( m_radiusSq > 0.0 );
        nc_assert_always( m_dy < 1e199 );
        nc_assert_always( m_dy >= 0.0 );
      }

      void toJSONDecodedItems(std::ostream& os) const
      {
        os << "\"name\":\"cyl\",\"r\":";
        streamJSON(os,m_input_radius_m);
        os << ",\"dy\":";
        streamJSON(os,m_dy);
      }

      void toCfgString(std::ostream& os) const
      {
        //fixme: should we really allow a default value for the radius of
        //cyl/sphere?
        os << "cyl;r="<<fmt(m_input_radius_m);
        if ( m_dy )
          os << ";dy="<<fmt(m_dy);
      }

      bool hasUnboundedDistToVolExit() const
      {
        return m_dy == 0.0;
      }

      //Fixme: would the usage of pointIsInside in the engine (or rather in the
      //constant source) not be better served by simply calling
      //distToVolumeEntry? Then we could also emit errors there and move the
      //particles forward, etc.? The constant source could even move itself
      //towards the geometry (or create a cloned moved copy of itself).
      bool pointIsInside( const Vector& v ) const
      {
        return ( ncsquare(v.x())+ncsquare(v.z()) <= m_radiusSq
                 && ( m_dy==0 || ncabs(v.y()) <= m_dy ) );
      }

      void distToVolumeExit( const NeutronBasket& nb,
                             Span<double> tgt ) const
      {
        nb.validateIfDbg();
        const std::size_t n = nb.nused;
        nc_assert( tgt.size() >= n);
        distToVolumeExitUnboundedImpl( nb.x, nb.y, nb.z,
                                       nb.ux, nb.uy, nb.uz,
                                       tgt.data(), n );
        if ( !m_dy )
          return;
        //check and apply slab limits as well:
        double tslab[basket_N];
        Utils::distToSlabExit( nb.y, nb.uy, tslab, n, m_dy );
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] = ncmin(tgt[i],tslab[i]);
      }

      void distToVolumeEntry( const NeutronBasket& nb,
                              Span<double> tgt,
                              std::size_t offset ) const
      {
        nb.validateIfDbg();
        nc_assert( tgt.size() >= nb.nused);
        nc_assert( offset < nb.nused);
        if ( !m_dy ) {
          distToVolumeEntryUnboundedImpl( nb.x + offset, nb.y + offset,
                                          nb.z + offset, nb.ux + offset,
                                          nb.uy + offset, nb.uz + offset,
                                          tgt.data() + offset, nb.nused-offset );
        } else {
          distToVolumeEntryBoundedImpl( nb.x + offset, nb.y + offset,
                                        nb.z + offset, nb.ux + offset,
                                        nb.uy + offset, nb.uz + offset,
                                        tgt.data() + offset, nb.nused-offset );
        }
      }


    private:
      void validateInputIfDbg( const double * x,
                               const double * y,
                               const double * z,
                               const double * ux,
                               const double * uy,
                               const double * uz,
                               std::size_t n ) const
      {
#ifndef NDEBUG
        for ( std::size_t i = 0; i < n; ++i ) {
          nc_assert( std::isfinite(x[i]) && std::isfinite(y[i]) && std::isfinite(z[i]) );
          nc_assert( std::isfinite(ux[i]) && std::isfinite(uy[i]) && std::isfinite(uz[i]) );
          nc_assert( ncabs(ncsquare(ux[i])+ncsquare(uy[i])+ncsquare(uz[i]) - 1.0) < 1e-13 );
        }
#else
        (void)x;
        (void)y;
        (void)z;
        (void)ux;
        (void)uy;
        (void)uz;
        (void)n;
#endif
      }

      void distToVolumeEntryBoundedImpl( const double * ncrestrict x,
                                         const double * ncrestrict y,
                                         const double * ncrestrict z,
                                         const double * ncrestrict ux,
                                         const double * ncrestrict uy,
                                         const double * ncrestrict uz,
                                         double * ncrestrict tgt,
                                         std::size_t n ) const ncnoexceptndebug//fixme: check all noexcept decorators
      {
        validateInputIfDbg( x, y, z, ux, uy, uz, n );

        // First find intersections with infinite cylinder:
        double twoA[basket_N], B[basket_N], C[basket_N], D[basket_N];
        calcCylIntersectionParams( x, z, ux, uz, twoA, B, C, D, n );
        double tmax[basket_N];//tmax here, using tgt itself for tmin

        for ( std::size_t i = 0; i < n; ++i ) {
          if ( twoA[i] == 0.0 ) {
              //ray parallel to cylinder axis:
            if ( C[i] <= 0 ) {
              //inside, entire ray is inside cylinder
              tgt[i]  = -kInfinity;
              tmax[i] =  kInfinity;
            } else {
              //no intersection, miss:
              tgt[i]=tmax[i]=-1.0;
            }
            continue;
          }
          if ( D[i]<= 0 ) {
            //miss:
            tgt[i] = tmax[i] = -1.0;
            continue;
          }
          const double sqrtD = std::sqrt(D[i]);
          const double invtwoA = 1.0 / twoA[i];
          const double tminval = ( -B[i] - sqrtD )*invtwoA;
          const double tmaxval = ( -B[i] + sqrtD )*invtwoA;
          //Constrain to forward direction only (this might leave tmin=tmax=0,
          //i.e. a miss):
          //fixme: maybe we do the ncmax thing later on anyway???
          tgt[i] = ncmax(0.0,tminval);
          tmax[i] = ncmax(0.0,tmaxval);
        }

        //Now we need to also handle the slab intersections, by narrowing the
        //above intersection:
        for ( std::size_t i = 0; i < n; ++i ) {
          if ( ncabs(y[i])==m_dy ) {
            //edge
            if ( y[i]*uy[i] > 0.0 ) {
              //leaving edge => miss (fixme: this check likely not needed?
              tgt[i] = tmax[i] = -1.;
              continue;
            }
          }
          if ( uy[i] == 0.0 ) {
            //parallel to slab plane
            if ( ncabs(y[i]) > m_dy ) {
              //no intersection possible
              tgt[i] = tmax[i] = -1.;
            } else {
              //entirely in plane, no additional narrowing needed.
            }
            continue;
          }
          //Find actual intersection, and constrain to positive values:
          const double t1 = -(y[i]+m_dy)/uy[i];
          const double t2 = (m_dy-y[i])/uy[i];
          tgt[i] = ncmax(0.0,ncmax(ncmin( t1, t2 ),tgt[i]));
          tmax[i] = ncmax(0.0,ncmin(ncmax( t1, t2 ),tmax[i]));
        }

        //For regular intersections, tgt[i] is now ready. But for misses or rays
        //already inside, we must convert to -1.0 or 0.0 respectively.
        for ( std::size_t i = 0; i < n; ++i ) {
          if ( !(tmax[i]>tgt[i]) ) {
            //miss
            tgt[i] = -1.0;
          } else if ( tmax[i]>0.0 && tgt[i]==0.0 ) {
            //inside
            tgt[i] = 0.0;
          }
        }
      }

      void distToVolumeEntryUnboundedImpl( const double * ncrestrict x,
                                           const double * ncrestrict y,
                                           const double * ncrestrict z,
                                           const double * ncrestrict ux,
                                           const double * ncrestrict uy,
                                           const double * ncrestrict uz,
                                           double * ncrestrict tgt,
                                           std::size_t n ) const ncnoexceptndebug//fixme: check all noexcept decorators
      {
        validateInputIfDbg( x, y, z, ux, uy, uz, n );

        double twoA[basket_N], B[basket_N], C[basket_N];
        //Note: using the tgt array for the discriminator D, to save a local
        //array.
        calcCylIntersectionParams( x, z, ux, uz, twoA, B, C, tgt, n );

        for ( std::size_t i = 0; i < n; ++i ) {
          if ( twoA[i] == 0.0 ) {
              //ray parallel to cylinder axis:
            tgt[i] = ( C[i] <= 0.0 ? 0.0 : -1.0 );
            continue;
          }
          if ( C[i] <= 0.0 ) {
              //Inside or on edge. If on edge and normal is not inwards, we
              //should return -1. Otherwise 0.:
            tgt[i] = ( ( C[i]==0.0 && B[i] >= 0.0 )
                       ? -1.0//on edge and headed out
                       : 0.0//not on edge, or headed in
                       );
            continue;
          }

          if ( tgt[i] <= 0 ) {//Remember tgt[i] holds the discriminant D
            tgt[i] = -1.0;
            continue;
          }
          const double sqrtD = std::sqrt(tgt[i]);//Remember tgt[i] holds
                                                 //the discriminant D

          const double tmin = ( -B[i] - sqrtD );
          const double tmax = ( -B[i] + sqrtD );
          if ( tmax < 0.0 ) {
            tgt[i] = -1.0;//intersection in backwards direction, miss
          } else if ( tmin > 0.0 ) {
            tgt[i] = tmin / twoA[i];//intersects twice
          } else {
            //Presumably we are on the edge if we can get here. Repeat edge
            //trick from above for safety:
            tgt[i] = (  B[i] < 0.0 ? 0.0 : -1.0 );
          }
        }
      }

      void calcCylIntersectionParams( const double * ncrestrict x,
                                      const double * ncrestrict z,
                                      const double * ncrestrict ux,
                                      const double * ncrestrict uz,
                                      double * ncrestrict twoA,
                                      double * ncrestrict B,
                                      double * ncrestrict C,
                                      double * ncrestrict D,
                                      std::size_t n ) const ncnoexceptndebug
      {
        //The condition |(x,z)+t*(ux,yz)|^2=r^2 leads to a second degree
        //equation for t: A*t^2+B*t+C = 0. This helper function calculates these
        //coefficients (2A,B,C,D=B^2-4AC)
        //
        //Note that: A = ux^2+uz^2 = 1-uy^2 (A=0 means ray parallel to cyl-axis)
        //           B/2 = x*ux + z*uz (whether or not direction is into
        //                              or out of cylinder)
        //           C = x^2+z^2-r^2 (so C<0 means "inside cylinder")
        //           D < 0 : no intersection
#ifndef NDEBUG
        for ( std::size_t i = 0; i < n; ++i ) {
          nc_assert( std::isfinite(x[i]) && std::isfinite(z[i]) );
          nc_assert( std::isfinite(ux[i]) && std::isfinite(uz[i]) );
          nc_assert( ncsquare(ux[i])+ncsquare(uz[i]) <= 1.0+1e-13 );
        }
#endif
        for ( std::size_t i = 0; i < n; ++i )
          twoA[i] = ux[i]*ux[i];
        for ( std::size_t i = 0; i < n; ++i )
          twoA[i] += uz[i]*uz[i];
        for ( std::size_t i = 0; i < n; ++i )
          twoA[i] *= 2;
        for ( std::size_t i = 0; i < n; ++i )
          B[i] = x[i]*ux[i];
        for ( std::size_t i = 0; i < n; ++i )
          B[i] += z[i]*uz[i];
        for ( std::size_t i = 0; i < n; ++i )
          B[i] *= 2;
        for ( std::size_t i = 0; i < n; ++i )
          C[i] = x[i]*x[i];
        for ( std::size_t i = 0; i < n; ++i )
          C[i] += z[i]*z[i];
        for ( std::size_t i = 0; i < n; ++i )
          C[i] -= m_radiusSq;
        for ( std::size_t i = 0; i < n; ++i )
          D[i] = twoA[i]*C[i];
        for ( std::size_t i = 0; i < n; ++i )
          D[i] *= -2.0;
        for ( std::size_t i = 0; i < n; ++i )
          D[i] += B[i]*B[i];
      }

      void distToVolumeExitUnboundedImpl( const double * ncrestrict x,
                                          const double * ncrestrict y,
                                          const double * ncrestrict z,
                                          const double * ncrestrict ux,
                                          const double * ncrestrict uy,
                                          const double * ncrestrict uz,
                                          double * ncrestrict tgt,
                                          std::size_t n ) const ncnoexceptndebug
      {
        validateInputIfDbg( x, y, z, ux, uy, uz, n );
#ifndef NDEBUG
        for ( std::size_t i = 0; i < n; ++i ) {
          nc_assert( x[i]*x[i]+z[i]*z[i] <= m_radiusSq*(1.0+1e-9) );
          nc_assert( m_dy==0 || ncabs(y[i]) <= m_dy*(1.0+1e-9) );
        }
#endif

        //To save place, using the tgt array for the D (discriminator) parameter:
        double twoA[basket_N], B[basket_N], C[basket_N];
        calcCylIntersectionParams( x, z, ux, uz, twoA, B, C, tgt, n );

        //tgt array now holds "D". Convert it into sqrt(D)-B:

        for ( std::size_t i = 0; i < n; ++i ) {
          //UB if not inside the volume, so we should always have solutions. The
          //abs call is just to avoid FPE in case of numerical fluctuations near
          //zero.
          tgt[i] = ncabs(tgt[i]);
        }
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] = std::sqrt(tgt[i]);
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] -= B[i];
        //Also clamp to [0,inf], since we know we are inside the cylinder, so
        //tmax>=0 apart from numerical fluctuations:
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] = ncmax(0.0,tgt[i]);

        //Final loop, accounting for degeneracies (unlikely to be vectorizable
        //due to branches):
        for ( std::size_t i = 0; i < n; ++i ) {
          nc_assert( C[i] <= 1e-9*m_input_radius_m );//we are inside
          C[i] = ncmin(C[i],0.0);//fluctuation guard
          nc_assert( !ncisnan(twoA[i]*C[i]) );
          if ( twoA[i]*C[i] == 0.0 ) ncunlikely {
            if ( !twoA[i] ) {
              tgt[i] = kInfinity;
              continue;
            } else if ( B[i] >= 0.0 ) {
              tgt[i] = 0.0;//edge, ray headed out
              continue;
            } else {
              //C[i] = 0 but ray headed in, proceed with standard code
            }
          }
          //Regular case! Remember tgt already holds sqrt(D)-B:
          nc_assert( twoA[i] > 0 );
          tgt[i] /= twoA[i];
        }
      }
    };

  }
}

#endif
