#ifndef NCrystal_MMC_Box_hh
#define NCrystal_MMC_Box_hh

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

    class Box {
      double m_dx;
      double m_dy;
      double m_dz;
    public:
      Box( Length dx, Length dy, Length dz )
        : m_dx(dx.get()), m_dy(dy.get()), m_dz(dz.get())
      {
        static_assert( Length::meter == 1.0, "" );
        nc_assert_always( dx.dbl() > 0.0 );
        nc_assert_always( dy.dbl() > 0.0 );
        nc_assert_always( dz.dbl() > 0.0 );
        nc_assert_always( dx.dbl() < 1e199 );
        nc_assert_always( dy.dbl() < 1e199 );
        nc_assert_always( dz.dbl() < 1e199 );
      }

      void toJSONDecodedItems(std::ostream& os) const
      {
        os << "\"name\":\"box\",\"dx\":";
        streamJSON(os,m_dx);
        os << ",\"dy\":";
        streamJSON(os,m_dy);
        os << ",\"dz\":";
        streamJSON(os,m_dz);
      }

      void toCfgString(std::ostream& os) const
      {
        os << "box;dx="<<fmt(m_dx);
        os << ";dy="<<fmt(m_dy);
        os << ";dz="<<fmt(m_dz);
      }

      bool pointIsInside( const Vector& v ) const
      {

        return ( ncabs(v.x()) <= m_dx
                 && ncabs(v.y()) <= m_dy
                 && ncabs(v.z()) <= m_dz );
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
        //The box is an intersection between 3 axis-aligned slabs.

        //"slab method".
        double t1[basket_N];
        double t2[basket_N];

        //use lowest() (~=-1.8e308) as a magic t-value to indicate box misses:
        constexpr double miss_value = std::numeric_limits<double>::lowest();

        //parameterise rays as (x+t*ux,y+t*uy,z+t*uz) and find first range of
        //[t1,t2] values intersecting the x-slab (handling ux=0 appropriately):
        for ( std::size_t i = 0; i < n; ++i ) {
          if ( ux[i] ) nclikely {
            const double invspeed(1.0/ux[i]);
            double tmp1 = (m_dx-x[i])*invspeed;
            double tmp2 = (-m_dx-x[i])*invspeed;
            t1[i] = ncmin( tmp1, tmp2 );
            t2[i] = ncmax( tmp1, tmp2 );
          } else {
            if ( ncabs(x[i]) > m_dx ) {
              //miss the slab AND box completely! Use t1==t2==miss_value to
              //indicate this.
              t1[i] = t2[i] = miss_value;
            } else {
              //entire ray is inside this slab:
              t1[i] = -kInfinity;
              t2[i] = kInfinity;
            }
          }
        }
        //Now look for intersections with the y-slab:
        for ( std::size_t i = 0; i < n; ++i ) {
          if ( t1[i]==miss_value && t2[i]==miss_value )
            continue;
          if ( uy[i] ) nclikely {
            const double invspeed(1.0/uy[i]);
            double tmp1 = (m_dy-y[i])*invspeed;
            double tmp2 = (-m_dy-y[i])*invspeed;
            t1[i] = ncmax( t1[i], ncmin(tmp1, tmp2) );
            t2[i] = ncmin( t2[i], ncmax(tmp1, tmp2) );
          } else {
            if ( ncabs(y[i]) > m_dy ) {
              //miss the slab AND box completely! Use t1==t2==miss_value to
              //indicate this.
              t1[i] = t2[i] = miss_value;
            } else {
              //entire ray is inside this slab, keep previous values unchanged.
            }
          }
        }
        //Now look for intersections with the z-slab:
        for ( std::size_t i = 0; i < n; ++i ) {
          if ( t1[i]==miss_value && t2[i]==miss_value )
            continue;
          if ( uz[i] ) nclikely {
            const double invspeed(1.0/uz[i]);
            double tmp1 = (m_dz-z[i])*invspeed;
            double tmp2 = (-m_dz-z[i])*invspeed;
            t1[i] = ncmax( t1[i], ncmin(tmp1, tmp2) );
            t2[i] = ncmin( t2[i], ncmax(tmp1, tmp2) );
          } else {
            if ( ncabs(z[i]) > m_dz ) {
              //miss the slab AND box completely! Use t1==t2==miss_value to
              //indicate this.
              t1[i] = t2[i] = miss_value;
            } else {
              //entire ray is inside this slab, keep previous values unchanged.
            }
          }
        }
        //Final results:
        for ( std::size_t i = 0; i < n; ++i ) {
          //if t1>=t2: this marks a miss, or a 0-length pathlength through the
          //surface => miss.
          //
          //if t2 == 0, we must have t1<0, i.e. we are on the surface but
          //heading out => miss.
          //
          //if t2<0 both intersections are behind us => miss.
          //
          //Otherwise we have a hit. If t1>=0 take that, otherwise take t2.
          //
          //But if we have a hit, we should also check if we are already
          //inside. That is the case if t1<0 and t2>0.
          if ( t1[i] >= t2[i] || t2[i] <= 0.0 ) {
            tgt[i] = -1.0;//miss
          } else {
            tgt[i] = ncmax( 0.0, t1[i] );
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
        Utils::distToSlabExit( x, ux, tgt, n, m_dx );
        double t2[basket_N];
        Utils::distToSlabExit( y, uy, t2, n, m_dy );
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] = ncmin( tgt[i], t2[i] );
        Utils::distToSlabExit( z, uz, t2, n, m_dz );
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] = ncmin( tgt[i], t2[i] );
      }

    };

  }
}

#endif
