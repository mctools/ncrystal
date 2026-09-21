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

////////////////////////////////////////////////////////////////////////////////
// Tests of VDOS::mergeGridsWithTol, which merges the points of two grids     //
// (used to add points needed for the low energy limit to the beta and alpha  //
// grids of a scattering kernel), including that the result is not affected   //
// by rounding errors in the position of points which mathematically are the  //
// same in the two grids.                                                     //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/vdos/NCVDOSUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>
#include <cstdint>
#include <algorithm>

namespace NC=NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  void testMergeGridsWithTol()
  {
    using NC::VectD;
    auto merge = []( const VectD& a, const VectD& b, double rtol = 0.01 )
    {
      auto r = NC::VDOS::mergeGridsWithTol( a, b, rtol );
      REQUIRE( NC::nc_is_grid( r ) );
      //All points of a are always in the result:
      for ( double v : a )
        REQUIRE( std::binary_search( r.begin(), r.end(), v ) );
      //Anything else in the result comes from b:
      for ( double v : r )
        REQUIRE( std::binary_search( a.begin(), a.end(), v )
                 || std::binary_search( b.begin(), b.end(), v ) );
      //and the extremes of both grids are always kept:
      REQUIRE( r.front() == NC::ncmin( a.front(), b.front() ) );
      REQUIRE( r.back() == NC::ncmax( a.back(), b.back() ) );
      return r;
    };

    //Points which are far from the neighbours in a are all kept:
    {
      auto r = merge( VectD{ 1.0, 2.0, 4.0, 8.0 }, VectD{ 1.5, 3.0, 6.0 } );
      REQUIRE( r == ( VectD{ 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0 } ) );
    }
    //Merging with an identical grid, or with a subset of the grid, adds nothing:
    {
      VectD a{ -3.0, -0.5, 0.0, 0.25, 1.0, 7.0 };
      REQUIRE( merge( a, a ) == a );
      REQUIRE( merge( a, VectD{ -0.5, 1.0 } ) == a );
      REQUIRE( merge( a, VectD{ 0.0, 0.25 } ) == a );
    }
    //A point in b which is close (within 1%) to both of its neighbours in a is
    //not kept, and it is kept when it is not (but note that it takes only being
    //far enough from one of the neighbours, see documentation):
    {
      VectD a{ 1.0, 1.006, 3.0 };
      REQUIRE( merge( a, VectD{ 1.003, 2.0 } ) == ( VectD{ 1.0, 1.006, 2.0, 3.0 } ) );
      REQUIRE( merge( a, VectD{ 1.003, 2.0 }, 0.001 )
               == ( VectD{ 1.0, 1.003, 1.006, 2.0, 3.0 } ) );
    }
    //Extremes of b outside the range of a are kept, even if close to a:
    {
      VectD a{ 1.0, 2.0 };
      REQUIRE( merge( a, VectD{ 0.9999, 1.5 } ) == ( VectD{ 0.9999, 1.0, 1.5, 2.0 } ) );
      REQUIRE( merge( a, VectD{ 1.5, 2.0001 } ) == ( VectD{ 1.0, 1.5, 2.0, 2.0001 } ) );
    }
    //Zeros and signs: zero is a single point, and points of different sign are
    //always far enough from each other, and from zero:
    {
      VectD a{ -1.0, 0.0, 1.0 };
      REQUIRE( merge( a, VectD{ 0.0, 1.0 } ) == a );
      REQUIRE( merge( a, VectD{ -1e-9, 0.0, 1e-9 } ) == ( VectD{ -1.0, -1e-9, 0.0, 1e-9, 1.0 } ) );
      //-0.999 and 0.999 are within 1% of a point in a, but far enough from the
      //other neighbour (0.0), so they are kept:
      VectD b{ -2.0, -0.999, -0.5, 0.999, 2.0 };
      auto r = merge( a, b );
      REQUIRE( r == ( VectD{ -2.0, -1.0, -0.999, -0.5, 0.0, 0.999, 1.0, 2.0 } ) );
    }
    //Mirror symmetry: merging negated and reversed grids gives the negated
    //and reversed result (since the rule is symmetric in the sign):
    {
      VectD a{ -6.0, -2.5, -0.01, 0.0, 0.3, 1.0, 5.0 };
      VectD b{ -6.5, -2.51, -1.0, -0.011, 0.001, 0.31, 0.5, 1.004, 4.0 };
      auto neg = []( const VectD& v ) {
        VectD r;
        for ( auto it = v.rbegin(); it != v.rend(); ++it )
          r.push_back( -*it );
        return r;
      };
      REQUIRE( neg( merge( a, b ) ) == merge( neg( a ), neg( b ) ) );
    }

    //The main point of the following tests: the same node of a lattice can be
    //calculated in different ways for the two grids, so points which
    //mathematically are identical can be bitwise identical, or differ by one or
    //a few ulps. This must not make any difference (previously, when the two
    //values were not bitwise identical, both would be kept if the point was
    //far away from the other neighbour in a).
    {
      VectD a{ 0.5, 0.6789802565272282, 0.9 };
      //a single point in b, which is an ulp from a point in a:
      for ( double x : { 0.6789802565272282, std::nextafter( 0.6789802565272282, 0.0 ),
                         std::nextafter( 0.6789802565272282, 1.0 ) } ) {
        REQUIRE( merge( a, VectD{ 0.5, x } ) == a );
      }
      //some ulps and up to 1e-12 (relative) from a point in a:
      double x = 0.6789802565272282;
      for ( int i = 0; i < 40; ++i ) {
        x = std::nextafter( x, 0.0 );
        REQUIRE( merge( a, VectD{ 0.5, x } ) == a );
      }
      REQUIRE( merge( a, VectD{ 0.5, 0.6789802565272282*(1.0-1e-12) } ) == a );
      REQUIRE( merge( a, VectD{ 0.5, 0.6789802565272282*(1.0+1e-12) } ) == a );
      //Same for negative values, and below/above the whole range of a:
      VectD an{ -0.9, -0.6789802565272282, -0.5 };
      REQUIRE( merge( an, VectD{ std::nextafter( -0.6789802565272282, 0.0 ), -0.5 } ) == an );
      REQUIRE( merge( an, VectD{ std::nextafter( -0.6789802565272282, -1.0 ), -0.5 } ) == an );
      //Clearly different values (here 1e-6 relative, still far below the
      //tolerance rtol so they are considered too close in the ratio test
      //to be "far" from that neighbour, but far from the other neighbour) are
      //kept as before:
      REQUIRE( merge( a, VectD{ 0.5, 0.6789802565272282*(1.0+1e-6) } ).size() == 4 );
    }

    //Whole grids: b is the same lattice as a (or a part of it), but with
    //nodes calculated differently, so with independent relative noise at the
    //level of rounding errors (up to 50 ulps). The result must then not depend
    //on the noise: same number of points, and the same points as without it.
    {
      struct Rng {
        std::uint64_t s = 88172645463325252ULL;
        double u()
        {
          s ^= s << 13;
          s ^= s >> 7;
          s ^= s << 17;
          return static_cast<double>( s % 2000001 ) / 1000000.0 - 1.0;
        }
      } rng;
      for ( int trial = 0; trial < 200; ++trial ) {
        //lattice with a few different spacings, containing 0 and negative
        //values, and with coarse and fine parts:
        VectD lat;
        const double dx = 0.01 + 0.03 * ( 1.0 + rng.u() );
        for ( int i = -80; i <= 120; ++i )
          lat.push_back( dx * i );
        VectD fine;
        for ( int i = 60; i <= 100; ++i )
          fine.push_back( dx * i / 8.0 );
        VectD a;//coarse lattice, without every 5th node
        for ( std::size_t i = 0; i < lat.size(); ++i )
          if ( i % 5 != 0 )
            a.push_back( lat[i] );
        //b: fine nodes plus every 5th (missing) coarse node, and also (as
        //duplicates of the nodes in a) a subset of the nodes of a:
        VectD bclean;
        for ( std::size_t i = 0; i < lat.size(); ++i )
          if ( i % 5 == 0 || i % 7 == 0 )
            bclean.push_back( lat[i] );
        for ( double v : fine )
          bclean.push_back( v );
        std::sort( bclean.begin(), bclean.end() );
        bclean.erase( std::unique( bclean.begin(), bclean.end() ), bclean.end() );
        //Zero is at a node in both: make sure it is exactly zero in both:
        const auto clean = merge( a, bclean );
        VectD bnoisy;
        for ( double v : bclean ) {
          double w = v;
          if ( v != 0.0 ) {
            const int nulp = static_cast<int>( 50.0 * std::abs( rng.u() ) );
            for ( int k = 0; k < nulp; ++k )
              w = std::nextafter( w, rng.u() < 0.0 ? -1e300 : 1e300 );
          }
          bnoisy.push_back( w );
        }
        REQUIRE( NC::nc_is_grid( bnoisy ) );
        const auto noisy = merge( a, bnoisy );
        REQUIRE( noisy.size() == clean.size() );
        for ( std::size_t i = 0; i < clean.size(); ++i )
          REQUIRE( NC::ncabs( noisy[i] - clean[i] )
                   <= 1e-13 * NC::ncmax( 1.0, NC::ncabs( clean[i] ) ) );
      }
    }
  }


}

int main()
{
  std::cout<<"Testing mergeGridsWithTol"<<std::endl;
  testMergeGridsWithTol();
  std::cout<<"Testing mergeGridsWithTol... done"<<std::endl;
  return 0;
}
