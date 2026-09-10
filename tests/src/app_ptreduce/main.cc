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

#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <iostream>

namespace NC=NCrystal;

#define REQUIRE(x) nc_assert_always(x)
#define REQUIREFLTEQ(x, y) nc_assert_always(::NC::floateq((x), (y)))

namespace {

  // Exercises reduction on small, degenerate, peaked, tail-heavy, and large
  // distributions. Each result is checked for size, ordering, endpoint
  // preservation, finiteness, and exact correspondence with an input point.
  // It also checks the no-op path when targetN is at least the input size.
  // The tests use deterministic data so failures are reproducible.
  void testReducePtsInDistribution()
  {
    using NC::VectD;
    using NC::vectAt;
    using NC::floateq;

#ifndef NDEBUG
    {
      nc_assert(std::isfinite(0.15));
      nc_assert(std::isfinite(1.0e-50));
      nc_assert(0.15 >= 0.0);
      nc_assert(1.0e-50 > 0.0);
    }
#endif

    auto check = [](const VectD& x, const VectD& y, std::size_t k) {
#ifndef NDEBUG
      {
        nc_assert(x.size() == y.size());
        nc_assert(x.size() >= 2);
        nc_assert(k >= 2);
        nc_assert(nc_is_grid(x));

        for (std::size_t i = 0; i < y.size(); ++i) {
          nc_assert(std::isfinite(vectAt(y, i)));
          nc_assert(vectAt(y, i) >= 0.0);
        }
      }
#endif

      const auto r = NC::reducePtsInDistribution( x, y, k );

      REQUIRE(r.first.size() == NC::ncmin(k, x.size()));
      REQUIRE(r.second.size() == r.first.size());
      REQUIRE(r.first.size() >= 2);
      REQUIRE(floateq(r.first.front(), x.front()));
      REQUIRE(floateq(r.first.back(), x.back()));
      REQUIRE(floateq(r.second.front(), y.front()));
      REQUIRE(floateq(r.second.back(), y.back()));

      for (std::size_t i = 1; i < r.first.size(); ++i) {
        REQUIRE(vectAt(r.first, i) > vectAt(r.first, i - 1));
      }

      std::size_t j = 0;
      for (std::size_t i = 0; i < r.first.size(); ++i) {
        while (j < x.size() &&
               (vectAt(x, j) != vectAt(r.first, i) ||
                vectAt(y, j) != vectAt(r.second, i)))
          ++j;

        REQUIRE(j < x.size());
        ++j;
      }
    };

    {
      VectD x{0.0, 1.0};
      VectD y{3.0, 3.0};
      check(x, y, 2);
      check(x, y, 7);
    }

    {
      VectD x{0.0, 1.0, 2.0};
      VectD y{1.0, 100.0, 1.0};
      check(x, y, 2);
    }

    {
      VectD x, y;
      for (std::size_t i = 0; i < 129; ++i) {
        x.push_back(double(i));
        y.push_back(1.0);
      }

      check(x, y, 2);
      check(x, y, 7);
      check(x, y, 32);
      check(x, y, x.size());
    }

    {
      VectD x, y;
      for (std::size_t i = 0; i < 257; ++i) {
        const double t = double(i);
        x.push_back(t);
        y.push_back(
                    1.0 +
                    100.0 * std::exp(-0.5 * (t - 42.0) * (t - 42.0) / 9.0) +
                    40.0 * std::exp(-0.5 * (t - 190.0) * (t - 190.0) / 16.0));
      }

      check(x, y, 2);
      check(x, y, 8);
      check(x, y, 32);
      check(x, y, 64);
    }

    {
      VectD x, y;
      for (std::size_t i = 0; i < 513; ++i) {
        x.push_back(double(i));

        double v = 1.0e-300;
        if (i == 37)
          v = 1.0;
        if (i == 401)
          v = 1.0e-12;

        y.push_back(v);
      }

      check(x, y, 2);
      check(x, y, 16);
      check(x, y, 64);
      check(x, y, 128);
    }

    {
      VectD x, y;
      for (std::size_t i = 0; i < 1001; ++i) {
        const double t = double(i);
        x.push_back(t);
        y.push_back(
                    0.5 +
                    0.2 * std::sin(0.07 * t) +
                    3.0 * std::exp(-0.5 * (t - 230.0) * (t - 230.0) / 25.0) +
                    8.0 * std::exp(-0.5 * (t - 760.0) * (t - 760.0) / 36.0));
      }

      check(x, y, 2);
      check(x, y, 10);
      check(x, y, 50);
      check(x, y, 100);
      check(x, y, 250);
    }

    {
      VectD x, y;
      for (std::size_t i = 0; i < 65; ++i) {
        x.push_back(double(i));
        y.push_back(i % 2 == 0 ? 1.0 : 1.0e-200);
      }

      check(x, y, 2);
      check(x, y, 8);
      check(x, y, 16);
    }
  }
}

int main() {
  std::cout<<"pt reduction testing start..."<<std::endl;
  testReducePtsInDistribution();
  std::cout<<"pt reduction testing done."<<std::endl;
  return 0;
}
