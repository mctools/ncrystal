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

#include "NCrystal/internal/vdos/NCVDOSKnlGrid.hh"
#include "NCrystal/internal/utils/NCSpan.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMsg.hh"//fixme: add lots of debug output?
#include "NCrystal/internal/fact_utils/NCFactoryJobs.hh"//fixme: use multithreaded gn merging (those gn of same binwidth in each thread)
#include <queue>
namespace NC=NCrystal;

namespace NCRYSTAL_NAMESPACE {

  //FIXME: FWD INCLUDE HACK FOR NOW:
  PairDD rangeXNexpMX(unsigned n, double eps, double accuracy = 1e-13 );

  namespace VDOS {
    namespace {

      // evalPWLSum: Evaluates the weighted sum of piecewise-linear functions on
      // the supplied grid.
      //
      // Each function is zero outside its own [x0, xmax()] interval.
      // Values inside an interval are obtained by linear interpolation.
      // `fs` contains the input functions and `ws` their corresponding weights.
      // The returned vector has one value for each point in `grid`.

      struct PWLFct {
        double x0; //x{i=0}
        double binWidth;//distance between x{i} and x{i+1}
        Span<const double> f;//values of f at the x{i} points. The size of the
                             //span encodes the number of points.
        double xmax() const { return x0 + (f.size()-1)*binWidth; }
        VectD dataHolder;//optional, so can hold its data if needed.
        Span<double> f_mutable()
        {
          //only possible when we hold our data
          nc_assert( dataHolder.size() == f.size() );
          return dataHolder;
        }


      };

      PWLFct pwlNarrowToPos(const PWLFct& p, double tol = 1e-3 )
      {
        nc_assert_always( p.binWidth > 0.0 );
        const double t = tol * p.binWidth;//threshold
        std::size_t i = ( p.x0 < t ? static_cast<std::size_t>
                          (std::ceil((t - p.x0) / p.binWidth)) : 0u );
        i = std::min<std::size_t>(i, p.f.size());
        PWLFct res;
        res.x0 = p.x0 + i * p.binWidth;
        res.binWidth = p.binWidth;
        res.dataHolder.assign(p.f.begin() + i, p.f.end());
        res.f = res.dataHolder;
        nc_assert_always( res.x0 >= t );
        nc_assert_always( res.f.size() >= 2 );
        return res;
      }

#if 1
      VectD evalPWLSum(Span<const PWLFct> fs,
                       Span<const double> grid,
                       Span<const double> ws = {})
      {
        VectD out(grid.size(), 0.0);

        const double* ncrestrict gridPtr = grid.data();
        double* ncrestrict outPtr = out.data();

        const bool weighted = !ws.empty();

        for (std::size_t n = 0; n < fs.size(); ++n) {
          const PWLFct& p = fs[n];
          const auto& f = p.f;

          const std::size_t lastBin = f.size() - 1;
          const double invbw = 1.0 / p.binWidth;
          const double weight = weighted ? ws[n] : 1.0;
          const double xmax = p.xmax();

          std::size_t g = static_cast<std::size_t>(
                                                   std::lower_bound(grid.begin(), grid.end(), p.x0)
                                                   - grid.begin());

          if (g == grid.size() || gridPtr[g] > xmax)
            continue;

          for (std::size_t i = 0; i < lastBin && g < grid.size(); ++i) {
            const double xLeft =
              p.x0 + static_cast<double>(i) * p.binWidth;

            const double xRight =
              xLeft + p.binWidth;

            const double y0 = vectAt(f, i);
            const double y1 = vectAt(f, i + 1);

            const double slope =
              (y1 - y0) * invbw;

            const double intercept =
              y0 - slope * xLeft;

            std::size_t end = g;

            while (end < grid.size() &&
                   gridPtr[end] < xRight) {
              ++end;
            }

            // This loop is deliberately simple so the compiler can
            // auto-vectorize it.
            for (std::size_t k = g; k < end; ++k) {
              outPtr[k] += weight *
                (intercept + slope * gridPtr[k]);
            }

            g = end;
          }

          // Handle the final sample at x == xmax.
          if (g < grid.size() && gridPtr[g] <= xmax) {
            const double y = vectAt(f, lastBin);

            std::size_t end = g;

            while (end < grid.size() &&
                   gridPtr[end] <= xmax) {
              ++end;
            }

            for (std::size_t k = g; k < end; ++k)
              outPtr[k] += weight * y;
          }
        }

        return out;
      }




#elif 0

      VectD evalPWLSum(Span<const PWLFct> fs,
                       Span<const double> grid,
                       Span<const double> ws = {})
      {
#ifndef NDEBUG
        nc_assert( ws.empty() || fs.size() == ws.size());

        for (auto& w : ws)
          nc_assert(std::isfinite(w));

        for (auto& p : fs) {
          nc_assert(std::isfinite(p.x0));
          nc_assert(std::isfinite(p.binWidth) && p.binWidth > 0.0);
          nc_assert(p.f.size() >= 2);
          nc_assert(std::isfinite(p.xmax()));

          for (auto& y : p.f)
            nc_assert(std::isfinite(y) && y >= 0.0);
        }

        for (auto& x : grid)
          nc_assert(std::isfinite(x));
#endif

        VectD out(grid.size(), 0.0);

        const bool weighted = !ws.empty();

        for (std::size_t n = 0; n < fs.size(); ++n) {
          const PWLFct& p = fs[n];
          const auto& f = p.f;

          const double xmax = p.xmax();
          const double invbw = 1.0 / p.binWidth;
          const double weight = weighted ? ws[n] : 1.0;

          const auto first =
            std::lower_bound(grid.begin(), grid.end(), p.x0);

          const auto last =
            std::upper_bound(first, grid.end(), xmax);

          const std::size_t begin =
            static_cast<std::size_t>(first - grid.begin());

          const std::size_t end =
            static_cast<std::size_t>(last - grid.begin());

          const std::size_t lastBin = f.size() - 1;

          for (std::size_t g = begin; g < end; ++g) {
            const double x = grid[g];

            const double q = (x - p.x0) * invbw;
            const std::size_t i = static_cast<std::size_t>(q);

            double y;

            if (i >= lastBin) {
              y = vectAt(f, lastBin);
            } else {
              const double t = q - static_cast<double>(i);
              const double y0 = vectAt(f, i);
              const double y1 = vectAt(f, i + 1);

              y = y0 + t * (y1 - y0);
            }

            out[g] += weight * y;
          }
        }

        return out;
      }

#else
      VectD evalPWLSum(Span<const PWLFct> fs,
                       Span<const double> grid,
                       Span<const double> ws = {})
      {
#ifndef NDEBUG
        nc_assert( ws.empty() || fs.size() == ws.size());

        for (auto& w : ws)
          nc_assert(std::isfinite(w));

        for (auto& p : fs) {
          nc_assert(std::isfinite(p.x0));
          nc_assert(std::isfinite(p.binWidth) && p.binWidth > 0.0);
          nc_assert(p.f.size() >= 2);
          nc_assert(std::isfinite(p.xmax()));

          for (auto& y : p.f)
            nc_assert(std::isfinite(y) && y >= 0.0);
        }

        for (auto& x : grid)
          nc_assert(std::isfinite(x));
#endif

        VectD invbws;
        invbws.reserve(fs.size());
        for (auto& p : fs) {
          const double invbw = 1.0 / p.binWidth;
          invbws.push_back(invbw);
        }

        VectD out;
        out.reserve(grid.size());

        for (auto& x : grid) {
          double z = 0.0;

          for (std::size_t n = 0; n < fs.size(); ++n) {
            const auto& p = fs[n];
            if (x < p.x0 || x > p.xmax())
              continue;
            const double q = (x - p.x0) * vectAt(invbws, n);
            const std::size_t i = static_cast<std::size_t>(q);
            const double w = ws.empty() ? 1.0 : ws[n];
            if (i + 1 >= p.f.size()) {
              z += w * p.f.back();
            } else {
              const double t = q - static_cast<double>(i);
              const double y0 = vectAt(p.f, i);
              const double y1 = p.f[i + 1];
              z += w * (y0 + t * (y1 - y0));
            }
          }
          out.push_back(z);
        }

        return out;
      }
#endif 

      // Calculate the weighted sum of several piecewise-linear functions
      // represented by uniformly spaced y-value arrays. Each input function may
      // have a different starting x-value and array length, but all must use
      // the same binwidth and their starting positions must be aligned to
      // integral binwidth offsets.  The output grid covers the combined support
      // of all input functions; values outside an individual function's support
      // do not contribute.
      std::pair<std::vector<double>, double>
      combinePWLFctsOfSameBinWidth( const std::vector<Span<const double>>& yarrays,
                                    Span<const double> x0values,
                                    double binWidth,
                                    Span<const double> weights,
                                    double alignment_tol = 1e-12 )
      {
        (void)alignment_tol;
        const std::size_t N = yarrays.size();
#ifndef NDEBUG
        nc_assert(N > 0);
        nc_assert(x0values.size() == N);
        nc_assert(weights.size() == N);
        nc_assert(std::isfinite(binWidth));
        nc_assert(binWidth > 0.0);
        nc_assert(std::isfinite(alignment_tol));
        nc_assert(alignment_tol >= 0.0);
        for (std::size_t i = 0; i < N; ++i) {
          nc_assert(yarrays[i].size() > 0);
          nc_assert(std::isfinite(x0values[i]));
          nc_assert(std::isfinite(weights[i]));
          for (std::size_t j = 0; j < yarrays[i].size(); ++j) {
            nc_assert(std::isfinite(yarrays[i][j]));
          }
        }
#endif
        // Use the first function's x0 as the reference origin.
        const double x0_reference = x0values[0];
        std::vector<std::int64_t> offsets(N);
        std::int64_t first_offset = 0;
        std::int64_t last_offset = 0;

        for (std::size_t i = 0; i < N; ++i) {
          const double offset = (x0values[i]-x0_reference)/binWidth;

          const std::int64_t rounded_offset
            = static_cast<std::int64_t>(std::llround(offset));

          if ( !( ncabs( offset - static_cast<double>(rounded_offset) )
                  <= static_cast<double>(alignment_tol) )) {
            //fixme
            NCRYSTAL_MSG("TKTEST FAIL offset="<<fmt(offset)<<" rounded: "<<rounded_offset);
            NCRYSTAL_MSG("x0values[i]="<<fmt(x0values[i])<<" x0_reference="<<fmt(x0_reference)
                         <<" binWidth="<<fmt(binWidth));
          }
          nc_assert( ncabs( offset - static_cast<double>(rounded_offset) )
                     <= static_cast<double>(alignment_tol) );

          offsets[i] = rounded_offset;

          if (i == 0) {
            first_offset = rounded_offset;
            last_offset = rounded_offset +
              static_cast<std::int64_t>(yarrays[i].size()) - 1;
          } else {
            if (rounded_offset < first_offset) {
              first_offset = rounded_offset;
            }

            const std::int64_t function_last =
              rounded_offset +
              static_cast<std::int64_t>(yarrays[i].size()) - 1;

            if (function_last > last_offset) {
              last_offset = function_last;
            }
          }
        }

        nc_assert(last_offset >= first_offset);

        const std::int64_t output_size_ll = last_offset - first_offset + 1;

        nc_assert(output_size_ll >= 0);
        nc_assert(output_size_ll <= 100000000);

        const std::size_t output_size =
          static_cast<std::size_t>(output_size_ll);

        std::vector<double> yout(output_size, 0.0);

        // The output grid starts at the smallest input-grid position.
        const double x0out = x0_reference + first_offset*binWidth;
        for (std::size_t i = 0; i < N; ++i) {
          const std::size_t start =
            static_cast<std::size_t>(offsets[i] - first_offset);
          const double weight = weights[i];
          for (std::size_t j = 0; j < yarrays[i].size(); ++j)
            yout[start + j] += weight * yarrays[i][j];
        }
        return std::make_pair(std::move(yout), x0out);
      }

      struct EquidistantGrid final {
        //fixme: just use the PWLFcts struct??
        double x0, binWidth;
        std::size_t npts;
        double x1() const { return x0 + binWidth*(npts-1); }
      };

#if 1

      // Merges several evenly spaced input grids into one sorted grid.
      // Input grids may have different offsets and bin widths.
      // Points closer than 10% of the smallest input bin width are merged.
      // The output is not necessarily evenly spaced.
      // Its first and last points are always the minimum and maximum input ranges.
      // Empty input returns an empty grid.
      VectD makeCommonGrid( Span<const EquidistantGrid> grids )
      {
#ifndef NDEBUG
        for ( auto& g : grids ) {
          nc_assert( std::isfinite(g.x0) );
          nc_assert( std::isfinite(g.binWidth) );
          nc_assert( g.binWidth > 0.0 );
          nc_assert( g.npts >= 2 );
          nc_assert( std::isfinite(g.x1()) );
        }
#endif

        if ( grids.empty() )
          return {};

        if ( grids.size() == 1 ) {
          auto& g = grids.front();
          return linspace(g.x0, g.x1(), g.npts);
        }

        //might prune before the end if memory grows too much:
        constexpr std::size_t maxPending = 200000;
        constexpr std::size_t pruneGap = 10;
        double xmin = kInfinity;
        double xmax = -kInfinity;
        double minbw = kInfinity;

        for ( auto& g : grids ) {
          xmin = ncmin(xmin, g.x0);
          xmax = ncmax(xmax, g.x1());
          minbw = ncmin(minbw, g.binWidth);
        }

        const double atol = 0.1 * minbw;
        nc_assert( atol > 1e-200 );
        VectD pts;
        pts.reserve(maxPending);

        auto prune = [&]()//fixme
        {
          if ( pts.empty() )
            return;

          std::sort(pts.begin(), pts.end());

          std::size_t n = 1;
          vectAt(pts, 0) = xmin;

          for ( std::size_t i = 1; i + 1 < pts.size(); ++i ) {
            const double x = vectAt(pts, i);
            const double prev = vectAt(pts, n - 1);
            if ( ncabs(x - prev) > atol && ncabs(xmax - x) > atol )
              vectAt(pts, n++) = x;
          }
          if ( xmax > xmin )
            vectAt(pts, n++) = xmax;
          pts.resize(n);
        };

        std::size_t sincePrune = 0;

        for ( auto& g : grids ) {
          for ( std::size_t i = 0; i < g.npts; ++i ) {
            const double x = ( i + 1 == g.npts
                               ? g.x1()
                               : g.x0 + g.binWidth * static_cast<double>(i) );
            pts.emplace_back(x);
          }

          ++sincePrune;
          if ( pts.size() > maxPending && sincePrune >= pruneGap ) {
            prune();
            sincePrune = 0;
          }
        }

        if ( sincePrune != 0 )
          prune();

        pts.front() = xmin;
        pts.back() = xmax;
        nc_assert( nc_is_grid(pts) );
        return pts;
      }
#elif 0
      // Merges sorted equidistant grids into one sorted grid.
      // Points within 10% of the smallest input bin width are merged.
      // The output starts at the minimum input point and ends at the maximum.
      // The output need not be equidistant.
      // Empty input produces an empty vector.
      // Each input grid must contain at least two finite points.
      VectD makeCommonGrid( Span<const EquidistantGrid> grids )
      {
#ifndef NDEBUG
        for ( auto& g : grids ) {
          nc_assert( std::isfinite(g.x0) );
          nc_assert( std::isfinite(g.binWidth) );
          nc_assert( g.binWidth > 0.0 );
          nc_assert( g.npts >= 2 );
          nc_assert( std::isfinite(g.x1()) );
        }
#endif

        if ( grids.empty() )
          return {};

        if ( grids.size() == 1 ) {
          auto& g = grids.front();
          return linspace(g.x0, g.x1(), g.npts);
        }

        struct Node {
          double x;
          std::size_t gi;
          std::size_t i;
        };

        struct NodeLess {
          bool operator()( const Node& a, const Node& b ) const
          {
            return a.x > b.x;
          }
        };

        double xmin = kInfinity;
        double xmax = -kInfinity;
        double minbw = kInfinity;

        for ( auto& g : grids ) {
          xmin = ncmin(xmin, g.x0);
          xmax = ncmax(xmax, g.x1());
          minbw = ncmin(minbw, g.binWidth);
        }

        const double atol = 0.1 * minbw;
        nc_assert( atol > 1e-200 );

        std::priority_queue<Node, std::vector<Node>, NodeLess> q;

        for ( std::size_t gi = 0; gi < grids.size(); ++gi )
          q.push({ grids[gi].x0, gi, 0 });

        VectD pts;
        pts.emplace_back(xmin);

        while ( !q.empty() ) {
          const Node p = q.top();
          q.pop();

          if ( p.x - pts.back() > atol && xmax - p.x > atol )
            pts.emplace_back(p.x);

          auto& g = grids[p.gi];
          const std::size_t i = p.i + 1;

          if ( i < g.npts ) {
            const double x = i == g.npts - 1
              ? g.x1()
              : g.x0 + g.binWidth * static_cast<double>(i);
            q.push({ x, p.gi, i });
          }
        }

        pts.emplace_back(xmax);
        nc_assert( nc_is_grid(pts) );
        return pts;
      }



#elif 0
      // Merges sorted equidistant grids into one sorted grid.
      // Points within 10% of the smallest input bin width are merged.
      // The output starts at the minimum input point and ends at the maximum.
      // The output need not be equidistant.
      // Empty input produces an empty vector.
      // Each input grid must contain at least two finite points.
      VectD makeCommonGrid( Span<const EquidistantGrid> grids )
      {
#ifndef NDEBUG
        for ( auto& g : grids ) {
          nc_assert( std::isfinite(g.x0) );
          nc_assert( std::isfinite(g.binWidth) );
          nc_assert( g.binWidth > 0.0 );
          nc_assert( g.npts >= 2 );
          nc_assert( std::isfinite(g.x1()) );
        }
#endif

        constexpr bool do_pre_coalescing = true;

        if ( grids.empty() )
          return {};

        if ( grids.size() == 1 ) {
          auto& g = grids.front();
          return linspace(g.x0, g.x1(), g.npts);
        }

        double xmin = kInfinity;
        double xmax = -kInfinity;
        double minbw = kInfinity;

        for ( auto& g : grids ) {
          xmin = ncmin(xmin, g.x0);
          xmax = ncmax(xmax, g.x1());
          minbw = ncmin(minbw, g.binWidth);
        }

        const double atol = 0.1 * minbw;
        nc_assert( atol > 1e-200 );

        struct Run {
          double x0, x1, bw;
          std::size_t npts;
        };

        VectD pts;
        std::vector<Run> runs;

        if ( do_pre_coalescing ) {
          for ( auto& g : grids ) {
            bool merged = false;

            for ( auto& r : runs ) {
              if ( r.bw != g.binWidth )
                continue;

              const double h = r.bw;
              const double ptol = ncmin(1e-5 * h, atol);
              const double phase = std::remainder(g.x0 - r.x0, h);

              if ( ncabs(phase) > ptol )
                continue;

              const double k = std::round((g.x0 - r.x0) / h);
              const double g0 = k;
              const double g1 = k +
                static_cast<double>(g.npts - 1);
              const double r1 = static_cast<double>(r.npts - 1);

              if ( g0 > r1 + 1.0 || g1 < -1.0 )
                continue;

              r.x0 = ncmin(r.x0, g.x0);
              r.x1 = ncmax(r.x1, g.x1());
              r.npts = static_cast<std::size_t>(
                                                std::round((r.x1 - r.x0) / h)) + 1;
              merged = true;
              break;
            }

            if ( !merged )
              runs.push_back({ g.x0, g.x1(), g.binWidth, g.npts });
          }
        }
        else {
          runs.reserve(grids.size());

          for ( auto& g : grids )
            runs.push_back({ g.x0, g.x1(), g.binWidth, g.npts });
        }

        struct Node {
          double x;
          std::size_t ri, i;
        };

        struct NodeLess {
          bool operator()( const Node& a, const Node& b ) const
          {
            return a.x > b.x;
          }
        };

        std::priority_queue<Node, std::vector<Node>, NodeLess> q;

        for ( std::size_t ri = 0; ri < runs.size(); ++ri )
          q.push({ vectAt(runs, ri).x0, ri, 0 });

        pts.emplace_back(xmin);

        while ( !q.empty() ) {
          const Node p = q.top();
          q.pop();

          if ( p.x - pts.back() > atol && xmax - p.x > atol )
            pts.emplace_back(p.x);

          auto& r = vectAt(runs, p.ri);
          const std::size_t i = p.i + 1;

          if ( i < r.npts ) {
            const double x = i == r.npts - 1
              ? r.x1
              : r.x0 + r.bw * static_cast<double>(i);
            q.push({ x, p.ri, i });
          }
        }

        pts.emplace_back(xmax);
        nc_assert( nc_is_grid(pts) );
        return pts;
      }

#else
      // Merges several equidistant grids into one sorted grid.
      // Points within 10% of the smallest input bin width are merged.
      // The output is not necessarily equidistant.
      // Its first and last points are the minimum and maximum input points.
      // Empty input produces an empty vector.
      // Each input grid must contain at least two finite points.
      VectD makeCommonGrid( Span<const EquidistantGrid> grids )
      {
#ifndef NDEBUG
        for ( auto& g : grids ) {
          nc_assert( std::isfinite(g.x0) );
          nc_assert( std::isfinite(g.binWidth) );
          nc_assert( g.binWidth > 0.0 );
          nc_assert( g.npts >= 2 );
          nc_assert( std::isfinite(g.x1()) );
        }
#endif

        constexpr bool do_pre_coalescing = true;
        constexpr std::size_t maxPending = 200000;
        constexpr std::size_t pruneGap = 10;

        if ( grids.empty() )
          return {};

        if ( grids.size() == 1 ) {
          auto& g = grids.front();
          return linspace(g.x0, g.x1(), g.npts);
        }

        struct Run {
          double x0, x1, bw;
          std::size_t npts;
        };

        double xmin = kInfinity;
        double xmax = -kInfinity;
        double minbw = kInfinity;
        std::size_t total = 0;

        for ( auto& g : grids ) {
          xmin = ncmin(xmin, g.x0);
          xmax = ncmax(xmax, g.x1());
          minbw = ncmin(minbw, g.binWidth);
          total += g.npts;
        }

        const double atol = 0.1 * minbw;
        nc_assert( atol > 1e-200 );

        std::vector<Run> runs;
        runs.reserve(grids.size());

        if ( do_pre_coalescing ) {
          for ( auto& g : grids ) {
            bool merged = false;

            for ( auto& r : runs ) {
              if ( r.bw != g.binWidth )
                continue;

              const double h = r.bw;
              const double ptol = ncmin(1e-5 * h, atol);
              const double phase = std::remainder(g.x0 - r.x0, h);

              if ( ncabs(phase) > ptol )
                continue;

              const double k = std::round((g.x0 - r.x0) / h);
              const double g0 = k;
              const double g1 = k +
                static_cast<double>(g.npts - 1);
              const double r1 = static_cast<double>(r.npts - 1);

              if ( g0 > r1 + 1.0 || g1 < -1.0 )
                continue;

              r.x0 = ncmin(r.x0, g.x0);
              r.x1 = ncmax(r.x1, g.x1());
              r.npts = static_cast<std::size_t>(
                                                std::round((r.x1 - r.x0) / h)) + 1;
              merged = true;
              break;
            }

            if ( !merged )
              runs.push_back({ g.x0, g.x1(), g.binWidth, g.npts });
          }
        }

        std::size_t runTotal = 0;

        if ( do_pre_coalescing ) {
          for ( auto& r : runs )
            runTotal += r.npts;
        }

        if ( !do_pre_coalescing || runTotal >= total ) {
          runs.clear();
          runs.reserve(grids.size());

          for ( auto& g : grids )
            runs.push_back({ g.x0, g.x1(), g.binWidth, g.npts });
        }

        VectD pts;
        pts.reserve(total <= maxPending ? total : maxPending);

        auto prune = [&]()
        {
          if ( pts.empty() )
            return;

          std::sort(pts.begin(), pts.end());

          std::size_t n = 1;
          vectAt(pts, 0) = xmin;

          for ( std::size_t i = 1; i + 1 < pts.size(); ++i ) {
            const double x = vectAt(pts, i);
            const double prev = vectAt(pts, n - 1);

            if ( x - prev > atol && xmax - x > atol )
              vectAt(pts, n++) = x;
          }

          if ( xmax > xmin )
            vectAt(pts, n++) = xmax;

          pts.resize(n);
        };

        std::size_t sincePrune = 0;

        for ( auto& r : runs ) {
          const std::size_t n = r.npts - 1;

          for ( std::size_t i = 0; i < n; ++i ) {
            const double x =
              r.x0 + r.bw * static_cast<double>(i);
            pts.emplace_back(x);
          }

          pts.emplace_back(r.x1);
          ++sincePrune;

          if ( pts.size() > maxPending &&
               (sincePrune >= pruneGap ||
                pts.size() > 2 * maxPending) ) {
            prune();
            sincePrune = 0;
          }
        }

        if ( sincePrune != 0 )
          prune();

        nc_assert( nc_is_grid(pts) );
        return pts;
      }


#endif


      bool detail_nearGrid(double a, double b, double tol)
      {
        return a != 0.0 && b != 0.0 &&
          ncabs(a - b) <= tol * ncmax(ncabs(a), ncabs(b));
      }

      double detail_sepGrid(double a, double b, double tol, double bw)
      {
        const double e = 1.0 +
          4.0 * std::numeric_limits<double>::epsilon();
        const double rel = tol * ncmax(ncabs(a), ncabs(b)) * e;
        const double nud = bw * 4.0 *
          std::numeric_limits<double>::epsilon();
        return ncmax(rel, nud);
      }

      void spaceOutGrid(VectD& g, double rtol);

      // Merge two finite, sorted, unique grids.
      // Endpoint values always take priority over averaging.
      VectD mergeGridsWithTol(const VectD& a, const VectD& b,
                              double rtol = 0.01 )
      {
        //fixme: no equal_rtol?!?
        nc_assert( rtol > 0.0 );
        nc_assert( rtol < 0.5 );
        nc_assert(nc_is_grid(a));
        nc_assert(nc_is_grid(b));
        nc_assert(!a.empty());
        nc_assert(!b.empty());
        nc_assert(std::isfinite(rtol) && rtol > 0.0);


        const double fac = 1.0 + rtol;
        const double lo = ncmin(a.front(), b.front());
        const double hi = ncmax(a.back(), b.back());

        VectD g;
        g.reserve(a.size() + b.size());

        auto add = [&g](double x)
        {
          if (g.empty() || x != g.back())
            g.push_back(x);
        };

        auto farEnough = [fac](double x, double y)
        {
          if (x == y || x == 0.0 || y == 0.0)
            return x != y;

          if ((x < 0.0) != (y < 0.0))
            return true;

          const double ax = ncabs(x);
          const double ay = ncabs(y);

          return ncmax(ax, ay) / ncmin(ax, ay) > fac;
        };

        std::size_t ia = 0;

        for (std::size_t ib = 0; ib < b.size(); ++ib) {
          const double x = vectAt(b, ib);

          while (ia < a.size() && vectAt(a, ia) < x) {
            add(vectAt(a, ia));
            ++ia;
          }

          bool keep = x == lo || x == hi;

          if (!keep && ia < a.size())
            keep = farEnough(x, vectAt(a, ia));

          if (!keep && ia > 0)
            keep = farEnough(x, vectAt(a, ia - 1));

          if (keep)
            add(x);
        }

        while (ia < a.size()) {
          add(vectAt(a, ia));
          ++ia;
        }

        return g;
      }





























      //   const double lo = ncmin(a.front(), b.front());
      //   const double hi = ncmax(a.back(), b.back());

      //   VectD g;
      //   g.reserve( a.size() + b.size() );

      //   auto itA = a.begin();
      //   auto itAE = a.end();
      //   auto itB = b.begin();
      //   auto itBE = b.end();
      //   //initial point:
      //   bool last_is_A;
      //   if ( *itA < *itB ) {
      //     last_is_A = true;
      //     g.push_back( *itA++ );
      //   } else {
      //     last_is_A = false;
      //     g.push_back( *itB++ );
      //   }

      //   auto acceptIfSpacedOut = [&g,rtol](double bb)
      //   {
      //     nc_assert(!g.empty());
      //     auto aa = g.back();
      //     bool ok = false;
      //     if ( bb>aa) {
      //       if ( aa < 0 ) {
      //         if ( aa < bb*(1.0+rtol) )
      //           ok = true;
      //       } else {
      //         if ( aa == 0.0 )
      //           ok = true;
      //         else if ( bb > (1.0+rtol)*aa )
      //           ok = true;
      //       }
      //     }
      //     if ( ok )
      //       g.push_back(bb);
      //     return ok;
      //   };

      //   while ( itA != itAE && itB !=itBE ) {
      //     if ( last_is_A ) {
      //       if ( *itA <= *itB ) {
      //         //two pts from A => accept without rtol check
      //         g.push_back( *itA++ );
      //         nc_assert( last_is_A == true );
      //       } else {
      //         if ( acceptIfSpacedOut( *itB++ ) )
      //           last_is_A = false;
      //       }
      //     } else {
      //       if ( *itB <= *itA ) {
      //         g.push_back( *itB++ );
      //         nc_assert( last_is_A == false );
      //       } else {
      //         if ( acceptIfSpacedOut( *itA++ ) )
      //           last_is_A = true;
      //       }
      //     }
      //   }
      //   //remainder:
      //   nc_assert( itA == itAE || itB == itBE );
      //   while ( itA != itAE ) {
      //     if ( !last_is_A ) {
      //       if ( acceptIfSpacedOut( *itA++ ) )
      //         last_is_A = true;
      //     } else {
      //       g.push_back( *itA++ );
      //     }
      //   }
      //   while ( itB != itBE ) {
      //     if ( last_is_A ) {
      //       if ( acceptIfSpacedOut( *itB++ ) )
      //         last_is_A = false;
      //     } else {
      //       g.push_back( *itB++ );
      //     }
      //   }
      //   g.back() = hi;
      //   return g;





      //   std::size_t i = 0;
      //   std::size_t j = 0;


      //   while (i < a.size() || j < b.size()) {
      //     double x;

      //     if (j == b.size() ||
      //         (i < a.size() && vectAt(a, i) <= vectAt(b, j)))
      //       x = vectAt(a, i++);
      //     else
      //       x = vectAt(b, j++);

      //     if (g.empty()) {
      //       g.push_back(x);
      //       continue;
      //     }

      //     // Exact equality must be handled separately because nearGrid()
      //     // deliberately does not regard zero as being within any tolerance.
      //     if (x == g.back())
      //       continue;

      //     if (detail_nearGrid(g.back(), x, equal_rtol)) {
      //       const double y = 0.5 * (g.back() + x);

      //       // Rounding can make the average equal to the preceding point.
      //       // Keeping x in that case preserves strict sorted uniqueness.
      //       if (g.size() == 1 ||
      //           y > vectAt(g, g.size() - 2))
      //         g.back() = y;
      //       else
      //         g.back() = x;
      //     } else {
      //       g.push_back(x);
      //     }
      //   }

      //   nc_assert(nc_is_grid(g));
      //   if (g.size() > 2)
      //     spaceOutGrid(g, rtol);

      //   // Preserve both domain endpoints.  If an endpoint was averaged with a
      //   // nearby non-endpoint, remove the non-endpoint instead.
      //   if (g.size() == 1) {
      //     if (lo != hi)
      //       g.push_back(hi);
      //     g.front() = lo;
      //     return g;
      //   }

      //   vectAt(g, 0) = lo;
      //   vectAt(g, g.size() - 1) = hi;

      //   if (g.size() > 2 &&
      //       detail_nearGrid(g.front(), vectAt(g, 1), equal_rtol))
      //     g.erase(g.begin() + 1);

      //   if (g.size() > 2 &&
      //       detail_nearGrid(vectAt(g, g.size() - 2), g.back(), equal_rtol))
      //     g.erase(g.end() - 2);

      //   return g;
      // }

#if 0
      //first slow attempt

      // Spread adjacent points that are too close according to rtol.
      // The first and last points are never changed.
      // Zero is never considered close to another point.
      // The nominal bin width is estimated from the current grid range.
      // If the domain cannot contain the requested separation, it is kept.
      void spaceOutGrid(VectD& g, double rtol)
      {
        nc_assert(nc_is_grid(g));
        nc_assert(std::isfinite(rtol) && rtol >= 0.0);

        if (g.size() < 3 || rtol == 0.0)
          return;

        const double bw = (g.back() - g.front()) /
          static_cast<double>(g.size() - 1);
        const std::size_t np = g.size() * 2 + 2;

        for (std::size_t pass = 0; pass < np; ++pass) {
          for (std::size_t i = 1; i + 1 < g.size(); ++i) {
            double& x = vectAt(g, i);
            const double l = vectAt(g, i - 1);
            const double r = vectAt(g, i + 1);

            if (detail_nearGrid(l, x, rtol)) {
              const double y = l + detail_sepGrid(l, x, rtol, bw);
              if (y < r)
                x = y;
            }

            if (detail_nearGrid(x, r, rtol)) {
              const double y = r - detail_sepGrid(x, r, rtol, bw);
              if (y > l)
                x = y;
            }
          }

          for (std::size_t i = g.size() - 1; i-- > 1;) {
            double& x = vectAt(g, i);
            const double l = vectAt(g, i - 1);
            const double r = vectAt(g, i + 1);

            if (detail_nearGrid(x, r, rtol)) {
              const double y = r - detail_sepGrid(x, r, rtol, bw);
              if (y > l)
                x = y;
            }

            if (detail_nearGrid(l, x, rtol)) {
              const double y = l + detail_sepGrid(l, x, rtol, bw);
              if (y < r)
                x = y;
            }
          }
        }
      }
#else
      // Move adjacent nonzero grid points apart when they are relatively close.
      // g must be a finite, strictly increasing grid with at least two points.
      // rtol is the relative separation tolerance and must be nonnegative.
      // The first and last points remain unchanged.
      // Zero is never considered close to another point.
      // Pairs that cannot be separated within the existing domain are retained.
      void spaceOutGrid(VectD& g, double rtol)
      {
#ifndef NDEBUG
        nc_assert(g.size() >= 2);
        nc_assert(std::isfinite(rtol) && rtol >= 0.0);
        for (const auto& x : g)
          nc_assert(std::isfinite(x));
        nc_assert(nc_is_grid(g));
#endif

        if (g.size() < 3 || rtol == 0.0)
          return;

        const std::size_t npts = g.size();
        const double bw = (g.back() - g.front()) /
          static_cast<double>(npts - 1);

        nc_assert(std::isfinite(bw) && bw > 0.0);

        for (std::size_t i = 1; i + 1 < npts; ++i) {
          double& x = vectAt(g, i);
          const double l = vectAt(g, i - 1);
          const double r = vectAt(g, i + 1);

          if (detail_nearGrid(l, x, rtol)) {
            const double y = l + detail_sepGrid(l, x, rtol, bw);
            if (y < r)
              x = y;
          }
        }

        for (std::size_t i = npts - 1; i-- > 1;) {
          double& x = vectAt(g, i);
          const double l = vectAt(g, i - 1);
          const double r = vectAt(g, i + 1);

          if (detail_nearGrid(x, r, rtol)) {
            const double y = r - detail_sepGrid(x, r, rtol, bw);
            if (y > l)
              x = y;
          }
        }
      }
#endif

      class OrderIntervalsByNearnessToZero final : private NoCopyMove {
        const double *m_x;
        std::size_t m_p, m_q, m_n;
      public:
        // Iterates over intervals [a,b] in a grid, ordered by |a| and with
        // positive a first in case of ties.
        explicit OrderIntervalsByNearnessToZero(Span<const double> x)
          : m_x(x.data()), m_p(0), m_q(0) {
          nc_assert( nc_is_grid(x) );
          m_n = x.size() - 1;
          const std::size_t z =
            std::lower_bound(x.begin(), x.end(), 0.0) - x.begin();
          m_p = z < m_n ? z : m_n;
          m_q = z ? z - 1 : m_n;
          if (m_q >= m_n)
            m_q = m_n;
        }
        bool hasMore() const { return m_p < m_n || m_q < m_n; }
        PairDD next()
        {
          nc_assert(hasMore());
          const bool takeP = ( m_p < m_n
                               && ( m_q >= m_n
                                    || ncabs(m_x[m_p]) <= ncabs(m_x[m_q]) ) );
          const std::size_t i = takeP ? m_p : m_q;
          nc_assert( i < m_n );
          PairDD res( m_x[i], m_x[i + 1] );
          if (takeP)
            ++m_p;
          else
            m_q = m_q ? m_q - 1 : m_n;
          return res;
        }
      };

      // Add points until the grid contains npts values.  New points are
      // inserted at centers of existing intervals closest to zero, if they are
      // wide enough. Might not succeed if not enough such intervals are
      // available (in that case, trying again with a lower rtol is advisable).
      void topOffGrid( VectD& g, std::size_t npts, double rtol = 0.1 )
      {
        nc_assert( npts > g.size() );
        g.reserve( npts );//Important to do this before initialising intervals
        {
          OrderIntervalsByNearnessToZero intervals(g);
          const double oneplusrtol = 1.0 + rtol;
          nc_assert_always( oneplusrtol > 1.0 );
          while ( g.size() < npts && intervals.hasMore() ) {
            double a,b;
            std::tie(a,b) = intervals.next();
            if ( a < 0.0 && b > 0.0 ) {
              //different signs => just ignore
              continue;
            }
            if ( a == 0.0 || b == 0.0
                 || ( ncmax(ncabs(a),ncabs(b))
                      > oneplusrtol*ncmin(ncabs(a),ncabs(b)) ) ) {
              g.push_back ( 0.5*(a+b) );
            }
          }
        }
        std::sort(g.begin(),g.end());
        nc_assert( nc_is_grid( g ) );
        if ( g.size() != npts ) {
          // nc_assert_always(false); //FIXME ok, but i would like to see if it ever happens
          nc_assert_always( rtol > 1e-4 );
          topOffGrid( g, npts, rtol*0.2 );
        }
        return;//fixme
#if 0
      // bool detail_nearGrid(double a, double b, double tol)
      // {
      //   return a != 0.0 && b != 0.0 &&
      //     ncabs(a - b) <= tol * ncmax(ncabs(a), ncabs(b));
      // }

        //FIXME: Easy one by TK for checking:
        nc_assert( nc_is_grid(g) );
        std::vector<std::pair<double,double>> gaps;//gapsize, newxval
        gaps.reserve(g.size());
        const double minrelgap = ncsquare(1.0+rtol) * 1.00001;
        nc_assert( minrelgap > 1.0 );
        for ( std::size_t i = 1; i < g.size(); ++i ) {
          double x0 = vectAt(g,i-1);
          double x1 = vectAt(g,i);
          double xsmall = ncmin(ncabs(x0),ncabs(x1));
          double xlarge = ncmax(ncabs(x0),ncabs(x1));
          const double relgap = ( xsmall ? xlarge/xsmall : kInfinity );
          if ( x0*x1 > 0.0//only consider x0, x1 of same sign, for simplicity
               && xlarge > minrelgap*xsmall ) {
            const double xnew = 0.5*(x0+x1);
#ifndef NDEBUG
            nc_assert( x0 < xnew );
            nc_assert( xnew < x1 );
            if ( x0 > 0.0 ) {
              nc_assert( x1 > (1.0+rtol)*xnew );
              nc_assert( xnew > (1.0+rtol)* x0 );
            } else {
              nc_assert( xnew < (1.0+rtol)* x1 );
              nc_assert( x0* 1 > (1.0+rtol)* xnew );
            }
#endif
            gaps.emplace_back( xlarge - xsmall, xnew );
            nc_assert( gaps.back().first > 0.0 );
          }
        }
        std::sort( gaps.begin(), gaps.end() );
        if ( g.size()+gaps.size() < npts )
          NCRYSTAL_THROW(BadInput,"Insufficient space in grid to"
                         " place additional points");
        auto it = gaps.begin();
        g.reserve(npts);
        while ( g.size() < npts ) {
          nc_assert( it < gaps.end() );
          g.push_back( (it++)->second );
          NCRYSTAL_MSG("Inserted "<<g.back());
        }
        std::sort(g.begin(),g.end());
        for ( auto& e : g )
          NCRYSTAL_MSG("Grid "<<e);
        nc_assert( nc_is_grid(g) );
        return;
#endif

        nc_assert(nc_is_grid(g));
        nc_assert(g.size() >= 2);
        nc_assert(npts >= 2);
        nc_assert(std::isfinite(rtol) && rtol >= 0.0);
        nc_assert(npts >= g.size());

        const double bw = (g.back() - g.front()) /
          static_cast<double>(npts - 1);

        while (g.size() < npts) {
          std::size_t k = g.size();
          double best = -1.0;

          for (std::size_t i = 0; i + 1 < g.size(); ++i) {
            const double d = vectAt(g, i + 1) - vectAt(g, i);

            if (d >= bw && d > best) {
              best = d;
              k = i;
            }
          }

          if (k == g.size()) {
            for (std::size_t i = 0; i + 1 < g.size(); ++i) {
              const double d = vectAt(g, i + 1) - vectAt(g, i);

              if (d > best) {
                best = d;
                k = i;
              }
            }
          }
          nc_assert(k < g.size());
          const double x = 0.5 *
            (vectAt(g, k) + vectAt(g, k + 1));
          nc_assert(!detail_nearGrid(x, vectAt(g, k), rtol));
          nc_assert(!detail_nearGrid(x, vectAt(g, k + 1), rtol));
          g.insert(g.begin() + static_cast<std::ptrdiff_t>(k + 1), x);
        }
        spaceOutGrid(g, rtol);
        nc_assert(g.size() == npts);
        nc_assert(nc_is_grid(g));
      }

      struct MergedGn {
        unsigned nmin, nmax;
        double binWidth;
        VectD yvals;
        double x0;
        double x1() const { return x0 + binWidth*(yvals.size()-1); }
      };


      // Adds the piecewise-linear function f evaluated at x_grid to target.
      // Values outside [f.x0, f.x1()] contribute zero.
      // f.yvals contains the values at equally spaced points beginning at f.x0.
      // x_grid must contain at least two finite, sorted, unique values.
      // target must have the same size as x_grid.
      // Input validation is performed in debug builds.
      void addEvaluationToGrid(const MergedGn& f, const VectD& x_grid, VectD& target)
      {
#ifndef NDEBUG
        nc_assert(x_grid.size() >= 2);
        nc_assert(target.size() == x_grid.size());
        nc_assert(f.yvals.size() >= 2);
        nc_assert(std::isfinite(f.x0));
        nc_assert(std::isfinite(f.binWidth));
        nc_assert(f.binWidth > 0.0);
        nc_assert(std::isfinite(f.x1()));

        for (const auto& v : f.yvals)
          nc_assert(std::isfinite(v));

        for (const auto& g : x_grid)
          nc_assert(std::isfinite(g));

        for (const auto& v : target)
          nc_assert(std::isfinite(v));

        nc_assert(nc_is_grid(x_grid));
#endif

        const std::size_t n = f.yvals.size();
        const double hi = f.x1();
        auto it = std::lower_bound(x_grid.begin(), x_grid.end(), f.x0);
        std::size_t j = static_cast<std::size_t>(it - x_grid.begin());

        for (; it != x_grid.end() && *it <= hi; ++it, ++j) {
          const double g = *it;

          if (g == hi) {
            vectAt(target, j) += f.yvals.back();
            continue;
          }

          const double q = (g - f.x0) / f.binWidth;
          const std::size_t k =
            q >= static_cast<double>(n - 1)
            ? n - 2
            : static_cast<std::size_t>(q);
          const double xk = f.x0 + f.binWidth * static_cast<double>(k);
          const double t = (g - xk) / f.binWidth;

          vectAt(target, j) +=
            vectAt(f.yvals, k) * (1.0 - t) +
            vectAt(f.yvals, k + 1) * t;
        }
      }

      double combinedGnFctWeight( VDOSGn::Order n )
      {
        // We want to assign higher weight to lower n values, and highest of all to
        // n=1. w(n)=1/n is unstable since sum of all weights does not converge, so
        // we go for ~1/n^2 instead. We also choose a form where the sum to n=inf is
        // normalised, and the parameter we tune is the weight of n=1. It is a bit
        // of trial and error, and likely not super crucial how this is setup. But
        // the following seemed to give sensible results.
#if 0
        return 1.0/(double(n.value()));
#endif

        constexpr double p1 = 0.3;//1st order phonon weight
        constexpr double a = 1.0/p1-2.0;
        constexpr double ap1 = a + 1.0;
        return ap1 / ((n.value()+a)*(n.value()+ap1));
      }



      // Calculates the trapezoidal integral of a non-negative piecewise-linear
      // function represented by x and y.
      // Removes as many points from the back as possible while keeping the
      // discarded integral at most frac times the original integral.
      // The vectors are modified in place; the original integral is returned.
      // frac defaults to 1e-6.
      void trimTail( VectD& x, VectD& y, double frac )
      {
#ifndef NDEBUG
        const auto npts = x.size();
        nc_assert(npts == y.size());
        nc_assert(npts >= 2);
        nc_assert(std::isfinite(frac) && frac >= 0.0);
        nc_assert( nc_is_grid(x) );
        for (auto& g : y) {
          nc_assert(std::isfinite(g));
          nc_assert(g >= 0.0);
        }
#endif
        double integral;
        {
          StableSumKahan sum;
          for (std::size_t i = 1; i < x.size(); ++i) {
            const double dx = vectAt(x, i) - vectAt(x, i - 1);
            sum.add( dx * (vectAt(y, i - 1) + vectAt(y, i)) );
          }
          integral = sum.sum()*0.5;
        }

        const double limit = frac * integral;
        double removed = 0.0;

        while (x.size() >= 2) {
          const auto i = x.size() - 2;
          const double dx = vectAt(x, i + 1) - vectAt(x, i);
          const double area = 0.5 * dx * (vectAt(y, i) + vectAt(y, i + 1));
          if (removed + area > limit)
            break;
          removed += area;
          x.pop_back();
          y.pop_back();
        }
      }
    }
  }
}

std::pair<NC::VectD,NC::VectD>
NC::VDOS::getCombinedGnFct( const VDOSGn& Gn )
{
  const auto nmax = Gn.maxOrder().value();
  const double invkT = 1.0 / Gn.kT();
#if 1
  std::vector<PWLFct> fs;
  VectD ws;
  fs.reserve(nmax);
  ws.reserve(nmax);
  for ( auto nm1 : ncrange(nmax) ) {
    const auto n = nm1+1;
    fs.emplace_back();
    fs.back().x0 = Gn.eRange(n).first*invkT;
    fs.back().binWidth = Gn.binWidth(n)*invkT;
    fs.back().f = Gn.getRawSpectrum(n);
    ws.push_back(combinedGnFctWeight( n ));
  };
      // struct PWLFct {
      //   double x0; //x{i=0}
      //   double binwidth;//distance between x{i} and x{i+1}
      //   Span<const double> f;//values of f at the x{i} points. The size of the
      //                        //span encodes the number of points.
      //   double xmax() const { return x0 + (f.size()-1)*binwidth; }
      // };

      // VectD evalPWLSum(Span<const PWLFct> fs,
      //                  Span<const double> ws,
      //                  Span<const double> grids)

  std::vector<EquidistantGrid> individual_grids;//fixme just use fs above!
  individual_grids.reserve(nmax);
  {
    unsigned n=0;
    for ( auto& f : fs ) {
      ++n;
      double bw = f.binWidth;
      std::size_t npts = f.f.size();
      //Speed up (a lot!) by already at this stage thinning out point candidates
      //at higher n values.
      while ( npts > 1000 && n > 2 ) {
        bw *= 2;
        npts /= 2;
      }
      while ( npts > 200 && n > 4 ) {
        bw *= 2;
        npts /= 2;
      }
      while ( npts > 100 && n > 10 ) {
        bw *= 2;
        npts /= 2;
      }
      while ( npts > 40 && n > 40 ) {
        bw *= 2;
        npts /= 2;
      }
      individual_grids.emplace_back();
      individual_grids.back().x0 = f.x0;
      individual_grids.back().binWidth = bw;
      individual_grids.back().npts = npts;
    }
  }

  VectD grid = makeCommonGrid( individual_grids );
  VectD vals = evalPWLSum(fs,grid,ws);
  std::pair<VectD,VectD> res;
  res.first = std::move(grid);
  res.second = std::move(vals);
  return res;
  (void)addEvaluationToGrid;//fixme
  (void)combinePWLFctsOfSameBinWidth;//fixme
#endif
#if 0





  // for ( auto nm1 : ncrange(Gn.maxOrder().value()) ) {
  //   NCRYSTAL_MSG(" G"<<(nm1+1)<<" binWidth: "<<fmt(Gn.binWidth(nm1+1)));
  // }

  //Note: Gn.binWidth() and Gn.eRange() return values in units of eV not
  //      kT (beta). We will do the work in these units, and only convert
  //      the final result to Beta.

  nc_assert_always( nmax >= 1 );

  //Group into ranges of identical binwidth (these can be merged without
  //interpolation and potentially in multiple threads):
  std::vector<MergedGn> mergedGns;
  mergedGns.reserve(16);
  mergedGns.emplace_back();
  mergedGns.back().nmin = 1;
  mergedGns.back().binWidth = Gn.binWidth(1);
  for ( unsigned n=2; n <= nmax; ++n ) {
    const double bw = Gn.binWidth(n);
    //bin width jumps in at least steps of 2, so checking if >1.5 is very
    //FP safe.
    // NCRYSTAL_MSG("TKTEST mergedGns bw check: "<<fmt(bw)<<" vs "<<fmt(mergedGns.back().binWidth));
//fixme was the other way:    nc_assert(floateq(bw,mergedGns.back().binWidth)
//fixme was the other way:              || mergedGns.back().binWidth > 1.99*bw);
    nc_assert(floateq(bw,mergedGns.back().binWidth)
              || bw > 1.99*mergedGns.back().binWidth);
    if ( bw > mergedGns.back().binWidth * 1.5 ) {
      mergedGns.back().nmax = n-1;
      mergedGns.emplace_back();
      mergedGns.back().nmin = n;
      mergedGns.back().binWidth = bw;
    }
  }
  mergedGns.back().nmax = nmax;
  {
    //Do the merging of all Gns at same binwidth multithreaded:
    FactoryJobs jobs;
    for ( auto ijob : ncrange(mergedGns.size()) ) {
      MergedGn* mgptr = &vectAt(mergedGns,ijob);
      jobs.queue([mgptr,&Gn]()
      {
        MergedGn& mg = *mgptr;
        const auto ngn = mg.nmax+1-mg.nmin;
        std::vector<Span<const double>> yarrays;
        VectD x0vals, weights;
        yarrays.reserve(ngn);
        x0vals.reserve(ngn);
        weights.reserve(ngn);
        for ( unsigned n = mg.nmin; n<=mg.nmax; ++n ) {
          yarrays.emplace_back(Gn.getRawSpectrum(n));
          x0vals.emplace_back(Gn.eRange(n).first);
          nc_assert( floateq( (Gn.eRange(n).second-Gn.eRange(n).first)/Gn.getRawSpectrum(n).size(), mg.binWidth) );
          nc_assert( floateq( mg.binWidth, Gn.binWidth(n) ) );
          weights.emplace_back(combinedGnFctWeight(n));
        }
        // NCRYSTAL_MSG("TKTEST Calling combinePWLFctsOfSameBinWidth for n="<<mg.nmin<<".."<<mg.nmax);

        std::tie(mg.yvals,mg.x0)
          = combinePWLFctsOfSameBinWidth( yarrays, x0vals,
                                          mg.binWidth, weights );
      });
    }
    jobs.waitAll();
  }
  // for ( auto& e : mergedGns )
  //   NCRYSTAL_MSG("TKTEST mergedGns nmin="<<e.nmin<<" nmax="<<e.nmax<<" bw="<<fmt(e.binWidth)<<" npts="<<e.yvals.size()
  //                <<" x0="<<e.x0<<" x1="<<e.x0+e.binWidth*(e.yvals.size()-1));

  VectD final_grid;
  //double min_binwidth = kInfinity;
  {
    std::vector<EquidistantGrid> grids;
    for ( auto& e : mergedGns ) {
      grids.emplace_back();
      grids.back().x0 = e.x0;
      grids.back().binWidth = e.binWidth;
      //min_binwidth = ncmin(min_binwidth,e.binWidth);
      grids.back().npts = e.yvals.size();
    }
    final_grid = makeCommonGrid( grids );
  }
  // NCRYSTAL_MSG("TKTEST Common grid has "<<final_grid.size()<<" pts from "<<fmt(final_grid.front())
  //              <<" to "<<fmt(final_grid.back()));



  VectD final_contrib;
  final_contrib.resize(final_grid.size(),0.0);
  for ( auto& f : mergedGns )
    addEvaluationToGrid(f, final_grid, final_contrib);


  //Almost done, but still need to convert to beta:

  //fixme: not here!
  // const double betaMax = gnexpn.upper_beta;
  // const double energyMax = betaMax * Gn.kT();

  for ( auto& e : final_grid )
    e *= invkT;

  std::pair<VectD,VectD> res;
  res.first = std::move(final_grid);
  res.second = std::move(final_contrib);
  return res;
#endif
}

std::pair<NC::VectD,NC::VectD>
NC::VDOS::determineAlphaBetaGridFromGn( const GnExpansion& gnexpn,
                                        unsigned nalpha, unsigned nbeta )//fixme switch order to nalpha,nbeta
{
  const unsigned requested_nalpha = nalpha;
  const unsigned requested_nbeta = nbeta;
  NCRYSTAL_MSG("TKTEST #beta: TOTAL = "<<requested_nbeta);
  NCRYSTAL_MSG("TKTEST #alpha: TOTAL = "<<requested_nalpha);


  //FIXME: before returning the final alpha/beta grids, go through and ensure
  //all pts are at least 1% apart, to avoid numerical issues in downstream algs.

  //For high quality performance in the low E limit, we dedicate some points
  //just for that purpose:
  VectD e0grid;
  {
    // NCRYSTAL_MSG("TKTEST orig nbeta="<<nbeta<<" nalpha="<<nalpha);
    unsigned npts_E0
      = std::max<unsigned>(12u,static_cast<unsigned>
                           (std::min<unsigned>(nalpha,nbeta)*0.2+0.5));
    VectD tmp;
    std::tie(e0grid,tmp) = setupE0ABGrid( gnexpn, npts_E0 );
    NCRYSTAL_MSG("TKTEST #beta: e0grid.size() = "<<e0grid.size());
    nc_assert_always(e0grid.size()>=5);
    nc_assert_always(nbeta > e0grid.size()+25);
    nc_assert_always(nalpha > e0grid.size()+20);
    nbeta -= static_cast<unsigned>(e0grid.size());
    nalpha -= static_cast<unsigned>(e0grid.size());
    nc_assert_always( nbeta > 30 && nbeta < 1000000 );
    nc_assert_always( nalpha > 20 && nalpha < 1000000 );
  }


  const double alpha2x = gnexpn.alpha2x;
  const double alphaMax = gnexpn.upper_alpha;
  const double betaMax = gnexpn.upper_beta;
  nc_assert_always(alpha2x > 0.0);
  nc_assert_always(alphaMax > 0.0);
  nc_assert_always(betaMax > 0.0);

  VectD bvals, gnprojvals;
  std::tie(bvals, gnprojvals) = getCombinedGnFct(gnexpn.Gn);//
  nc_assert_always( bvals.size() >= 5 );
  nc_assert_always( bvals.front() < 0.0 );
  nc_assert_always( bvals.back() > 0.0 );

  Span<const double> bvals_view(bvals);
  Span<const double> gnprojvals_view(gnprojvals);
  {
    //Discard pts outside [-betaMax,betaMax], but occasionally keep one point
    //extra, to avoid an edge-effects due to a sudden extrapolation towards 0 in
    //the edge region.
    std::size_t i = 0;
    while ( bvals_view[i] < -betaMax )
      ++i;
    if ( i > 0 && bvals_view[i] > -betaMax )
      --i;//keep one point going over the edge.

    bvals_view = bvals_view.subspan(i);
    gnprojvals_view = gnprojvals_view.subspan(i);
    if ( bvals_view.back() > betaMax ) {
      //Do the same for the upper limit, although this is expected to happen
      //only extremely rarely in usual operations.
      auto newsize = bvals_view.size();
      nc_assert_always( bvals_view.size() >= 5 );
      while ( newsize > 2
              && ncmin(bvals_view[newsize-1],
                       bvals_view[newsize-2]) > betaMax ) {
        --newsize;
      }
      bvals_view = bvals_view.subspan(0,newsize);
      gnprojvals_view = gnprojvals_view.subspan(0,newsize);
    }

    nc_assert_always( bvals_view.size() >= 5 );
    nc_assert_always( bvals_view.front() < 0.0 );
    nc_assert_always( bvals_view.back() > 0.0 );
    nc_assert_always( bvals_view.size() == gnprojvals_view.size() );
  }

  //Reduce number of points:
  std::tie(bvals, gnprojvals)
    = reducePtsInDistribution( bvals_view, gnprojvals_view, nbeta );
  NCRYSTAL_MSG("TKTEST #beta: npts_from_Gn_shape = "<<nbeta);
  nc_assert_always( bvals.size() >= 2 );
  nc_assert_always( bvals.size() <= nbeta );
  nc_assert_always( bvals.size() == gnprojvals.size() );
  nc_assert_always( bvals.front() < 0.0 );
  nc_assert_always( bvals.back() > 0.0 );

  //FIXME: WE SHOULD INCLUDE BETA=0 ALWAYS!!! WE SHOULD ALSO ALWAYS INCLUDE E0
  //PTS IN BOTH ALPHA AND BETA GRIDS.

  //Alpha is easier, since we know the behaviour of the formulas
  //f(x,n)=exp(-x)*x^n/n!. Investigations show that a powspace with p=2 gives
  //consistent relative errors over the entire range. This is the same as
  //linearly spacing in sqrt(x) values, and since the f(x,n) values approach
  //gaussians with sigma=sqrt(n) at high x values, this is perhaps not
  //surprising in hindsight, at least not at high x values.
  //
  //However, it might be a good idea to use p=2.5 or p=3.0 which reduces the
  //precision a bit at high x-values and improves it at low x values where we
  //might be more interested in a good behaviour. FIXME revisit!
  //
  //fixme: put alphagrid in separate function if we are not using the beta
  //       points!

  const double xmax = alpha2x*alphaMax;
  nc_assert_always(alpha2x>0.0&&std::isfinite(alpha2x));
  const double x2alpha = 1.0/alpha2x;
  nc_assert_always(x2alpha>0.0&&std::isfinite(x2alpha));

  //Always lead with 0.0 + a small geomspace, to get around abysmal log-interpolation
  //artifacts for G1 which is linear at low x. But then transition into the
  //powspace which works well for the rest of the range..



  constexpr double alphapow = 2.0;//fixme 2.5? 3.0?
  constexpr double alphapow_low = 4.0;//fixme tune
  unsigned nlow = static_cast<unsigned>(nalpha*0.08+4.5);
  unsigned nmid = static_cast<unsigned>(nalpha*0.12+10.5);
  const unsigned nzero = 1;//always begin with 0
  //Fixme: initially I simply had const nlow+nmid and:
#if 0
    nc_assert_always(nalpha > nlow+nmid+nzero+10);
#else
  while ( !(nalpha > nlow+nmid+nzero+5) ) {
    nlow -= 1;
    nmid -= 1;
  }
  nc_assert_always(nlow >= 3);
  nc_assert_always(nmid >= 8);
  nc_assert_always(nalpha > nlow+nmid+nzero+5);
#endif

  const auto npow = nalpha-(nlow+nmid+nzero);
  const double xpow0 = ncmin(2.5,0.2*xmax);
  const double xmid0 = ncmin(0.15,0.2*xpow0);
  const double xlow0 = ncmin( 0.005, 0.1*xmid0 );//0.005 since the fall-back to
                                                //linear interpolation over
                                                //[0.0,xlow0] is exactly the
                                                //~linear behaviour exhibited by
                                                //the n=1 curve at low x!

  auto xlow = powspace(xlow0,xmid0,nlow+1,alphapow_low);
  xlow.pop_back();
  //auto xmid = linspace(xmid0,xpow0,nmid+1);
  auto xmid = powspace(xmid0,xpow0,nmid+1,1.5);//FIXME linspace, pow??
  xmid.pop_back();
  auto xpow = powspace( xpow0, xmax, npow, alphapow );

  NCRYSTAL_MSG("TKTEST xlow0*x2alpha = "<<xlow0*x2alpha);
  NCRYSTAL_MSG("TKTEST xmid0*x2alpha = "<<xmid0*x2alpha);
  NCRYSTAL_MSG("TKTEST xpow0*x2alpha = "<<xpow0*x2alpha);

  VectD avals;
  avals.reserve( requested_nalpha );
  avals.push_back(0.0);
  for ( auto& x : xlow )
    avals.push_back( x*x2alpha );
  NCRYSTAL_MSG("TKTEST #alpha: n-zero = "<<1);
  NCRYSTAL_MSG("TKTEST #alpha: n-xlow = "<<xlow.size());
  NCRYSTAL_MSG("TKTEST #alpha: n-xmid = "<<xmid.size());
  NCRYSTAL_MSG("TKTEST #alpha: n-xpow = "<<xpow.size());

  for ( auto& x : xmid )
    avals.push_back( x*x2alpha );
  for ( auto& x : xpow )
    avals.push_back( x*x2alpha );
  nc_assert_always(floateq(avals.back(),alphaMax));
  avals.back() = alphaMax;//force exact
  nc_assert_always(nc_is_grid(avals));//fixme always
  nc_assert_always(avals.size()==nalpha);

  //Before returning we need to merge the e0grid points into both alpha and beta
  //grids, ensure we avoid pts too close, and to ensure the exact final
  //count. For alpha, it is possible for extreme values of alpha2x (e.g. in
  //low-temp uranium) to get a gap between the last e0 pt and the first xlow
  //point. If this happens, we try to fill it out.
  nc_assert( !e0grid.empty() );

  //  NCRYSTAL_MSG("TKTEST e0grid.back() = "<<e0grid.back()<<" avals.front() = "<<avals.front());
  if ( e0grid.back() < xlow.front()*x2alpha ) {
    auto nn = static_cast<std::size_t>(e0grid.size() * 0.1+0.5);
    e0grid.resize( e0grid.size()-nn );
    auto tmp = linspace( e0grid.back(), xlow.front()*x2alpha, nn + 2 );
    for ( auto e : Span<const double>(tmp).subspan(1,nn) )
      e0grid.push_back( e );
  }
  {
    //FIXME TOPOFFGRID should go back and look for the last points removed in initial pt reduction
    //    const double equal_rtol = 1e-4;//collapse pts equal
    const double rtol = 0.01;//avoid pts too close to each other
    const std::size_t avals_size_withoute0 = avals.size();
    avals = mergeGridsWithTol( avals, e0grid, rtol );
    NCRYSTAL_MSG("TKTEST #alpha: e0pts actually merged = "<<avals.size()-avals_size_withoute0);
    // NCRYSTAL_MSG("TKTEST alpha ntopoff:="<<(long)requested_nalpha - (long)avals.size());
    NCRYSTAL_MSG("TKTEST #alpha: top off pts = "<<requested_nalpha-avals.size());
    if ( avals.size() < requested_nalpha )
      topOffGrid(avals, requested_nalpha, rtol );
    // NCRYSTAL_MSG("TKTEST beta ntopoff pre-merge:="<<(long)requested_nbeta - (long)bvals.size());
    const std::size_t bvals_size_withoute0 = bvals.size();
    bvals = mergeGridsWithTol( bvals, e0grid, rtol );

    NCRYSTAL_MSG("TKTEST #beta: e0pts actually merged = "<<bvals.size()-bvals_size_withoute0);
    nc_assert_always( bvals.size() >= bvals_size_withoute0 );

    // NCRYSTAL_MSG("TKTEST beta ntopoff:="<<(long)requested_nbeta - (long)bvals.size());
    NCRYSTAL_MSG("TKTEST #beta: top off pts = "<<requested_nbeta-bvals.size());

  {
    VectD tmp( bvals );
    nc_assert( nc_is_grid( tmp ) );
    const double kT = gnexpn.Gn.kT();
    for ( auto& e : tmp )
      e *= kT;
    nc_assert( nc_is_grid( tmp ) );
  }
    if ( bvals.size() < requested_nbeta )
      topOffGrid(bvals, requested_nbeta, rtol );

  {
    VectD tmp( bvals );
    nc_assert( nc_is_grid( tmp ) );
    const double kT = gnexpn.Gn.kT();
    for ( auto& e : tmp )
      e *= kT;
    nc_assert( nc_is_grid( tmp ) );
  }

  }

  // NCRYSTAL_MSG("TKTEST final nbeta="<<bvals.size()<<" nalpha="<<avals.size());

  nc_assert_always( avals.size() == requested_nalpha );
  nc_assert_always( bvals.size() == requested_nbeta );


  return {avals, bvals};
}

std::pair<unsigned,unsigned> NC::VDOS::gridDimFromLux( unsigned vdoslux )
{
  const unsigned override_nbins = ncgetenv_int("HACK_NBINS");
  if ( override_nbins )
    return { override_nbins, override_nbins };

  nc_assert_always( vdoslux <= 5 );

  switch ( vdoslux ) {
  case 0: return {  119, 238 };//fixme: should we lower and just accept that vdoslux=0 is crap, and only useful for quick unit tests?
  case 1: return { 178, 356 };
  case 2: return { 267, 533 };//was 200x400
  default:
    //  case 3: return { 400, 800 };//was 400x800
  case 3: return { 280, 650 };//was 400x800
  case 4: return { 600, 1200 };
    //FIXME: case 5: return { 900, 1800 };
  case 5: return { 2000, 4000 };//fixme: just for now, for a solid ref
    //FIXME: case 6: return { 2000, 4000 };
  }

#if 0
  const unsigned nbeta = 100*(1<<vdoslux);//100 * 2^(vdoslux).
  const unsigned nalpha = nbeta/2;
#else
  const unsigned scale = (1<<vdoslux);//2^(vdoslux).
  const unsigned nbeta = 125*scale;
  const unsigned nalpha = 40*scale;
  // const unsigned nbeta = 200*scale;
  // const unsigned nalpha = 25*scale;
  // const unsigned nbeta = 151*scale;
  // const unsigned nalpha = 33*scale;
#endif
  return { nalpha, nbeta };
}

std::pair<NC::VectD,NC::VectD>
NC::VDOS::setupE0ABGrid( const GnExpansion& gnexpn, unsigned npts )
{
  nc_assert_always(npts >= 5);
  const VDOSGn& Gn = gnexpn.Gn;
  const double alpha2x = gnexpn.alpha2x;
  nc_assert_always(alpha2x>0.0);
  const double x2alpha = 1.0 / alpha2x;
  nc_assert_always(x2alpha>0.0);
  const auto nmax = Gn.maxOrder().value();
  const double kT = Gn.kT();
  const double invkT = 1.0/kT;
  std::vector<PWLFct> fcts;
  fcts.reserve(16);

  double minus_log_nfactorial = 0.0;// accumulate -ln(n!)
  for ( auto nm1 : ncrange(nmax) ) {
    const auto n = nm1+1;
    minus_log_nfactorial -= std::log(static_cast<double>(n));
    if ( n>1 ) {
      //check if we can break already
      const double betamax = Gn.eRange(n).second*invkT;
      if ( betamax <= 0.0 )
        break;
      constexpr double relcontriblvl = 1e-8;
      constexpr double accuracy = 1e-11;
      auto xrange = rangeXNexpMX(n, relcontriblvl, accuracy );
      if ( betamax < xrange.first * x2alpha )
        break;
    }

    //Collect contributions on the line alpha=beta, but modified with a factor
    //of sqrt(beta) to account for the relative width of the phase space at
    //different beta values.

    //shave off negative beta:
    fcts.emplace_back();
    auto& f = fcts.back();
    f.x0 = Gn.eRange(n).first * invkT;
    f.binWidth = Gn.binWidth(n)*invkT;
    nc_assert_always(f.binWidth>0.0);
    f.f = Gn.getRawSpectrum(n);
    //Discard non-positive values:
    f = pwlNarrowToPos(f);
    nc_assert_always( f.x0 > 0.0 && f.f.size() >= 2 );
    auto f_fmut = f.f_mutable();

    //Now add alpha and phase space factors:
    for (std::size_t i = 0; i < f.f.size(); ++i) {
      //note alpha=beta on the E->0 phasespace, so x=beta*alpha2x
      const double beta = f.x0 + i * f.binWidth;
      //relative phasespace width is proportional to sqrt(b) as E->0
      double factor = std::sqrt(beta);
      //Add also alpha factor: exp(-x)*x^n/n!:
      const double x = alpha2x*beta;//alpha = beta in the E->0 limit
      factor *= std::exp(-x + static_cast<double>(n)*std::log(x) + minus_log_nfactorial);
      f_fmut[i] *= factor;
    }
  }

  //Find common grid for all of these functions:
  VectD grid;
  {
    std::vector<EquidistantGrid> allgrids;
    allgrids.reserve(fcts.size());
    for ( auto& f : fcts ) {
      allgrids.emplace_back();
      allgrids.back().x0 = f.x0;
      allgrids.back().binWidth = f.binWidth;
      nc_assert_always(allgrids.back().binWidth>0.0);
      allgrids.back().npts = f.f.size();
    }
    grid = makeCommonGrid( allgrids );
    while ( grid.back() >= gnexpn.upper_beta )
      grid.pop_back();
    nc_assert_always( grid.size() >= 10 );
    //ensure we have 0.0 in this:
    if ( grid.front() < 1e-3 * fcts.front().binWidth ) {
      grid.front() = 0.0;
    } else {
      //A bit expensive:
      grid.insert(grid.begin(), 0.0);
    }
  }
  nc_assert_always(grid.size() > 5 && grid.front() == 0.0 );
  nc_assert( nc_is_grid(grid) );

  //Add up all contributions on the common grid (here with no ad-hoc weights,
  //since we actually incorporated the real weight factors into the fct values
  //just above):
  auto contrib = evalPWLSum(fcts,grid);

  //Do not keep extreme tail points with negligible impact:
  trimTail( grid, contrib, 1e-9 );

  nc_assert_always( grid.size() == contrib.size() );
  if ( npts < grid.size() )
    std::tie(grid, contrib) = reducePtsInDistribution( grid, contrib, npts );
  nc_assert_always( grid.size() <= npts );
  nc_assert_always( grid.size() == contrib.size() );
  nc_assert_always( grid.front() >= 0.0 );
  nc_assert_always( grid.back() <= gnexpn.upper_beta );

  std::pair<VectD,VectD> res;
  res.first = std::move(grid);
  res.second = std::move(contrib);
  return res;
}
