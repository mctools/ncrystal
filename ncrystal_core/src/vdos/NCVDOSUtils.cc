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

#include "NCrystal/internal/vdos/NCVDOSUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/phys_utils/NCKinUtils.hh"

namespace NC=NCrystal;

NC::PairDD NC::VDOS::rangeXNexpMX(unsigned n, double eps, double accuracy ) {
  //FIXME: cache high-res results for lowest n=200 orders for a limited set of
  //eps values? That would remove another source of irreproducibilities. Or, we
  //could always use a high accuracy, and build up a cache of already returned
  //values.

  nc_assert(eps>0.0&&eps<1.0&&eps>1e-200&&n>0&&accuracy>0&&accuracy<=1e-2);

  //The function f(x) = x^n*exp(-x) peaks at x=n and falls off on both
  //sides. Returns the two solutions to f(x)= f(n)*eps, describing the central
  //range around x=n where the function is higher than eps times the peak
  //value.

  //Must solve for x:
  //  x^n*exp(-x) = eps*[n^n*exp(-n)]
  //Raise to 1/n power and get:
  //  x*exp(-x/n) = eps^(1/n)*n*exp(-1)
  //<=> (x/n)*exp(-x/n) =  (1/e)*eps^(1/n) = k
  //
  //Which can be solved numerically for x/n:

  const double n_dbl = static_cast<double>(n);
  const double k = kInvE * std::pow(eps,1.0/n_dbl);
  auto f = [k](double y) { return y*std::exp(-y)-k; };
  return { n_dbl*findRoot2( f, 0.0,   1.0, accuracy ),
           n_dbl*findRoot2( f, 1.0, 700.0, accuracy ) };
}

NC::Rectangle NC::VDOS::findABExtentWithinKB( const Rectangle& r,
                                              double E_div_kT )
{
  nc_assert(std::isfinite(r.x0()));
  nc_assert(std::isfinite(r.x1()));
  nc_assert(std::isfinite(r.y0()));
  nc_assert(std::isfinite(r.y1()));
  nc_assert(std::isfinite(E_div_kT));
  const auto& alphaRange = r.xRange();
  const auto& betaRange = r.yRange();
  nc_assert(alphaRange.first >= 0.0);
  nc_assert(alphaRange.second > alphaRange.first);
  nc_assert(betaRange.second > betaRange.first);
  nc_assert(E_div_kT > 0.0);

  const double a0 = alphaRange.first;
  const double a1 = alphaRange.second;
  const double b0 = betaRange.first;
  const double b1 = betaRange.second;
  const double e = E_div_kT;

  double al = a0;
  double ah = a1;

  // beta^-(alpha) < b1.
  const auto lim1 = getAlphaLimits(e, b1);
  if (!(lim1.first < lim1.second))
    return {};

  // For a fixed alpha, the phase-space beta interval must overlap
  // the rectangle's beta interval.  The condition involving the
  // rectangle's upper edge b1 is:
  //
  //     beta^-(alpha) < b1
  //
  // getAlphaLimits(e, b1) gives the two alpha values where the
  // horizontal line beta=b1 crosses the phase-space boundary.
  //
  // If b1 < 0, beta^-(0)=0 is above b1, so alpha must be at
  // least lim1.first.  If b1 >= 0, beta^-(0)<=b1, so there
  // is no lower alpha restriction from b1.  In both cases,
  // alpha must not exceed lim1.second.
  if (b1 < 0.0)
    al = ncmax(al, lim1.first);
  ah = ncmin(ah, lim1.second);

  // b0 < beta^+(alpha).
  if (b0 >= 0.0) {
    const auto lim0 = getAlphaLimits(e, b0);
    al = ncmax(al, lim0.first);
  }

  if (!(al < ah))
    return {};

  // beta^-(alpha) has its minimum at alpha = e.
  const double amin = ncmin(ncmax(e, al), ah);
  const double betaLo = ncmax(b0, getBetaMinus(e, amin));
  const double betaHi = ncmin(b1, getBetaPlus(e, ah));
  if ( !(betaLo<betaHi) )
    return {};//unlikely, except due to numerical imprecision.
  return { al, ah, betaLo, betaHi };
}

namespace NCRYSTAL_NAMESPACE {
  namespace VDOS {
    namespace {

      static constexpr auto sogFailMsg = "Unable to space out grid";

      double positiveLowerBound(double b, double r)
      {
        double x = b / r;
        while (r * x > b) {
          const double next = std::nextafter(x, 0.0);
          if (next == x)
            NCRYSTAL_THROW(BadInput,sogFailMsg);
          x = next;
        }
        if (!std::isfinite(x) || x <= 0.0)
          NCRYSTAL_THROW(BadInput,sogFailMsg);
        return x;
      }

      double positiveUpperBound(double a, double r)
      {
        const double x = r * a;
        if (!std::isfinite(x) || x <= 0.0)
          NCRYSTAL_THROW(BadInput,sogFailMsg);
        return x;
      }

      template<bool IsPos>
      void spaceOutPass(Span<double> g, bool toRight, double rtol)
      {
#ifndef NDEBUG
        nc_assert(std::isfinite(rtol));
        nc_assert(rtol >= 0.0 && rtol < 1.0);
        nc_assert(IsPos ? g.front()>0.0 : g.back() < 0.0 );
        for (std::size_t i = 0; i < g.size(); ++i)
          nc_assert(std::isfinite(g[i]) &&
                    (IsPos ? g[i] > 0.0 : g[i] < 0.0));
#endif

        if (g.size() < 2)
          return;

        const double r = 1.0 + rtol;

        if (r == 1.0)
          return;

        if (toRight) {
          for (std::size_t i = 1; i < g.size(); ++i) {
            const double a = g[i - 1];
            const double b = g[i];

            if (IsPos) {
              if (b < r * a)
                g[i] = positiveUpperBound(a, r);
            } else {
              if (a > r * b) {
                // Reflect the negative problem around zero:
                //
                //   a <= r*b
                //   -a >= r*(-b)
                //
                // The new negative b is therefore the reflection of the
                // largest positive x satisfying r*x <= -a.
                g[i] = -positiveLowerBound(-a, r);
              }
            }
          }
        } else {
          for (std::size_t i = g.size() - 1; i != 0; --i) {
            const double a = g[i - 1];
            const double b = g[i];

            if (IsPos) {
              if (b < r * a) {
                // Choose x conservatively so that r*x <= b even after
                // the multiplication is rounded.
                g[i - 1] = positiveLowerBound(b, r);
              }
            } else {
              if (a > r * b) {
                // Reflect the positive upper-bound operation around zero.
                const double x = positiveUpperBound(-b, r);
                g[i - 1] = -x;
              }
            }
          }
        }
      }
    }
  }
}

void NC::VDOS::spaceOutGrid(Span<double> g, double rtol)
{
  nc_assert(nc_is_grid(g));//checks: sorted, unique, finite, size>=2
  nc_assert(std::isfinite(rtol));
  nc_assert(rtol >= 0.0 && rtol < 1.0);

  if (rtol == 0.0)
    return;

  double* const itB = g.data();
  double* const itE = itB + g.size();
  double* const itZero = std::lower_bound(itB, itE, 0.0);

  double* const itNegE = itZero;
  double* itPosB = itZero;

  if (itZero != itE && *itZero == 0.0) {
    *itZero = 0.0;// -0.0 -> 0.0
    ++itPosB;
  }

  Span<double> neg(itB, itNegE), pos(itPosB, itE);

  if (neg.size() >= 2) {
    const double origFirst = neg.front();
    spaceOutPass<false>(neg, false, rtol);
    if (neg.front() != origFirst) {
      neg.front() = origFirst;
      const double origLast = neg.back();
      spaceOutPass<false>(neg, true, rtol);
      if (neg.back() != origLast)
        NCRYSTAL_THROW(BadInput,sogFailMsg);
    }
  }
  if (pos.size() >= 2) {
    const double origLast = pos.back();
    spaceOutPass<true>(pos, true, rtol);
    if (pos.back() != origLast) {
      pos.back() = origLast;
      const double origFirst = pos.front();
      spaceOutPass<true>(pos, false, rtol);
      if (pos.front() != origFirst)
        NCRYSTAL_THROW(BadInput,sogFailMsg);
    }
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace VDOS {
    namespace {
      class OrderIntervalsByNearnessToZero final : private NoCopyMove {
        const double *m_x;
        std::size_t m_p, m_q, m_n;
        static double nearestEndpointAbs(const double *x, std::size_t i)
        {
          return ncmin(ncabs(x[i]), ncabs(x[i + 1]));
        }
      public:
        // Iterates over intervals [a,b] in a grid, ordered by the magnitude of
        // their nearest endpoint to zero, with the positive-side interval first
        // in case of ties.
        explicit OrderIntervalsByNearnessToZero(Span<const double> x)
          : m_x(x.data()), m_p(0), m_q(0)
        {
          nc_assert( nc_is_grid(x) );
          m_n = x.size() - 1;
          const std::size_t z =
            std::lower_bound(x.begin(), x.end(), 0.0) - x.begin();
          m_p = z < m_n ? z : m_n;
          m_q = z ? z - 1 : m_n;
          if (z > m_n) {
            //negative only grid
            m_q = m_n - 1;
          }
        }
        bool hasMore() const { return m_p < m_n || m_q < m_n; }
        PairDD next()
        {
          nc_assert(hasMore());
          const bool takeP = ( m_p < m_n
                               && ( m_q >= m_n
                                    || nearestEndpointAbs(m_x, m_p)
                                    <= nearestEndpointAbs(m_x, m_q) ) );
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
    }
  }
}

void NC::VDOS::topOffGrid( VectD& g, std::size_t npts, double rtol )
{
  nc_assert( nc_is_grid(g) );
  if ( g.size() >= npts )
    return;
  g.reserve( npts );//Important to do this before initialising intervals
  {
    OrderIntervalsByNearnessToZero intervals(g);
    const double oneplusrtol = 1.0 + rtol;
    nc_assert( oneplusrtol > 1.0 );
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
        const double mid = a + 0.5*(b-a);//stable since sign(a)==sign(b)
        if (a < mid && mid < b)
          g.push_back( mid );
      }
    }
  }
  std::sort(g.begin(),g.end());
  nc_assert(g.size() <= npts);
  nc_assert(nc_is_grid(g));
}

NC::VectD NC::VDOS::makeCommonGrid(Span<const EquidistantGrid> grids )
{
  struct Entry {
    double x;
    std::size_t grid;
    std::size_t point;
  };

  struct EntryGreater {
    bool operator()(const Entry& a, const Entry& b) const
    {
      return a.x > b.x;
    }
  };

#ifndef NDEBUG
  for ( auto& g : grids ) {
    nc_assert( std::isfinite(g.x0) );
    nc_assert( std::isfinite(g.binWidth) );
    nc_assert( g.binWidth > 0.0 );
    nc_assert( g.npts >= 2 );
    nc_assert( std::isfinite(g.x1()) );
    nc_assert( g.x0 + g.binWidth > g.x0 );
  }
#endif

  if ( grids.empty() )
    return {};

  if ( grids.size() == 1 ) {
    auto& g = grids.front();
    return linspace(g.x0, g.x1(), g.npts);
  }

  //NB: We might be able to deal with grids.size()==2 even more efficiently, if
  //it ever became very useful.

  double xmin(kInfinity), xmax(-kInfinity), minbw(kInfinity);
  for ( auto& g : grids ) {
    xmin = ncmin(xmin, g.x0);
    xmax = ncmax(xmax, g.x1());
    minbw = ncmin(minbw, g.binWidth);
  }
  const double atol = 0.1 * minbw;
  nc_assert( atol > 0.0 );
  std::vector<Entry> heap;
  heap.reserve(grids.size());
  for ( std::size_t j = 0; j < grids.size(); ++j ) {
    nc_assert( grids[j].npts >= 2);
    heap.push_back({grids[j].x0, j, 0});
  }
  EntryGreater greater;
  std::make_heap(heap.begin(), heap.end(), greater);
  auto siftDown = [&heap, &greater]()
  {
    const std::size_t n = heap.size();
    const Entry item = heap.front();
    std::size_t p = 0;
    while ( true ) {
      const std::size_t l = 2 * p + 1;
      if ( l >= n )
        break;
      std::size_t c = l;
      const std::size_t r = l + 1;
      if ( r < n && greater(vectAt(heap, l), vectAt(heap, r)) )
        c = r;
      if ( !greater(item, vectAt(heap, c)) )
        break;
      vectAt(heap, p) = vectAt(heap, c);
      p = c;
    }
    vectAt(heap, p) = item;
  };

  VectD pts;
  pts.reserve(12000);
  pts.emplace_back(xmin);

  while ( !heap.empty() ) {
    Entry item = heap.front();
#ifndef NDEBUG
    const auto oldPoint = item.point;
#endif
    const double prev = pts.back();
    nc_assert( item.x >= prev );
    nc_assert( xmax >= item.x );
    const auto& g = grids[item.grid];
    if ( item.x - prev > atol && xmax - item.x > atol )
      pts.emplace_back(item.x);
    if ( item.point == g.npts - 1 ) {
      std::pop_heap(heap.begin(), heap.end(), greater);
      heap.pop_back();
      continue;
    }
    ++item.point;
    nc_assert(item.point > oldPoint);
    nc_assert(item.point < g.npts);
    item.x = ( item.point == g.npts - 1
               ? g.x1()
               : g.x0 + g.binWidth * static_cast<double>(item.point) );
    heap.front() = item;
    siftDown();
    nc_assert( heap.empty()
               || heap.front().point < grids[heap.front().grid].npts );
  }
  if ( xmax > xmin )
    pts.emplace_back(xmax);
  nc_assert( nc_is_grid(pts) );
  return pts;
}

NC::VDOS::EquidistantGrid
NC::VDOS::coverEquidistantGrids( const EquidistantGrid& g1,
                                 const EquidistantGrid& g2 )
{
#ifndef NDEBUG
  auto valid = [](const EquidistantGrid& g) {
    nc_assert(g.npts > 0);
    nc_assert(g.binWidth > 0.0);
    nc_assert(std::isfinite(g.x0));
    nc_assert(std::isfinite(g.binWidth));
    nc_assert(std::isfinite(g.x1()));
  };
  valid(g1);
  valid(g2);
  nc_assert(g1.binWidth == g2.binWidth);
  nc_assert(intervalsOverlap(g1.x0, g1.x1(), g2.x0, g2.x1()));
#endif

  const double x0 = ncmin(g1.x0, g2.x0);

  auto lastIndex = [x0](const EquidistantGrid& g) {
    const double q = (g.x0 - x0) / g.binWidth;
    nc_assert(std::isfinite(q));
    nc_assert(q >= 0.0);
    const double first = std::ceil(q);
    nc_assert(first
              < static_cast<double>(std::numeric_limits<std::size_t>::max()));
    const std::size_t i0 = static_cast<std::size_t>(first);
    nc_assert(g.npts - 1 <= std::numeric_limits<std::size_t>::max() - i0);
    return i0 + g.npts - 1;
  };

  const std::size_t last = std::max<std::size_t>(lastIndex(g1), lastIndex(g2));

  //Given a double with...:
  static_assert(std::numeric_limits<double>::radix == 2, "");
  static_assert(std::numeric_limits<double>::digits == 53, "");
  static_assert(std::numeric_limits<std::size_t>::digits >= 54, "");
  //...we can represent all integers up to 2^53 exactly:
  constexpr std::size_t szdblmax = 9007199254740992ULL;
  (void)szdblmax;
  nc_assert( last+1 <= szdblmax );
  return EquidistantGrid{x0, g1.binWidth, last + 1};
}
