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
#include "NCrystal/internal/utils/NCIter.hh"
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
               : g.xAt(item.point) );
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
    //Tolerance, since (for our VDOS-Gn usage) g.x0-x0 is often an exact integer
    //number of bins (grids on the same lattice), and we should not get a
    //different result depending on which side of the integer rounding errors
    //put us:
    const double first = std::ceil(q - 1e-6);
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

NC::VDOS::PWLFct NC::VDOS::pwlNarrowToPos( const PWLFct& p, double tol )
{
  if ( p.f.size() < 2 ) {
    nc_assert( p.f.size() == 0 );
    return {};
  }
  nc_assert( p.binWidth > 0.0 );
  const double t = tol * p.binWidth;//threshold
  std::size_t i = ( p.x0 < t ? static_cast<std::size_t>
                    (std::ceil((t - p.x0) / p.binWidth)) : 0u );
  i = std::min<std::size_t>(i, p.f.size());

  // Correct possible floating-point rounding in the index calculation, and
  // ensure that results are consistent with what p.xAt(..) provides.
  while ( i > 0 && p.xAt(i - 1) >= t )
    --i;
  while ( i < p.f.size() && p.xAt(i) < t )
    ++i;

  if ( i >= p.f.size() || p.f.size() - i < 2 ) {
    // There are not two grid points at or above t:
    return {};
  }

  PWLFct res;
  res.x0 = p.xAt(i);
  res.binWidth = p.binWidth;
  res.dataHolder.assign(p.f.begin() + i, p.f.end());
  res.f = res.dataHolder;
  nc_assert_always( res.x0 >= t );
  nc_assert_always( res.f.size() >= 2 );
  nc_assert( res.dataHolder.size() == res.f.size() );
  return res;
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    NCRYSTAL_FMADISPATCH_ATTR
    void pwlSumAccumulateSegment( double* out, const double* gridPtr,
                                  std::size_t g, std::size_t end,
                                  double y0, double slope, double xLeft,
                                  double ylo, double yhi, double weight )
    {
      for ( std::size_t k = g; k < end; ++k ) {
        const double v = ncclamp( std::fma(slope,gridPtr[k]-xLeft,y0),
                                  ylo, yhi );
        out[k] = std::fma( weight, v, out[k] );
      }
    }
  }
}

NC::VectD NC::VDOS::evalPWLSum( Span<const PWLFct> fs,
                                Span<const double> grid,
                                Span<const double> ws )
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
    const double xmax = p.x1();

    // Grid points which are meant to be at the endpoints of the function (e.g.
    // because they are nodes of another function with the same endpoint) can be
    // a few ulps outside, due to rounding errors in the calculation of one or
    // the other. Since the function value at an endpoint might be very
    // different from the (zero) value just outside it, we treat points within a
    // tiny tolerance of the endpoints as being at the endpoints.
    const double xtol = 1e-9 * p.binWidth;

    std::size_t g = static_cast<std::size_t>
      (std::lower_bound(grid.begin(), grid.end(), p.x0 - xtol) - grid.begin());

    if (g == grid.size() || gridPtr[g] > xmax + xtol)
      continue;

    // For consistency, the bin edges must be calculated exactly like the node
    // positions (p.xAt(i)). Otherwise a grid point at a node can end up in the
    // neighbouring bin and get a value affected by rounding errors instead of
    // exactly the value at the node.
    double xLeft = p.x0;
    for (std::size_t i = 0; i < lastBin && g < grid.size(); ++i) {
      const double xRight = p.xAt(i + 1);
      const double y0 = vectAt(f, i);
      const double y1 = vectAt(f, i + 1);
      const double slope = (y1 - y0) * invbw;
      const double ylo = ncmin(y0, y1);
      const double yhi = ncmax(y0, y1);
      std::size_t end = g;
      while (end < grid.size() && gridPtr[end] < xRight)
        ++end;

      // The result is exactly y0 at x=xLeft, and the clamping ensures that the
      // result never falls outside [y0,y1] (or below 0 for non-negative data)
      // due to rounding errors.
      pwlSumAccumulateSegment( outPtr, gridPtr, g, end,
                               y0, slope, xLeft, ylo, yhi, weight );
      g = end;
      xLeft = xRight;
    }

    // Handle the final sample at x == xmax.
    if (g < grid.size() && gridPtr[g] <= xmax + xtol) {
      const double y = vectAt(f, lastBin);
      std::size_t end = g;
      while (end < grid.size() && gridPtr[end] <= xmax + xtol)
        ++end;
      for (std::size_t k = g; k < end; ++k)
        outPtr[k] = std::fma( weight, y, outPtr[k] );
    }
  }

  return out;
}

void NC::VDOS::trimTailByIntegral( VectD& x, VectD& y, double frac )
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

  while (x.size() > 2) {
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

void NC::VDOS::trimEquidistantGridUpperEdge(EquidistantGrid& g, double xmax)
{
  //A tiny tolerance ensures that a node which mathematically is at xmax is
  //always kept, rather than randomly dropped or kept depending on rounding
  //errors in xmax or in the node position:
  xmax += 1e-9 * g.binWidth;
#ifndef NDEBUG
  nc_assert(g.npts >= 2);
  nc_assert(std::isfinite(g.x0));
  nc_assert(std::isfinite(g.binWidth));
  nc_assert(g.binWidth > 0.0);
  nc_assert(std::isfinite(g.x1()));
  nc_assert(std::isfinite(xmax));
  EquidistantGrid expected = g;
  auto check = [&expected,&g] {
    return ( expected.x0 == g.x0
             && expected.binWidth == g.binWidth
             && expected.npts == g.npts );
  };
  nc_assert(check());
  while ( expected.npts > 2 && expected.x1()>xmax )
    --expected.npts;
#endif
  if ( g.npts <= 2 || g.x1() <= xmax ) {
    nc_assert(check());
    return;
  }

  if ( xmax <= g.x0 + g.binWidth ) {
    g.npts = 2;
    nc_assert(check());
    return;
  }

  const std::size_t orig_npts = g.npts;
  const double q = (xmax - g.x0) / g.binWidth;
  nc_assert( std::isfinite(q) && q >= 0.0 );
  if (q < 1.0) {
    g.npts = 2;
  } else if (q >= static_cast<double>(orig_npts - 1)) {
    nc_assert( g.npts == orig_npts );
  } else {
    const double k = std::floor(q);
    nc_assert(std::isfinite(k) && k >= 1.0);
    nc_assert(k < static_cast<double>(orig_npts - 1));
    const std::size_t npts = static_cast<std::size_t>(k) + 1;
    g.npts = orig_npts < npts ? orig_npts : npts;
  }

  // Correct rounding errors:
  while (g.npts > 2 && g.x1() > xmax)
    --g.npts;
  while (g.npts < orig_npts) {
    ++g.npts;
    if (g.x1() > xmax) {
      --g.npts;
      break;
    }
  }
  nc_assert(check());
}

NC::VectD NC::VDOS::mergeGridsWithTol( const VectD& a, const VectD& b,
                                       double rtol )
{
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
  auto nearlyEqual = [](double x, double y)
  {
    return ncabs(x - y) <= 1e-11 * ncmax(ncabs(x), ncabs(y));
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
    //A point which is the same as one in a, up to rounding errors, is not
    //kept. Notably the same node of a lattice might have been calculated in
    //different ways for the two grids. Without this, whether they end up
    //bitwise identical (and thus one of them is removed by add) or an ulp
    //apart (in which case both would be kept) would depend on rounding errors.
    if (!keep) {
      if (ia < a.size() && nearlyEqual(x, vectAt(a, ia)))
        continue;
      if (ia > 0 && nearlyEqual(x, vectAt(a, ia - 1)))
        continue;
    }
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

//Estimate where a locally-Gaussian-like curve on the equidistant grid
//x0+i*binwidth crosses yval, given the immediate bracket
//(idxLo,idxHi=idxLo+1) with min(spec[idxLo],spec[idxHi]) < yval <=
//max(...). High-order Gn spectra approach a Gaussian shape via the
//central limit theorem acting on repeated self-convolution, i.e.
//ln(spec) is close to QUADRATIC (not linear) in energy in the tails
//where this function is used -- so this fits a quadratic (ordinary
//least squares) to ln(spec) vs x over a small window extending up to
//nExtra points beyond the bracket on each side (clipped to the array
//bounds, and to spec>0 points, since non-positive values can't enter a
//log-space fit), and solves it for the root nearest the bracket. Using
//several points symmetrically, rather than trusting only the single
//point on each side of the naive bracket, reduces sensitivity to a
//last-ULP-level fluctuation landing on any one of those points
//(confirmed empirically in app_gnerange: with the quadratic model,
//injecting a 1e-3 relative perturbation at the two bracket points moves
//the estimated crossing by only ~1e-6 relative on a synthetic Gaussian
//test case). The result is always clamped to stay within the window's
//own x-extent, so an ill-conditioned or unlucky fit can never
//extrapolate far from the region it was fitted to, and the function
//falls back to the plain 2-point bracket (log-linear, or if that also
//fails, linear) interpolation whenever the windowed fit is degenerate
//(fewer than 5 usable points, non-finite/non-positive discriminant, or
//a fit that is not really quadratic):
double NC::VDOS::estimateSpectrumCrossing( double x0, double binwidth,
                                           Span<const double> spec,
                                           std::size_t idxLo, std::size_t idxHi,
                                           double yval, std::size_t nExtra )
{
  nc_assert( idxHi == idxLo + 1 );
  nc_assert( ncmin(spec[idxLo],spec[idxHi]) < yval
            && yval <= ncmax(spec[idxLo],spec[idxHi]) );
  auto xAt = [x0,binwidth](std::size_t i) { return VDOS::equidistantGridPoint(x0,binwidth,i); };

  auto twoPointFallback = [&]( bool logspace )
  {
    const double v0 = spec[idxLo];
    const double v1 = spec[idxHi];
    const double t = ( logspace && v0 > 0.0 && v1 > 0.0
                      ? (std::log(yval)-std::log(v0))/(std::log(v1)-std::log(v0))
                      : (yval-v0)/(v1-v0) );
    return nclerp( xAt(idxLo), xAt(idxHi), ncclamp(t,0.0,1.0) );
  };
  if ( !( spec[idxLo] > 0.0 && spec[idxHi] > 0.0 ) )
    return twoPointFallback(false);//can't do log-space at all

  const std::size_t wlo = ( idxLo >= nExtra ? idxLo-nExtra : 0 );
  const std::size_t whi = std::min<std::size_t>( idxHi+nExtra, spec.size()-1 );

  //Fit ln(spec) = a + b*xc + c*xc^2, xc=x-xmid (centred on the bracket
  //midpoint for conditioning), via the normal equations for an
  //ordinary least-squares quadratic fit:
  const double xmid = 0.5*( xAt(idxLo) + xAt(idxHi) );
  double S0(0.0), S1(0.0), S2(0.0), S3(0.0), S4(0.0);
  double T0(0.0), T1(0.0), T2(0.0);
  for ( auto i : ncrange(wlo,whi+1) ) {
    if ( !(spec[i]>0.0) )
      continue;
    const double x = xAt(i)-xmid;
    const double y = std::log(spec[i]);
    const double x2 = x*x;
    S0 += 1.0; S1 += x; S2 += x2; S3 += x2*x; S4 += x2*x2;
    T0 += y; T1 += x*y; T2 += x2*y;
  }
  if ( S0 < 5.0 )
    return twoPointFallback(true);
  //Solve [S0 S1 S2; S1 S2 S3; S2 S3 S4]*[a;b;c] = [T0;T1;T2] via
  //Cramer's rule:
  const double D  = S0*(S2*S4-S3*S3) - S1*(S1*S4-S3*S2) + S2*(S1*S3-S2*S2);
  if ( !( ncabs(D) > 0.0 ) )
    return twoPointFallback(true);
  const double Da = T0*(S2*S4-S3*S3) - S1*(T1*S4-S3*T2) + S2*(T1*S3-S2*T2);
  const double Db = S0*(T1*S4-S3*T2) - T0*(S1*S4-S3*S2) + S2*(S1*T2-T1*S2);
  const double Dc = S0*(S2*T2-T1*S3) - S1*(S1*T2-T1*S2) + T0*(S1*S3-S2*S2);
  const double a = Da/D, b = Db/D, c = Dc/D;
  const double target = std::log(yval) - a;//solve c*xc^2+b*xc-target=0
  double xcross;
  if ( !std::isfinite(c) || ncabs(c) < 1e-8*ncabs(b) ) {
    //Effectively linear (flat curvature in this window):
    if ( !( std::isfinite(b) && ncabs(b) > 0.0 ) )
      return twoPointFallback(true);
    xcross = target/b;
  } else {
    const double disc = b*b + 4.0*c*target;
    if ( !(disc >= 0.0) )
      return twoPointFallback(true);
    const double sq = std::sqrt(disc);
    const double r1 = (-b+sq)/(2.0*c);
    const double r2 = (-b-sq)/(2.0*c);
    //Both roots solve the fitted quadratic; pick whichever is nearest
    //the bracket (i.e. smallest in the centred coordinate), since the
    //other root is some unrelated, far-away crossing of the same
    //parabola:
    xcross = ( ncabs(r1) < ncabs(r2) ? r1 : r2 );
  }
  const double result = xcross + xmid;
  if ( !std::isfinite(result) )
    return twoPointFallback(true);
  return ncclamp( result, xAt(wlo), xAt(whi) );
}

void NC::VDOS::applyCrossingTaper( Span<double> spec, double xcross,
                                   bool risingEdge, double halfwidth )
{
#ifndef NDEBUG
  nc_assert_always( std::isfinite(xcross) );
  nc_assert_always( halfwidth > 0.0 );
#endif
  const double lo = xcross - halfwidth;
  const double hi = xcross + halfwidth;
  const double ilo_d = ncmax( 0.0, std::floor(lo) );
  const double ihi_d = ncmin( double(spec.size()-1), std::ceil(hi) );
  if ( ihi_d < ilo_d )
    return;
  const std::size_t ilo = static_cast<std::size_t>(ilo_d);
  const std::size_t ihi = static_cast<std::size_t>(ihi_d);
  const double invWidth = 1.0/(hi-lo);
  for ( auto i : ncrange(ilo,ihi+1) ) {
    const double t = ncclamp( ( double(i) - lo ) * invWidth, 0.0, 1.0 );
    //Quintic smootherstep (Ken Perlin): 0 and 1 derivatives vanish at both
    //ends, so the taper introduces no kink at the window edges either:
    const double s = t*t*t*(t*(t*6.0-15.0)+10.0);
    spec[i] *= ( risingEdge ? s : 1.0-s );
  }
}

NC::PairDD NC::VDOS::estimateGnErange( double egrid_lower, double egrid_binwidth,
                                       Span<const double> spec,
                                       double relcontriblvl )
{
#ifndef NDEBUG
  nc_assert_always( spec.size() >= 2 );
  nc_assert_always( relcontriblvl > 0.0 && relcontriblvl < 1.0 );
  nc_assert_always( std::isfinite(egrid_lower) && egrid_binwidth > 0.0 );
  nc_assert_always( *std::min_element(spec.begin(),spec.end()) >= 0.0 );
#endif
  auto xAt = [egrid_lower,egrid_binwidth](std::size_t i)
  { return equidistantGridPoint(egrid_lower,egrid_binwidth,i); };
  const double spec_max = *std::max_element( spec.begin(), spec.end() );
  const double threshold = relcontriblvl * spec_max;
  PairDD erange( xAt(0), xAt(spec.size()-1) );

  constexpr std::size_t nExtra = 2;//extra points on each side of the
                                   //immediate bracket to include in the
                                   //windowed log-linear fit.

  for ( auto e : enumerate(spec) ) {
    if ( e.val >= threshold ) {
      erange.first = ( e.idx == 0 )
        ? xAt(0)
        : estimateSpectrumCrossing( egrid_lower, egrid_binwidth, spec,
                            e.idx-1, e.idx, threshold, nExtra );
      break;
    }
  }
  for ( std::size_t i = spec.size(); i > 0; --i ) {
    if ( spec[i-1] >= threshold ) {
      const double x = ( i == spec.size() )
        ? xAt(i-1)
        : estimateSpectrumCrossing( egrid_lower, egrid_binwidth, spec,
                           i-1, i, threshold, nExtra );
      erange.second = ncmin( erange.second, x );
      break;
    }
  }
  nc_assert( erange.second >= erange.first );
  return erange;
}

double NC::VDOS::estimateFFTConvolutionNoiseFloor( double peak, std::size_t n,
                                                   double safetyFactor )
{
#ifndef NDEBUG
  nc_assert_always( peak >= 0.0 );
  nc_assert_always( n >= 1 );
  nc_assert_always( safetyFactor >= 0.0 );
#endif
  return ( safetyFactor * peak * std::numeric_limits<double>::epsilon()
          * std::sqrt( static_cast<double>(n) ) );
}
