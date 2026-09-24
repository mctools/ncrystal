#ifndef NCrystal_VDOSUtils_hh
#define NCrystal_VDOSUtils_hh

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

#include "NCrystal/internal/utils/NCRect.hh"
#include "NCrystal/internal/utils/NCSpan.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace VDOS {

    ///////////////////////////////////////////////////
    // Various utilities needed for VDOS processing. //
    ///////////////////////////////////////////////////

    //Interval where f(x) = x^n*exp(-x) is above eps*fpeak.
    PairDD rangeXNexpMX(unsigned n, double eps, double accuracy = 1e-13 );

    // Estimate where a tabulated, locally-Gaussian-like spectrum on the
    // equidistant grid x0+i*binwidth crosses yval, given the immediate
    // bracket (idxLo,idxHi=idxLo+1) with min(spec[idxLo],spec[idxHi]) <
    // yval <= max(...). Fits a quadratic (ordinary least squares) to
    // ln(spec) vs x over a window extending nExtra points beyond the
    // bracket on each side (clipped to the array bounds and to spec>0
    // points), solved for the root nearest the bracket, so the estimate
    // depends on several neighbouring points rather than trusting only the
    // single point on each side of the bracket -- much less sensitive to a
    // last-ULP-level fluctuation landing on any one of those points than a
    // plain two-point interpolation (confirmed empirically in app_gnerange:
    // injecting a 1e-3 relative perturbation at the bracket points moves
    // the estimate by only ~1e-6 relative on a synthetic Gaussian case).
    // The result is always clamped to the window's own x-extent, and falls
    // back to a plain two-point (log-linear, or linear) interpolation
    // whenever the windowed fit is degenerate (fewer than 5 usable points,
    // a non-finite/non-positive discriminant, or a fit that is not really
    // quadratic). Used by estimateGnErange below (with nExtra=2), and by
    // NCVDOSGn.cc's produceNewOrderByConvolutionImpl (with a larger nExtra)
    // to make the per-order truncation edge itself far less sensitive to
    // which side of a threshold a single noisy point happens to fall on --
    // see docs/claude_session_vdos_fma_reprod.md.
    double estimateSpectrumCrossing( double x0, double binwidth,
                                     Span<const double> spec,
                                     std::size_t idxLo, std::size_t idxHi,
                                     double yval, std::size_t nExtra );

    // Estimate the interval [x0,x1] outside of which a tabulated Gn spectrum
    // (density values at the equidistant grid points egrid_lower+i*
    // egrid_binwidth, i=0..spec.size()-1) is everywhere below
    // relcontriblvl*max(spec) (0<relcontriblvl<1). Takes the grid as
    // (lower,binwidth) rather than a materialised array, since it only ever
    // needs a handful of points from a small local window around each
    // crossing (computed on demand via equidistantGridPoint) -- avoiding an
    // O(spec.size()) allocation on every call, which matters since this is
    // called ~twice per phonon order from VDOSGn::eRange. See app_gnerange
    // and docs/claude_session_vdos_fma_reprod.md for the reproducibility
    // investigation this is part of.
    PairDD estimateGnErange( double egrid_lower, double egrid_binwidth,
                             Span<const double> spec,
                             double relcontriblvl );

    // Estimate a numerical noise floor for a spectrum of length n and peak
    // magnitude peak, as produced by FastConvolve's FFT-based linear
    // convolution (NCFastConvolve.hh). FFT round-off accumulates roughly as
    // peak*eps*sqrt(n) (a standard bound for a length-n discrete Fourier
    // transform); safetyFactor scales this into a value comfortably above
    // where convolution round-off is actually observed in practice
    // (calibrated in app_gnconvnoisefloor against a Neumaier-summed O(n^2)
    // reference convolution, both on synthetic spectra and on real Gn
    // spectra generated live from production VDOSEval/FastConvolve). NOT
    // YET wired into production: intended for raising the truncation
    // cutoff in NCVDOSGn.cc's produceNewOrderByConvolutionImpl above the
    // convolution's own noise floor, so the truncation edge there is not
    // decided by an isolated point that only clears a bare relative-to-peak
    // threshold by chance round-off. See
    // docs/claude_session_vdos_fma_reprod.md.
    double estimateFFTConvolutionNoiseFloor( double peak, std::size_t n,
                                             double safetyFactor = 8.0 );

    // Returns the intersection between the provided Rectangle in the alpha-beta
    // plane, and the kinematically available phasespace for a neutron of a
    // given E/kT (i.e. the set of points satisfying (alpha-beta)^2 <
    // 4*E_div_kT*alpha). The returned rectangle is empty when the intersection
    // has no area. Boundary-only intersections are therefore excluded, but of
    // course the usual caveats of floating point arithmetic applies.
    Rectangle findABExtentWithinKB( const Rectangle&, double E_div_kT );

    // In-place spacing out of the positive and negative parts of a finite,
    // strictly increasing grid. For adjacent nonzero points a < b with the same
    // sign, the points are considered too close when (0.0<=rtol<1.0):
    //
    //   b - a < rtol * min(abs(a), abs(b)).
    //
    // The first and last point of both positive and negative subgrids remain
    // fixed, while interior points are moved in place as needed. If there is
    // not enough room between fixed points to space out intermediate points, an
    // error is raised. In case of any errors, the grid might be left in an
    // invalid state.
    void spaceOutGrid(Span<double>, double rtol);

    // Attempt to add points to a grid of size < npts until it contains npts
    // values. At most one point is inserted into each interval of grid
    // points. Intervals not crossing 0 are considered in increasing order of
    // the magnitude of their nearest endpoint to zero, with positive-side
    // intervals preferred on ties. Intervals having an endpoint equal to zero
    // qualify regardless of rtol; other non-crossing intervals qualify when
    // their endpoint magnitudes differ by more than a factor of 1 + rtol. On
    // returning, the grid might still contain fewer than npts values if not
    // enough suitable original intervals are available (the function can of
    // course be called again if needed, possibly with a lower rtol).
    void topOffGrid( VectD& g, std::size_t npts, double rtol = 0.1 );

    // Position of point i of an equidistant grid: x0+i*binWidth, but calculated
    // with a single rounding (using an explicit std::fma). For consistency, all
    // node positions of equidistant grids and piecewise linear functions must
    // be calculated with this function.
    double equidistantGridPoint( double x0, double binWidth, std::size_t i );

    // The makeCommonGrid function merges several evenly spaced input grids into
    // a single grid.
    //
    // * Input grids may have different offsets (x0), point spacing (binWidth),
    //   and number of points (npts).
    // * Points closer than 10% of the smallest input bin width are merged.
    // * The output is not necessarily evenly spaced and thus requires
    //   description in a full vector of points.
    // * Its first and last points are always the min and max of input domains.
    struct EquidistantGrid final {
      double x0, binWidth;
      std::size_t npts;
      double xAt( std::size_t i ) const;//x{i}=x0+i*binWidth (safely)
      double x1() const { return xAt( npts - 1 ); }
    };
    VectD makeCommonGrid( Span<const EquidistantGrid> );


    // Merge two grids a and b into one. Taking all points of a and those pts in
    // b that are either an extreme point (highest or lowest of both grids) or
    // far enough from pts in a. Far enough means different signs (+/-/0), or
    // that the ratio of their absolute values is larger than 1+rtol.
    VectD mergeGridsWithTol( const VectD& a, const VectD& b,
                             double rtol = 0.01 );

    // Combines two overlapping grids with common binWidth into a single grid.
    // The returned grid starts at the lower origin of the grids
    // (i.e. x0=min(g1.x0,g2.x0)), and has just enough points for both original
    // domains to be covered by the new interval [x0,x1].
    EquidistantGrid coverEquidistantGrids( const EquidistantGrid&,
                                           const EquidistantGrid& );

    struct PWLFct final : private MoveOnly {
      // Utility struct representing a piecewise-linear function with values f
      // at the equidistant points x0+i*binWidth, and zero outside [x0,x1()]. It
      // is a move-only class, since it can optionally own the data that f
      // refers to.
      double x0 = 0.0; //x{i=0}
      double binWidth = 0.0;//distance between x{i} and x{i+1}
      Span<const double> f;//values of f at the x{i} points. The size of the
                           //span encodes the number of points.
      double xAt( std::size_t i ) const;//x{i}=x0+i*binWidth (safely)
      double x1() const { return xAt( f.size()-1 ); }
      VectD dataHolder;//optional, so can hold its data if needed.

      //Move semantics (re-points f if it refers to dataHolder):
      PWLFct() = default;
      PWLFct( PWLFct&& ) noexcept;
      PWLFct& operator=( PWLFct&& ) noexcept;
      Span<double> f_mutable();//edit data (must own data!)
    private:
      void moveDataFrom( PWLFct& ) noexcept;
    };

    // Narrow a PWL function to its part at positive x (new x0 will be above
    // tol*binWidth to zero). Returns an empty (PwlFct{}) function in case there
    // are not two grid points above the threshold. Otherwise the returned
    // function owns its data (in its dataHolder).
    PWLFct pwlNarrowToPos( const PWLFct&, double tol = 1e-3 );

    // Evaluates the weighted sum of functions on the supplied grid. Optionally
    // weights can be applied to each function (unit weights if empty).
    // For stability, grid points within a tiny tolerance (1e-9*binWidth) of the
    // endpoints of a function are considered to be exactly at the endpoint.
    VectD evalPWLSum( Span<const PWLFct> fs,
                      Span<const double> grid,
                      Span<const double> weights = {} );

    // Given a non-negative piecewise-linear function defined by by x and y,
    // this function removes (in-place) as many points from the back as possible
    // while keeping the discarded integral at most frac times the total.
    void trimTailByIntegral( VectD& x, VectD& y, double frac );

    // Remove points from the top of the grid until g.x1() <= xmax, but always
    // keeping at least 2 points. A tolerance of 1e-9*binWidth is applied, so a
    // node which mathematically is at xmax is kept regardless of rounding
    // errors.
    void trimEquidistantGridUpperEdge( EquidistantGrid& g, double xmax );

  }
}

////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace VDOS {
    inline PWLFct::PWLFct( PWLFct&& o ) noexcept
      : x0(o.x0), binWidth(o.binWidth), f(o.f)
    {
      moveDataFrom( o );
    }

    inline PWLFct& PWLFct::operator=( PWLFct&& o ) noexcept
    {
      if ( this != &o ) {
        x0 = o.x0;
        binWidth = o.binWidth;
        f = o.f;
        moveDataFrom( o );
      }
      return *this;
    }

    inline void PWLFct::moveDataFrom( PWLFct& o ) noexcept
    {
      const bool owns = ( !o.dataHolder.empty()
                          && o.f.data() == o.dataHolder.data() );
      dataHolder = std::move( o.dataHolder );
      o.dataHolder.clear();
      if ( owns )
        f = Span<const double>( dataHolder );
      o.f = Span<const double>();
    }

    inline Span<double> PWLFct::f_mutable()
    {
      nc_assert( dataHolder.size() == f.size() );
      return dataHolder;
    }

    inline double equidistantGridPoint( double x0, double binWidth,
                                        std::size_t i )
    {
      return std::fma( binWidth, static_cast<double>(i), x0 );
    }

    inline double PWLFct::xAt( std::size_t i ) const
    {
      return equidistantGridPoint( x0, binWidth, i );
    }

    inline double EquidistantGrid::xAt( std::size_t i ) const
    {
      return equidistantGridPoint( x0, binWidth, i );
    }

  }
}
#endif
