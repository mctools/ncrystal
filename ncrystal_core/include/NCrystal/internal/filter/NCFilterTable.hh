#ifndef NCrystal_FilterTable_hh
#define NCrystal_FilterTable_hh

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

#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/internal/utils/NCStrView.hh"

//FIXME: Check whether some of the utilities here (e.g. reducePoints and
//evalTable) duplicate functionality available elsewhere in NCrystal.

namespace NCRYSTAL_NAMESPACE {

  namespace Filter {

    //Piecewise linear tables of the total cross section (scattering +
    //absorption, in barn per atom) of an isotropic material, as a function of
    //the neutron wavelength (in Angstrom). The intended use is to attenuate a
    //beam passing through a material (e.g. a filter or a window), without
    //caring about what happens to the scattered neutrons.
    //
    //The table covers wavelengths from 0 to wlmax. Beyond wlmax, the last
    //segment can be extrapolated linearly (clamped at 0). The value at 0 is the
    //limit of the cross section for wavelength -> 0.
    //
    //The table is created by evaluating the cross section at a dense set of
    //wavelengths (log-spaced, plus points just below and above each known
    //discontinuity), and then selecting a subset of these points such that
    //linear interpolation reproduces the cross section at all the dense points
    //within the tolerance. The result is verified at the midpoints between all
    //neighbouring dense points, and points are added to the dense sample where
    //needed, until the tolerance is met everywhere. Finally, the table is
    //verified at random points in each segment.
    //
    //The tolerance is relative to max(xs,1e-12 barn per atom). Without such a
    //(negligible) absolute floor, cross sections vanishing at wavelength 0
    //(e.g. without absorption and inelastic scattering) could not be
    //tabulated.
    //
    //Known discontinuities (Bragg edges, and the boundaries of the energy
    //domains of the processes) are represented by two table points with the
    //same wavelength: the first with the cross section just below the
    //discontinuity, the second with the cross section just above it.
    //
    //Features narrower than the dense sampling, which are not known in advance,
    //can not be found reliably by sampling. Therefore a warning is emitted for
    //physics processes not known to be free of such features (e.g. a future
    //process with absorption resonances). Exceptions are thrown if the cross
    //section is not finite or negative, if it does not converge for wavelength
    //-> 0, or if an unexpected discontinuity is found.

    enum class ReductionAlgo { Greedy, DouglasPeucker };

    struct TableParams {
      double wlmax = 500.0;//Angstrom
      double tol = 1e-3;//Max relative error of the cross section (see below)
      unsigned ndense = 0;//Log-spaced dense sample points (plus discontinuities). 0: auto (see autoNDense)
      bool discontinuities = true;//Include points at the known discontinuities
      ReductionAlgo algo = ReductionAlgo::Greedy;
    };

    struct Table {
      VectD wl;//Angstrom
      VectD xs;//barn/atom
      double numberDensity;//atoms/Angstrom^3
      std::size_t nDiscontinuities;//Known discontinuities (clusters) in range
      std::size_t nDense;//Number of dense sample points (after refinement)
      std::size_t nRefine;//Number of refinement iterations
      std::vector<std::string> processes;//Names of the (leaf) physics processes
    };

    Table createTable( const MatCfg&, const TableParams& );

    //The same, for a general function (e.g. for testing), with the
    //discontinuities given explicitly (numberDensity will be 0 and processes
    //empty):
    Table createTable( const std::function<double(double)>& xs_of_wl,
                       VectD discontinuities, const TableParams& );

    //The lowest wavelength of the log-spaced dense sample (Angstrom). Below
    //it, the table is based on the limit at wavelength 0 and refinement:
    constexpr double denseWavelengthMin() { return 1e-9; }

    //The number of log-spaced dense sample points used for a given tolerance
    //and wlmax, when TableParams::ndense is 0:
    unsigned autoNDense( double tol, double wlmax );

    //Reduce a dense sample (x must be sorted, but may contain pairs of
    //identical values to represent discontinuities) to a subset of points,
    //such that linear interpolation reproduces all y values within the
    //tolerance tol, relative to max(|y|,abs_floor). The indices of the
    //selected points are returned.
    std::vector<std::size_t> reducePoints( const VectD& x, const VectD& y,
                                           double tol, ReductionAlgo,
                                           double abs_floor = 0.0 );

    //Evaluate a table at a given wavelength: linear interpolation, where a
    //pair of identical wavelengths (a discontinuity) uses the second value for
    //wavelengths >= the discontinuity, and linear extrapolation of the last
    //segment (clamped at 0) beyond the table:
    double evalTable( const double* wl, const double* xs, std::size_t n,
                      double wavelength );

    //JSON query: ["filtertable",CFGSTR,OPTION1,OPTION2,...] where the options
    //are strings like "tol=1e-3", "wlmax=500", "ndense=0" (auto),
    //"discontinuities=1" and "algo=greedy" (or "algo=dp").
    //FIXME: Use streamJSONHugeArray for the arrays when it becomes available.
    using Query = SmallVector<StrView,8>;
    void JSONQuery( std::ostream&, const Query& );

  }

}

#endif
