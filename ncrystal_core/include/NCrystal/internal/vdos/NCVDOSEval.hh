#ifndef NCrystal_VDOSEval_hh
#define NCrystal_VDOSEval_hh

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

#include "NCrystal/interfaces/NCSABData.hh"

namespace NCRYSTAL_NAMESPACE {

  class VDOSEval final : private MoveOnly {
  public:

    // Class which provides evaluations of the phonon density function f(E),
    // based on the tabulated points in a VDOSData instance. It also provides
    // derivation of derived quantitied, including those determined by
    // integrations with f(E).
    //
    // The f(E) function is defined by linear interpolation within the provided
    // grid. For energies below the first grid point, the density values will be
    // assumed to scale as the square of the energy - an approximation which
    // should be suitable for solid materials. For energies above the last grid
    // point, densities will be assumed to be 0.
    //
    // The VDOSEval constructor assumes (and will error if not) that the
    // provided VDOS parameterisation has already been regularised (e.g. by the
    // regulariseVDOSGrid function further down).

    VDOSEval( const VDOSData& );
    ~VDOSEval();

    //Basic parameters (npts from emin to emax, npts_extended from 0 to emax):
    struct GridInfo { double emin, emax; unsigned npts, npts_extended; };
    GridInfo getGridInfo() const;

    Temperature temperature() const { return m_temperature; }
    double kT() const { return m_kT; }
    AtomMass elementMassAMU() const { return m_elementMassAMU; }

    //Evaluate f(E) with proper interpolation, extrapolation and normalisation:
    double eval(double energy) const;

    //f(E) is automatically normalised [0,emax]. For debugging purposes, the
    //original value of its integral over that region can be accessed:
    double originalIntegral() const { return m_originalIntegral; }

    //Estimate various parameters of interest through relevant numerical
    //integrations:
    double calcEffectiveTemperature() const;
    double calcGamma0() const;

    //Evaluate the \bar G_1 function (Sjolander1958, eq. II.29) which is used as
    //seed of the phonon expansion procedure of Sjolander1958 in the VDOSGn
    //class. The gamma0 value passed in must be the result of a call to
    //calcGamma0() above:
    double evalG1Symmetric( double energy, double gamma0 ) const;
    double evalG1Asymmetric( double energy, double gamma0 ) const;

    //For energy transfers beyond 200kT it is not numerically safe in usual
    //double precision to calculate G1Asymmetric by calculating G1Symmetric and
    //applying a detailed balance factor manually. Rather the evalG1Asymmetric
    //method should be used separately for both:
    static constexpr double numericallySafeG1SymmetricELimitInUnitsOfKT = 200.0;

    //Helper function for evaluating G1Asymmetric at both -energy and +energy,
    //in the safest + efficient way, and returning
    //{G1asym(-energy),G1asym(+energy)}. Assumes energy to be positive!
    PairDD evalG1AsymmetricAtEPair( double energy, double gamma0 ) const;

    //Estimate mean-squared-displacement for atoms of a given mass. The gamma0
    //value passed in must be the result of a call to calcGamma0() above:
    double getMSD( double gamma0 ) const;
    double getMSD() const { return getMSD( calcGamma0() ); }


    ///////////////////////////////////////////////////////////
    // Enable verbose output (default is disabled unless the //
    // NCRYSTAL_DEBUG_PHONON environment variable is set.    //
    ///////////////////////////////////////////////////////////

    static void enableVerboseOutput(bool status = true);
    static bool verboseOutputEnabled();

  private:
    VectD m_density;
    double m_emin, m_emax, m_k, m_binwidth, m_invbinwidth;
    double m_kT;
    Temperature m_temperature;
    AtomMass m_elementMassAMU;
    double m_originalIntegral;
    unsigned m_nptsExtended;
    template <class Fct, class TStableSum>
    void integrateBinsWithFunction(Fct, TStableSum& ) const;
  };

  //////////////////////////////////////////////////////////////////////////////////
  //
  // Utility function for reparameterising a user-provided VDOS to conform to the
  // form required by VDOSEval and VDOSGn.
  //
  // Accepts a VDOS parameterisation (number of egrid points must be same as
  // the number of points in the density vector, or it must be just 2 - which
  // is shorthand for linspace(egrid[0],egrid[1],density.size()).
  //
  // Returns regularised VDOS parameterisation: {regular_egrid,regular_density}.
  //
  // A regular VDOS is characterized by having equidistantly spaced egrid points
  // from emin to emax (both positive), which would coincide with zero if the
  // grid was extended downwards by a whole number of bins. The returned egrid
  // is always specified in the compact form with just two entries, {emin,emax}.
  //
  std::pair<VectD,VectD> regulariseVDOSGrid( const VectD& egrid, const VectD& density );

  //Check if grid is already regular within stated tolerance. If it is, returns
  //high-precision corrected emax value, otherwise returns 0.0.
  double checkIsRegularVDOSGrid( const VectD& egrid, const VectD& density, double tolerance=1e-6 );
  double checkIsRegularVDOSGrid( const PairDD& egrid_range, const VectD& density, double tolerance=1e-6 );

}

#endif
