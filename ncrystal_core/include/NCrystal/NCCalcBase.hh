#ifndef NCrystal_CalcBase_hh
#define NCrystal_CalcBase_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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

/////////////////////////////////////////////////////////////////////
// Base class for all calculators, providing common infrastructure //
// including random sampling of various of distributions.          //
/////////////////////////////////////////////////////////////////////

#include "NCrystal/NCRandom.hh"
#include <string>
#include <vector>
#include <limits>
#include <cmath>

namespace NCrystal {


  class RandomBase;

  //Base class for ref-counted objects, destructors of which should be protected
  //rather than public.

  class CalcBase  : public RCBase {
  public:
    CalcBase(const char * calculator_type_name);
    const char * getCalcName() const { return m_name.c_str(); }
    void setCalcName(const char * n) { m_name = n; }

    //Can be used to control and access the random number stream (but note that
    //often the setDefaultRandomGenerator(..) in NCRandom.hh will be an easier
    //way to do this globally):
    void setRandomGenerator(RandomBase* rg);
    RandomBase* getRandomGenerator() { return m_randgen.obj(); }

  protected:
    //Access random numbers in derived classes (uniformly distributed in
    //(0,1]. Note that the method is const, although obviously the random stream
    //will be changed as a result of calling this:
    double rand() const { return m_randgen.obj() ? m_randgen.obj()->generate() : initDefaultRand(); }
    //Convenience:
    double randIsotropicScatterAngle() const;
    double randIsotropicScatterMu() const;
    void randIsotropicDirection(double (&)[3]) const;//result will be unit vector
    void randDirectionGivenScatterAngle(double theta, const double(&in)[3], double(&out)[3]) const;//outdir will be unit vector
    void randDirectionGivenScatterAngle(double costh, double sinth, const double(&in)[3], double(&out)[3]) const;//outdir will be unit vector
    void randNorm(double&g1) const;//sample single value from unit Gaussian
    void randNorm(double&g1, double&g2) const;//sample two independent values from unit Gaussian.
    void registerSubCalc(CalcBase*);//when implementation owns other CalcBase instances, they must be registered
  private:
    std::vector<CalcBase*> m_subcalcs;
    std::string m_name;
    mutable RCHolder<RandomBase> m_randgen;
    double initDefaultRand() const;
  protected:
    virtual ~CalcBase();
  };


  //Utility functions for converting between neutron wavelength [Aa] and kinetic
  //energy [eV], and for providing infinity:
  double wl2ekin( double wl );
  double ekin2wl( double ekin );
  const double infinity = std::numeric_limits<double>::infinity();
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCrystal {

  //The constant 8.1804... in the functions wl2ekin and ekin2wl is based on the
  //equation "h^2 * c^2 / (2.0*m)", using CODATA Internationally recommended
  //2014 values of the fundamental physical constants
  //(http://physics.nist.gov/cuu/Constants/Table/allascii.txt):
  //
  //  h = 4.135667662e-15 [Ev*s] <- Planck constant in
  //  c = 299792458.0e10 [Aa/s]  <- speed of light in vacuum
  //  m = 939.5654133e6 [eV]     <- neutron mass energy equivalent
  //
  //  h^2 * c^2 / (2.0*m) = 0.081804209605330899

  inline double wl2ekin( double wl)
  {
    //Aangstrom to eV
    return wl ? 0.081804209605330899 / ( wl * wl ) : std::numeric_limits<double>::infinity();
  }

  inline double ekin2wl( double ekin)
  {
    //eV to Aangstrom
    return ekin ? sqrt( 0.081804209605330899 / ekin ) : std::numeric_limits<double>::infinity();
  }

}

#endif
