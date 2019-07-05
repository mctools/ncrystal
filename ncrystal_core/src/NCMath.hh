#ifndef NCrystal_Math_hh
#define NCrystal_Math_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2019 NCrystal developers                                   //
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

#include "NCrystal/NCDefs.hh"
#include <vector>

namespace NCrystal {

  const double constant_c  = 299792458e10;// speed of light in Aa/s
  const double constant_dalton2kg =  1.660539040e-27; // amu to kg (source: NIST/CODATA 2018)
  const double constant_dalton2eVc2 =  931494095.17; // amu to eV/c^2 (source: NIST/CODATA 2018)
  const double constant_avogadro = 6.022140857e23; // mol^-1 (source: NIST/CODATA 2018)
  const double constant_boltzmann = 8.6173303e-5;  // eV/K
  const double const_hhm = 4.144249671718981e-3; // hbar^2/neutron_mass
  const double const_neutron_mass = 1.674927471e-24; //gram
  const double const_neutron_atomic_mass = 1.00866491588; //atomic unit
  const double const_ekin_2200m_s = 0.02529886 ; //eV, neutron kinetic energy at 2200m/s.
  const double constant_planck = 4.135667662e-15 ;//[eV*s]
  const double constant_hbar = constant_planck*kInv2Pi; //[eV*s]

  //Our own min/max/abs to make sure only double versions are used:
  double ncmin(double, double);
  double ncmax(double, double);
  double ncclamp(double val, double low, double up);//returns val, clamped to interval [low,up]
  double ncabs(double);
  bool ncisnan(double);
  bool ncisinf(double);
  double nccopysign(double,double);//std::copysign from C++11, with fallback code for c++98.

  //Error function (only in cmath from C++11 onwards):
  double ncerf(double);

  //std::is_sorted (only in algorithm from C++11 onwards). Only support std::vector<double>
  bool ncis_sorted(std::vector<double>::const_iterator itb, std::vector<double>::const_iterator ite);

  //Generate numbers from thermal (Maxwell) distributions, passing T[Kelvin] or sqrt(T[Kelvin]):
  //(nb. genThermalNeutronEnergy doesn't actually depend on the particle being a neutron)
  double genThermalNeutronVelocity(double sqrt_temp_kelvin, double rand);//neutron velocity [m/s]
  double genThermalNeutronEnergy(double temp_kelvin, double rand);//neutron energy [eV]
  double genThermalNeutronWavelength(double sqrt_temp_kelvin, double rand);//neutron wavelength [Aa]
  double genThermalY(double rand);//Generate y=v/sqrt(kT/m) [dimensionless]

  //sinus/cosine options (indicated approximate timings from 2014 thinkpad with gcc 6.3.1):
  void sincos(double A,double&cosA, double& sinA);//slow cos(A) and sin(A) for any A   [40ns/call]
  void sincos_mpi8pi8(double A, double& cosA, double& sinA);//Requires -pi/8<=A<=pi/8  [ 6ns/call]
  void sincos_mpi2pi2(double A, double& cosA, double& sinA);//Requires -pi/2<=A<=pi/2  [ 9ns/call]
  void sincos_mpipi(double A, double& cosA, double& sinA);//Requires -pi<=A<=pi        [14ns/call]
  void sincos_0pi32(double A, double& cosA, double& sinA);//Requires 0<=A<=3pi/2.      [13ns/call]
  void sincos_02pi(double A, double& cosA, double& sinA);//Requires 0<=A<=2pi          [14.5ns/call]
  void sincos_0pi(double A, double& cosA, double& sinA);//same as sincos_0pi32, but asserts that 0<=A<=pi
  void sincos_mpi8pi8_lowres(double A, double& cosA, double& sinA);//As sincos_mpi8pi8 but absolute precision is at worst 5e-8 [3ns/call].
  double cos_02pi(double A);//Requires 0<=A<=2pi          [11ns/call]
  double cos_mpipi(double A);//Requires -pi<=A<=pi        [9.7ns/call]
  double cos_mpi2pi2(double A);//Requires -pi/2<=A<=pi/2  [7.6ns/call]
  double cos_mpi8pi8(double A);//Requires -pi/8<=A<=pi/8  [4.7ns/call]
  double sin_02pi(double A);//Requires 0<=A<=2pi          [8.8ns/call]
  double sin_mpipi(double A);//Requires -pi<=A<=pi        [8.3ns/call]
  double sin_mpi2pi2(double A);//Requires -pi/2<=A<=pi/2  [6.4ns/call]
  double sin_mpi8pi8(double A);//Requires -pi/8<=A<=pi/8  [4.7ns/call]

  //Approximations for exponential function (errors is 0.7e-10 or better):
  double exp_negarg_approx(double x);//error smaller than 0.7e-10, negative arguments only (slightly faster)
  double exp_approx(double x);//error smaller than 0.7e-10
  double exp_smallarg_approx(double x);//7th order Taylor expansion
  double atan_smallarg_approx(double x);//9th order Taylor expansion. Error smaller than 1e-5 for |x|<0.442, 1e-3 for |x|<0.68, 1e-2 for |x|<0.85.
  double atan_approx(double x);//calling atan_smallarg_approx when |x|<0.442 and falling back to std::atan and exact results otherwise.
  double expm1_smallarg_approx(double x);//7th order Taylor expansion

  //numerical integration
  void gauleg_10_ord(const double x1, const double x2, std::vector<double>& x, std::vector<double>& w);
  void gauloba_10_ord(const double x1, const double x2, std::vector<double>& x, std::vector<double>& w);
  double gaulobatto_grid(unsigned bz,double lower_beta, double upper_beta, std::vector<double>& beta, std::vector<double>& totwgt);
  double gaulobatto_grid(const std::vector<double>& old_beta, std::vector<double>& beta, std::vector<double>& totwgt);

  //vector operation
  double trapz(const std::vector<double> &y,const std::vector<double> &x);
  void flip(const std::vector<double> & arr, std::vector<double> & arr_dst, bool opposite_sign=false);
  std::vector<double> logspace(double start, double stop, unsigned num);
  std::vector<double> linspace(double start, double stop, unsigned num);
  double simpsons_irregular(const std::vector<double> &y, const std::vector<double> &x); //odd number of points is preferred
  void concatenate(std::vector<double>& arr_front,const std::vector<double>& arr_back, unsigned skip_pos=0);

  //misc:
  bool isPrime(unsigned n);//simple O(n^0.5) implementation
  bool intervalsOverlap(double a0, double b0, double a1, double b1);
  bool intervalsDisjoint(double a0, double b0, double a1, double b1);
  bool valueInInterval(double a, double b, double x);

  class Fct1D {
  public:
    //A very basic function object (probably will be replaced with c++11 lambdas
    //once c++98 support is dropped).
    Fct1D(){};
    virtual ~Fct1D();
    virtual double eval(double x) const = 0;
  };

  //estimate f'(x) derivatives numerically, based on Taylor expansions and
  //Richardson's extrapolation formulas, by evaluating the function at points
  //within |h| of x. In the absense of numerical uncertainties, the error one
  //the returned estimates should be O(h^order), but estimates based on higher
  //orders or very low values of h might suffer from numerical uncertainties.
  double estimateDerivative(const Fct1D*, double x, double h = 1e-4, unsigned order=6);//order must be 4 or 6
  double estimateSingleSidedDerivative(const Fct1D*, double x, double h, unsigned order=4);//order must be 1,2,3 or 4, sign of h gives evaluation side.

  //Root finding:
  double findRoot(const Fct1D*f,double a, double b, double acc = 1e-13);

  class CosSinGridGen {
  public:
    //Helper class for looping through a set of equidistant points and getting
    //the cos/sin values with the use of sin/cos addition formulas, rather than
    //trigonometric calls (except for every 128 points, to prevent numerical
    //drifting). Use like this to loop over n points uniformly spaced between
    //offset and offset+(n-1)*delta (thus to loop from phimin to phimax, put
    //offset=phimin and delta = (phimax-phimin)/(n-1):
    //
    //
    //   CosSinGridGen grid(n,offset,delta);
    //   do {
    //     std::cout <<grid.current_cosval()<< " "<<grid.current_sinval()<<std::endl;
    //   } while (grid.step());
    //
    //If the small_angles flag is enabled, offset and delta must be in [0,3pi/2]
    //and [0,pi/8] respectively, which allows faster initialisation.
    CosSinGridGen( unsigned n, double offset, double delta, bool small_angles = true );
    ~CosSinGridGen(){}
    bool step();
    double current_cosval() const { return m_c; }
    double current_sinval() const { return m_s; }
  private:
    static const unsigned recalcfreq = 128;
    double m_c, m_s;
    double m_cd, m_sd;
    unsigned m_left, m_recalc;
    double m_phimax, m_negdelta;
  };

  class StableSum {
  public:
    //Numerically stable summation, based on Neumaier's
    //algorithm (doi:10.1002/zamm.19740540106).
    StableSum();
    ~StableSum();
    void add(double x);
    double sum() const;
  private:
    double m_sum, m_correction;
  };

}


////////////////////////////
// Inline implementations //
////////////////////////////

#if __cplusplus >= 201103L

//Notice that the usage of e.g. std::min rather than std::fmin, and the order
//switching of parameters passed to those functions is rather important, as it
//has to be just like this in order to end with one instruction. See
//https://godbolt.org/g/v4dyuW).
inline double NCrystal::ncmin(double a, double b) { return std::min(b,a); }
inline double NCrystal::ncmax(double a, double b) { return std::max(b,a); }
inline double NCrystal::ncabs(double a) { return std::abs(a); }
inline bool NCrystal::ncisnan(double a) { return std::isnan(a); }
inline bool NCrystal::ncisinf(double a) { return std::isinf(a); }
#else
inline double NCrystal::ncmin(double a, double b) { return a < b ? a : b; }
inline double NCrystal::ncmax(double a, double b) { return a > b ? a : b; }
inline double NCrystal::ncabs(double a) { return a < 0 ? -a : a; }
inline bool NCrystal::ncisnan(double a) {
#  if defined(__GNUC__) || defined(__clang__)
  return __builtin_isnan(a);
#  else
  return a!=a;//in principle only guaranteed for IEEE numbers and might break with -fastmath
#  endif
}
inline bool NCrystal::ncisinf(double a) {
#  if defined(__GNUC__) || defined(__clang__)
  return __builtin_isinf(a);
#  else
  return ncabs(a) > std::numeric_limits<double>::max();
#  endif
}
#endif

inline double NCrystal::nccopysign(double x,double y)
{
#if __cplusplus >= 201103L
  return std::copysign(x,y);
#else
  return ncabs(x)*(y<0?-1.0:1.0);
#endif
}


inline double NCrystal::ncclamp(double val, double low, double up)
{
  nc_assert(up>=low);
  //NB: C++17 has std::clamp, not sure if it is more efficient. The version here
  //was checked to produce just two instructions (do NOT change argument order
  //here, or it will generate more instructions).
  return ncmin(ncmax(val,low),up);
}

#if __cplusplus >= 201103L
  inline double NCrystal::ncerf(double x) { return std::erf(x); }
#elif defined(__GNUC__) || defined(__clang__)
  //oups - no guarantee for erf in C++98, but clang and gcc seems to provide a
  //builtin (not disabled by -fno-builtin):
  inline double NCrystal::ncerf(double x) { return __builtin_erf(x); }
#else
  //double oups. Attempt fallback and hope it came from somewhere by e.g. the
  //cmath header due to C99 support:
  inline double NCrystal::ncerf(double x) { return erf(x); }
#endif

inline double NCrystal::genThermalNeutronVelocity(double sqrt_temp_kelvin, double rand)
{
  return (90.791228875348494 * sqrt_temp_kelvin) * genThermalY(rand);
}

inline double NCrystal::genThermalNeutronEnergy(double temp_kelvin, double rand)
{
  double tmp = genThermalY(rand);
  return (4.3086714999999998e-5 * temp_kelvin) * tmp * tmp;
}

inline double NCrystal::genThermalNeutronWavelength(double sqrt_temp_kelvin, double rand)
{
  //ncmax essentially to avoid division by zero when rand=0:
  return 43.572864563953402 / ( sqrt_temp_kelvin * ncmax(genThermalY(rand),1.0e-300) );
}

inline void NCrystal::sincos(double A,double&cosA, double& sinA)
{
  cosA = std::cos(A);
  sinA = std::sin(A);
}

inline void NCrystal::sincos_mpi8pi8_lowres(double A, double& cosA, double& sinA) {
  nc_assert(ncabs(A)<=kPi*0.12500001);
  double mx2 = -A*A;
  sinA = A*(1.0 + mx2 * ( 1.66666666666666666666666666666666666666666667e-1 // - x^3 / 3!
                + mx2 * ( 8.33333333333333333333333333333333333333333333e-3 // + x^5 / 5!
                + mx2 * ( 1.98412698412698412698412698412698412698412698e-4 // - x^7 / 7!
                           ))));
  cosA = 1.0 + mx2 * ( 0.5 // - x^2 / 2!
             + mx2 * ( 4.16666666666666666666666666666666666666666667e-2 // + x^4 / 4!
             + mx2 * ( 1.38888888888888888888888888888888888888888889e-3 // - x^6 / 6!
                       )));
}

inline void NCrystal::sincos_mpipi(double A, double& cosA, double& sinA) {
  nc_assert(A>=-kPi&&A<=kPi);
  //Use abs/min/copysign tricks to actually do the evaluation in [-pi/2,pi/2].
  double Aabs = ncabs(A);
  nc_assert(Aabs<=kPi);
  sincos_mpi2pi2(ncmin(Aabs,kPi-Aabs),cosA,sinA);
  cosA = nccopysign(cosA,kPiHalf-Aabs);
  sinA = nccopysign(sinA,A);
}

inline void NCrystal::sincos_0pi32(double A, double& cosA, double& sinA) {
  //Like sincos_mpipi, but without need to deal with negative A's. Works until 3pi/2.
  nc_assert(A>=0.0&&A<=1.5*kPi);
  sincos_mpi2pi2(ncmin(A,kPi-A),cosA,sinA);
  cosA = nccopysign(cosA,kPiHalf-A);
}

inline void NCrystal::sincos_0pi(double A, double& cosA, double& sinA) {
  //Like sincos_mpipi, but without need to deal with negative A's. Works until 3pi/2.
  nc_assert(A>=0.0&&A<=kPi);
  sincos_mpi2pi2(ncmin(A,kPi-A),cosA,sinA);
  cosA = nccopysign(cosA,kPiHalf-A);
}

inline void NCrystal::sincos_02pi(double A, double& cosA, double& sinA) {
  nc_assert(A>=0.0&&A<=k2Pi);
  //Inlining sincos_mpipi code for A->A-kPi, with a final sign flip for both
  //output values.
  A -= kPi;
  double Aabs = ncabs(A);
  nc_assert(Aabs<=kPi);
  sincos_mpi2pi2(ncmin(Aabs,kPi-Aabs),cosA,sinA);
  cosA = nccopysign(cosA,Aabs-kPiHalf);
  sinA = nccopysign(sinA,-A);
}

inline double NCrystal::cos_02pi(double A) {
  nc_assert(A>=0.0&&A<=k2Pi);
  A-=kPi;
  return nccopysign(cos_mpipi(A),ncabs(A)-kPiHalf);
}

inline double NCrystal::sin_02pi(double A) {
  nc_assert(A>=0.0&&A<=k2Pi);
  A -= kPi;
  double Aabs = ncabs(A);
  double Sabs = sin_mpi2pi2(ncmin(Aabs,kPi-Aabs));
  return nccopysign(Sabs,-A);
}

inline bool NCrystal::intervalsOverlap(double a0, double b0, double a1, double b1)
{
  nc_assert(a0<=b0&&a1<=b1);
  return a1<=b0 && a0 <= b1;
}

inline bool NCrystal::intervalsDisjoint(double a0, double b0, double a1, double b1)
{
  nc_assert(a0<=b0&&a1<=b1);
  return a1 > b0 || a0 > b1;
}

inline bool NCrystal::valueInInterval(double a, double b, double x)
{
  nc_assert(b>=a);
  nc_assert(bool((x-a)*(x-b) <= 0.0) == bool(x>=a&&x<=b));
  return (x-a)*(x-b) <= 0.0;
}

inline double NCrystal::exp_smallarg_approx( double x )
{
  //7th order taylor expansion
  return 1.0+x*(1+x*(0.5+x*(0.16666666666666666666666666666666666667+x*(0.04166666666666666666666666666666666667
           +x*(0.00833333333333333333333333333333333333+x*(0.00138888888888888888888888888888888889
           +x*0.00019841269841269841269841269841269841))))));
}

inline double NCrystal::expm1_smallarg_approx( double x )
{
  //7th order taylor expansion
  return x*(1+x*(0.5+x*(0.16666666666666666666666666666666666667+x*(0.04166666666666666666666666666666666667
           +x*(0.00833333333333333333333333333333333333+x*(0.00138888888888888888888888888888888889
           +x*0.00019841269841269841269841269841269841))))));
}

inline double NCrystal::exp_negarg_approx( double x )
{
  nc_assert(x<=0.0);
  if (x<-0.1) {
    //Bring into range where approximation below is good, by using:
      // exp(x) = exp(x/n)^n with n=2^8
    double y = exp_negarg_approx(x*0.00390625);
    y*=y; y*=y; y*=y; y*=y;
    y*=y; y*=y; y*=y; y*=y;
    return y;
  }
  //7th order taylor expansion:
  return exp_smallarg_approx(x);
}

inline double NCrystal::exp_approx(double x)
{
  return x>0.0 ? 1.0/exp_negarg_approx(-x) : exp_negarg_approx(x);
}

inline double NCrystal::atan_smallarg_approx(double x) {
  //Taylor expansion: x-x^3/3+x^5/5-x^7/7+x^9/9 (next omitted term is -x^11/11):
  double x2 = x*x;
  return x*(1.0+x2*(-0.333333333333333333333333333333333333+x2*(0.2
       +x2*(-0.142857142857142857142857142857142857+0.111111111111111111111111111111111111*x2))));
}

inline double NCrystal::atan_approx(double x)
{
  //arctan, valid to 5 significant digits everywhere
  double x2 = x*x;
  if (x2<0.195364) {
    //Taylor expansion: x-x^3/3+x^5/5-x^7/7+x^9/9 (next omitted term is -x^11/11):
    return x*(1.0+x2*(-0.333333333333333333333333333333333333+x2*(0.2
       +x2*(-0.142857142857142857142857142857142857+0.111111111111111111111111111111111111*x2))));
  }
  return std::atan(x);
}

inline NCrystal::CosSinGridGen::CosSinGridGen( unsigned n, double offset, double delta, bool small_angles )
  : m_left(n-1),
    m_recalc(((recalcfreq-1)+(n/recalcfreq)*recalcfreq)-n),
    m_phimax(offset+(n-1)*delta),
    m_negdelta(-delta)
{
  nc_assert(n>0);
  if (small_angles) {
    sincos_0pi32(offset,m_c,m_s);
    sincos_mpi8pi8(delta,m_cd,m_sd);
  } else {
    sincos(offset,m_c,m_s);
    sincos(delta,m_cd,m_sd);
  }
}

inline bool NCrystal::CosSinGridGen::step() {
  if (!m_left)
    return false;
  --m_left;
  if ((m_left+m_recalc)%recalcfreq) {
    //Advance using addition formulas (aka rotations) instead of expensive
    //trigonometric calls:
    //       cos(val+delta) = cos(val)*cos(delta)-sin(val)*sin(delta)
    //       sin(val+delta) = cos(val)*sin(delta)+sin(val)*cos(delta)
    double c = m_c*m_cd-m_s*m_sd;
    m_s = m_c*m_sd+m_s*m_cd;
    m_c = c;
  } else {
    //For numerical stability we combat drifting values by occasionally
    //recalculating precise values:
    double v = m_phimax+m_negdelta*m_left;
    sincos(v,m_c,m_s);
  }
  return true;
}

inline NCrystal::StableSum::StableSum()
  : m_sum(0.0), m_correction(0.0)
{
}

inline NCrystal::StableSum::~StableSum()
{
}

inline void NCrystal::StableSum::add( double x )
{
  double t = m_sum + x;
  m_correction += ncabs(m_sum)>=ncabs(x)  ? (m_sum-t)+x : (x-t)+m_sum;
  m_sum = t;
}

inline double NCrystal::StableSum::sum() const
{
  return m_sum + m_correction;
}

#endif
