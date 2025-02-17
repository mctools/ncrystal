#ifndef NCrystal_Math_hh
#define NCrystal_Math_hh

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

#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/internal/utils/NCSpan.hh"
#include "NCrystal/internal/utils/NCRomberg.hh"
#include <functional>

namespace NCRYSTAL_NAMESPACE {

  //NB: For physics constants, also see NCDefs.hh.

  //Primary constants [NB: Some replicated in Python interface!]:
  constexpr double constant_c  = 299792458e10;// speed of light in Aa/s
  constexpr double constant_dalton2eVc2 =  931494095.17; // amu to eV/c^2 (source: NIST/CODATA 2018)
  constexpr double constant_avogadro = 6.02214076e23; // mol^-1 (source: NIST/CODATA 2018)
  constexpr double constant_dalton2gpermol = constant_dalton2kg*constant_avogadro*1000.0; // dalton to gram/mol
  //NB: constant_dalton2gmol is almost but not quite unity (cf. https://doi.org/10.1007/s00769-013-1004-9)

  constexpr double const_neutron_atomic_mass = const_neutron_mass_amu; // [amu]//obsolete name
  constexpr double constant_planck = 4.135667662e-15 ;//[eV*s]
  constexpr double constant_gas_R = 8.31446261815324;// [J/(K*mol) = m^3*Pa/(K*mol) = kg*m^2/(K*mol*s^2)]. Value is exact definition.

  //Derived values:
  constexpr double const_neutron_mass_evc2 = const_neutron_mass_amu * constant_dalton2eVc2 / (constant_c*constant_c);// [ eV/(Aa/s)^2 ]
  constexpr NeutronEnergy const_ekin_2200m_s = NeutronEnergy{0.5 * const_neutron_mass_evc2 * 2200.0e10 * 2200.0e10}; //neutron kinetic energy at 2200m/s [eV]
  constexpr double constant_hbar = constant_planck*kInv2Pi; //[eV*s]
  constexpr double const_hhm = constant_hbar*constant_hbar/const_neutron_mass_evc2;// hbar^2/neutron_mass [ eV*Aa^2 ]
  constexpr double constant_ekin2v = constexpr_sqrt(2.0/const_neutron_mass_evc2);//multiply this with sqrt(ekin[eV]) to get velocity in Aa/s

  //Our own min/max/abs to make sure only double versions are used:
  double ncmin(double, double);
  double ncmax(double, double);
  double ncclamp(double val, double low, double up);//returns val, clamped to interval [low,up]
  double ncclamp(double val, PairDD low_and_up);
  double ncabs(double);
  bool ncisnan(double);
  bool ncisinf(double);
  bool ncisnanorinf(double);

  //Versions with multiple parameters:
  double ncmin(double, double, double);
  double ncmax(double, double, double);
  double ncmin(double, double, double, double);
  double ncmax(double, double, double, double);

  //Various constexpr versions (not efficient for runtime usage):
  //NB: See NCDefs.hh for max/gdc/lcm functions:
  template<class TInt>
  constexpr TInt ncconstexpr_ispow2( TInt );
  template<class TInt>
  constexpr TInt ncconstexpr_roundupnextpow2( TInt );
  constexpr unsigned ncconstexpr_log10ceil(unsigned);

  //Check that span contains values that could be a grid. I.e. is non-empty,
  //sorted, no duplicated values, no NaN/inf's.
  bool nc_is_grid(Span<const double>);

  //sinus/cosine options (indicated approximate timings from 2014 thinkpad with gcc 6.3.1):
  void sincos(double A,double&cosA, double& sinA);//slow cos(A) and sin(A) for any A   [40ns/call]  NB: See also sincos_fast below!!!!
  void sincos_mpi256pi256(double A, double& cosA, double& sinA);//Requires -pi/256<=A<=pi/256 [ 3ns/call]
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

  //Functions which remaps arguments and calls sincos_02pi for a potential
  //factor of 3 speedup. To be conservative we do not retrofit this into the
  //sincos(..)  function above (yet).:
  struct sincosval_t { double sin, cos; };
  sincosval_t sincos_2pix( double );//Returns { sin(2*pi*x), cos(2*pi*x) }
  sincosval_t sincos_fast( double );//Returns { sin(x), cos(x) }

  //Approximations for exponential function (errors is 0.7e-10 or better):
  double exp_negarg_approx(double x);//error smaller than 0.7e-10, negative arguments only (slightly faster)
  double exp_approx(double x);//error smaller than 0.7e-10
  double exp_smallarg_approx(double x);//7th order Taylor expansion
  double atan_smallarg_approx(double x);//9th order Taylor expansion. Error smaller than 1e-5 for |x|<0.442, 1e-3 for |x|<0.68, 1e-2 for |x|<0.85.
  double atan_approx(double x);//calling atan_smallarg_approx when |x|<0.442 and falling back to std::atan and exact results otherwise.
  double expm1_smallarg_approx(double x);//7th order Taylor expansion

  //Evaluate erfc(a)-erfc(b) in a relatively numerically safe
  //manner and with as few actual calls to std::erfc as possible:
  double erfcdiff(double a, double b);

  //Returns exp(b)*erfc(x), which has the advantage that when b ~= x^2, it
  //cancels out the exp(-x^2) behaviour of erfc at high x, thus enabling
  //meaningful (rescaled) evaluations of erfc tails (which would otherwise be
  //loss due to numerical precision issues):
  double erfc_rescaled(double x, double b);

  //Evenly spaced points (like Numpy equivalent functions):
  VectD linspace(double start, double stop, unsigned num);
  VectD logspace(double start, double stop, unsigned num);
  VectD geomspace(double start, double stop, unsigned num);

  //misc:
  constexpr double constexpr_sqrt(double);//compile time sqrt
  inline constexpr double ncsquare( double x ) noexcept { return x*x; }
  inline constexpr double nccube( double x ) noexcept { return x*x*x; }
  bool isPrime(unsigned n);//simple O(n^0.5) implementation
  bool intervalsOverlap(double a0, double b0, double a1, double b1);
  bool intervalsOverlap( const PairDD&, const PairDD& );
  bool intervalsDisjoint(double a0, double b0, double a1, double b1);

  //Quick check that a<=x<=b (do not use if a or b might be infinite):
  bool valueInInterval(double a, double b, double x);
  bool valueInInterval( const PairDD& ab, double x);

  class Fct1D {
  public:
    //A very basic function object based on dynamic polymorphism (we can replace
    //with c++11 lambdas now c++98 support has been dropped).
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

  template <class Func>
  double findRoot2(Func&& f,double a, double b, double acc = 1e-13);

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
    static constexpr unsigned recalcfreq = 128;
    double m_c, m_s;
    double m_cd, m_sd;
    unsigned m_left, m_recalc;
    double m_phimax, m_negdelta;
  };

  class StableSum {
  public:
    //Numerically stable summation, based on Neumaier's
    //algorithm (doi:10.1002/zamm.19740540106).
    void add(double x);
    double sum() const;
  private:
    double m_sum = 0.0, m_correction = 0.0;
  };

  template<typename T, typename Container>
  inline bool hasValue(const Container & c, const T & value)
  {
    using std::begin;
    using std::end;
    auto itE = end(c);
    return std::find(begin(c), itE, value) != itE;
  }

  VectD::const_iterator findClosestValInSortedVector(const VectD& v, double value);

  //Test equality of floating point numbers, within relative and absolute tolerances:
  bool floateq(double a, double b, double rtol=1.0e-6, double atol=1.0e-6);

  //Numerical integration (see NCRomberg.hh for more flexible Romberg integration):
  template <class Func>
  double integrateSimpsons(const Func& f, double a, double b, unsigned n);
  template <class Func>
  double integrateTrapezoidal(const Func& f, double a, double b, unsigned n);
  template <class Func>
  double integrateRomberg17(Func&& f, double a, double b);
  template <class Func>
  double integrateRomberg33(Func&& f, double a, double b);
  template <class Func>
  double integrateRombergFlex(Func&& f, double a, double b,
                              double prec=1e-12, unsigned minlvl = 3, unsigned maxlvl = 10 );


  //Reduce pts on curve by removing points that are least important for overall shape:
  std::pair<VectD,VectD> reducePtsInDistribution( const VectD& x, const VectD& y, std::size_t targetN );

  //Vector utilities:
  inline void vectorAppend(VectD& v1, const VectD& v2);//appends contents of v2 to v1
  template<class TVector, class Func>
  inline TVector vectorTrf(const TVector&, const Func&);//create new vector of same type with function applied to all elements
  //Access vector contents with .at() in debug builds and [] in optimised builds
  template <class TVector>
  typename TVector::value_type& vectAt(TVector& v, typename TVector::size_type idx);
  template <class TVector>
  const typename TVector::value_type& vectAt(const TVector& v, typename TVector::size_type idx);

  //Hashing utils (std::size_t is chosen for compatibility with std::hash):
  typedef std::size_t HashValue;
  template <class T>
  HashValue calcHash(const T&);
  template <class T>
  void hash_combine(HashValue& seed, const T&);
  template <class TContainer>
  HashValue hashContainer(const TContainer& v);

}


////////////////////////////
// Inline implementations //
////////////////////////////

//Notice that the usage of e.g. std::min rather than std::fmin, and the order
//switching of parameters passed to those functions is rather important, as it
//has to be just like this in order to end with minimal set of instructions (in
//any clang or gcc-trunk April-2020: just 1 per std::min/std::max call, no more):
inline double NCrystal::ncmin(double a, double b) { return std::min(b,a); }
inline double NCrystal::ncmax(double a, double b) { return std::max(b,a); }
inline double NCrystal::ncmax(double a, double b, double c) { return std::max(c,std::max(b,a)); }
inline double NCrystal::ncmin(double a, double b, double c) { return std::min(c,std::min(b,a)); }
inline double NCrystal::ncmin(double a, double b, double c, double d) { return std::min(std::min(d,c),std::min(b,a)); }
inline double NCrystal::ncmax(double a, double b, double c, double d) { return std::max(std::max(d,c),std::max(b,a)); }
inline double NCrystal::ncclamp(double val, PairDD low_and_up) { return ncclamp(val,low_and_up.first,low_and_up.second); }

inline constexpr unsigned NCrystal::ncconstexpr_log10ceil( unsigned val )
{
  return val < 10u ? 1u : 1u + ncconstexpr_log10ceil( val / 10u );
}

namespace NCRYSTAL_NAMESPACE {
  namespace detail {
    template<class TInt>
    inline static constexpr TInt ncconstexpr_ispow2_helper( TInt a, TInt k )
    {
      return a <= k ? (a==k) : ncconstexpr_ispow2_helper(a,k*2);
    }
  }
}
template<class TInt>
inline constexpr TInt NCrystal::ncconstexpr_ispow2( TInt a )
{
  return detail::ncconstexpr_ispow2_helper<TInt>(a, 1);
}

template<class TInt>
constexpr TInt NCrystal::ncconstexpr_roundupnextpow2( TInt a )
{
  return ncconstexpr_ispow2(a) ? a : ncconstexpr_roundupnextpow2( a + 1 );
}

inline double NCrystal::ncabs(double a) { return std::abs(a); }
inline bool NCrystal::ncisnan(double a) { return std::isnan(a); }
inline bool NCrystal::ncisinf(double a) { return std::isinf(a); }
inline bool NCrystal::ncisnanorinf(double a) { return std::isnan(a) || std::isinf(a); }

inline double NCrystal::ncclamp(double val, double low, double up)
{
  nc_assert(up>=low);
  //NB: C++17 has std::clamp, not sure if it is more efficient. The version here
  //was checked to produce just two instructions (do NOT change argument order
  //here, or it will generate more instructions).
  return ncmin(ncmax(val,low),up);
}

inline bool NCrystal::intervalsOverlap( const PairDD& x, const PairDD& y )
{
  return intervalsOverlap( x.first, x.second, y.first, y.second );
}

inline NCrystal::sincosval_t NCrystal::sincos_2pix( double x )
{
  //NB: We used to simply return PairDD but that triggers annoying gcc warning
  //on arm, like discussed on https://stackoverflow.com/questions/77729813
  //. Moving to a custom struct solves the issue.
  sincosval_t res;
  sincos_02pi((x-std::floor(x))*k2Pi,res.cos,res.sin);
  return res;
}

inline NCrystal::sincosval_t NCrystal::sincos_fast( double x )
{
  return sincos_2pix( x * kInv2Pi );
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

inline void NCrystal::sincos_mpi256pi256(double A, double& cosA, double& sinA)
{
  nc_assert(ncabs(A)<=kPi*0.003906250001);
  sincos_mpi8pi8_lowres(A,cosA,sinA);
}

inline void NCrystal::sincos_mpipi(double A, double& cosA, double& sinA) {
  nc_assert(A>=-kPi&&A<=kPi);
  //Use abs/min/copysign tricks to actually do the evaluation in [-pi/2,pi/2].
  double Aabs = ncabs(A);
  nc_assert(Aabs<=kPi);
  sincos_mpi2pi2(ncmin(Aabs,kPi-Aabs),cosA,sinA);
  cosA = std::copysign(cosA,kPiHalf-Aabs);
  sinA = std::copysign(sinA,A);
}

inline void NCrystal::sincos_0pi32(double A, double& cosA, double& sinA) {
  //Like sincos_mpipi, but without need to deal with negative A's. Works until 3pi/2.
  nc_assert(A>=0.0&&A<=1.5*kPi);
  sincos_mpi2pi2(ncmin(A,kPi-A),cosA,sinA);
  cosA = std::copysign(cosA,kPiHalf-A);
}

inline void NCrystal::sincos_0pi(double A, double& cosA, double& sinA) {
  //Like sincos_mpipi, but without need to deal with negative A's. Works until 3pi/2.
  nc_assert(A>=0.0&&A<=kPi);
  sincos_mpi2pi2(ncmin(A,kPi-A),cosA,sinA);
  cosA = std::copysign(cosA,kPiHalf-A);
}

inline void NCrystal::sincos_02pi(double A, double& cosA, double& sinA) {
  nc_assert(A>=0.0&&A<=k2Pi);
  //Inlining sincos_mpipi code for A->A-kPi, with a final sign flip for both
  //output values.
  A -= kPi;
  double Aabs = ncabs(A);
  nc_assert(Aabs<=kPi);
  sincos_mpi2pi2(ncmin(Aabs,kPi-Aabs),cosA,sinA);
  cosA = std::copysign(cosA,Aabs-kPiHalf);
  sinA = std::copysign(sinA,-A);
}

inline double NCrystal::cos_02pi(double A) {
  nc_assert(A>=0.0&&A<=k2Pi);
  A-=kPi;
  return std::copysign(cos_mpipi(A),ncabs(A)-kPiHalf);
}

inline double NCrystal::sin_02pi(double A) {
  nc_assert(A>=0.0&&A<=k2Pi);
  A -= kPi;
  double Aabs = ncabs(A);
  double Sabs = sin_mpi2pi2(ncmin(Aabs,kPi-Aabs));
  return std::copysign(Sabs,-A);
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

inline bool NCrystal::valueInInterval( const NCrystal::PairDD& ab, double x)
{
  return valueInInterval(ab.first,ab.second,x);
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

inline bool NCrystal::floateq(double a, double b, double rtol, double atol)
{
  nc_assert(!std::isnan(a));
  nc_assert(!std::isnan(b));
  if ( std::isinf(a) || std::isinf(b) ) {
    return a == b;
  } else {
    return ncabs(a-b) <= (0.5 * rtol) * ( ncabs(a) + ncabs(b) ) + atol;
  }
}

template <class Func>
inline double NCrystal::integrateSimpsons(const Func& f, double a, double b, unsigned n) {
  nc_assert(b>a&&n>1);
  if (n%2==0)
    n+=1;//ensure odd
  StableSum res;
  const double dx = (b-a)/(n-1);
  res.add(f(a)*0.5);
  const unsigned nm1(n-1);
  for (unsigned i = 1; i < nm1; i+=2)
    res.add(f(a + i * dx)*2.0);
  for (unsigned i = 2; i < nm1; i+=2)
    res.add(f(a + i * dx));
  res.add(f(b)*0.5);
  return res.sum()*dx*(2.0/3.0);
}

template <class Func>
inline double NCrystal::integrateTrapezoidal(const Func& f, double a, double b, unsigned n) {
  nc_assert(b>a&&n>1);
  StableSum res;
  const double dx = (b-a)/(n-1);
  res.add(f(a)*0.5);
  const unsigned nm1(n-1);
  for (unsigned i = 1; i < nm1; ++i)
    res.add(f(a + i * dx));
  res.add(f(b)*0.5);
  return res.sum()*dx;
}

template <class Func>
inline double NCrystal::integrateRomberg17(Func&& f, double a, double b)
{
  auto callf = [&f]( double x ) -> double { return f(x); };
  struct R17 final : public Romberg {
    decltype(callf) m_f;
    R17(decltype(callf)&& ff) : m_f(std::move(ff)) {}
    double evalFunc(double x) const override { return m_f(x); }
    bool accept(unsigned, double, double,double,double) const override { return true; }
  };
  return R17(std::move(callf)).integrate(a,b);
}

template <class Func>
inline double NCrystal::integrateRomberg33(Func&& f, double a, double b)
{
  auto callf = [&f]( double x ) -> double { return f(x); };
  struct R33 final : public Romberg {
    decltype(callf) m_f;
    R33(decltype(callf)&& ff) : m_f(std::move(ff)) {}
    double evalFunc(double x) const override { return m_f(x); }
    bool accept(unsigned lvl, double, double,double,double) const override { return lvl>4; }
  };
  return R33(std::move(callf)).integrate(a,b);
}

template <class Func>
inline double NCrystal::integrateRombergFlex(Func&& f, double a, double b, double prec, unsigned minlvl, unsigned maxlvl )
{
  auto callf = [&f]( double x ) -> double { return f(x); };
  struct RFlex final : public Romberg {
    decltype(callf) m_f;
    double m_prec;
    unsigned m_minlvl;
    unsigned m_maxlvl;
    RFlex(decltype(callf)&& ff, double pp,unsigned minl, unsigned maxl) : m_f(std::move(ff)), m_prec(pp), m_minlvl(minl), m_maxlvl(maxl) {}
    double evalFunc(double x) const override { return m_f(x); }
    bool accept(unsigned level, double prev_est, double est, double, double) const override
    {
      return level>=m_minlvl && (level>=m_maxlvl || floateq(prev_est,est,m_prec,0.0));
    }
  };
  nc_assert(maxlvl>=minlvl);
  nc_assert(minlvl>=3);
  nc_assert(prec<1.0&&prec>=0.0);
  return RFlex(std::move(callf),prec,minlvl,maxlvl).integrate(a,b);
}

template <class Func>
inline double NCrystal::findRoot2(Func&& f,double a, double b, double acc)
{
  auto callf = [&f]( double x ) -> double { return f(x); };
  struct FctWrap final : public Fct1D {
    decltype(callf) m_f;
    FctWrap(decltype(callf)&& ff) : m_f(std::move(ff)) {};
    double eval(double x) const final { return m_f(x); }
  } fwrap(std::move(callf));
  return findRoot(&fwrap,a,b,acc);
}

inline void NCrystal::vectorAppend( NCrystal::VectD& v1, const NCrystal::VectD& v2 )
{
  v1.reserve(v1.size()+v2.size());
  v1.insert( v1.end(), v2.begin(), v2.end() );
}

template<class TVector, class Func>
inline TVector NCrystal::vectorTrf(const TVector& input, const Func& f)
{
  TVector out;
  out.reserve(input.size());
  for ( auto e : input )
    out.emplace_back( f(e) );
  return out;
}

template <class TVector>
inline typename TVector::value_type& NCrystal::vectAt(TVector& v, typename TVector::size_type idx)
{
#ifndef NDEBUG
  return v.at(idx);
#else
  return v[idx];
#endif
}

template <class TVector>
inline const typename TVector::value_type& NCrystal::vectAt(const TVector& v, typename TVector::size_type idx)
{
#ifndef NDEBUG
  return v.at(idx);
#else
  return v[idx];
#endif
}

namespace NCRYSTAL_NAMESPACE {

  template <class T>
  inline HashValue calcHash(const T& t)
  {
    return std::hash<T>{}(t);
  }

  template <class T>
  inline void hash_combine(HashValue& seed, const T& v)
  {
    //Based on boost::hash_combine and discussions at stackoverflow. We could
    //easily pick another algorithm here:
    seed ^= calcHash(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
  }

  template <class TContainer>
  inline HashValue hashContainer(const TContainer& v)
  {
    HashValue seed(0);
    for ( auto e : v )
      hash_combine(seed,e);
    return seed;
  }

}

#endif
