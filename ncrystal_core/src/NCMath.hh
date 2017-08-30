#ifndef NCrystal_Math_hh
#define NCrystal_Math_hh

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

#include "NCrystal/NCException.hh"
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace NCrystal {

  const double constant_boltzmann = 8.6173303e-5;  //eV/K
  const double const_hhm = 4.144249671718981e-3; //constant_hbar*constant_hbar/constant_neutron_mass
  const double const_neutron_mass = 1.674927471e-24; //gram
  const double const_neutron_atomic_mass = 1.00866491588; //atomic unit
  const double const_ekin_2200m_s = 0.02529886 ; //eV, neutron kinetic energy at 2200m/s.

  //Our own min/max/abs to make sure only double versions are used:
  double ncmin(double, double);
  double ncmax(double, double);
  double ncabs(double);
  bool ncisnan(double);
  bool ncisinf(double);
  //Error function (only in cmath from C++11 onwards):
  double ncerf(double);

  //Generate numbers from thermal (Maxwell) distributions, passing T[Kelvin] or sqrt(T[Kelvin]):
  //(nb. genThermalNeutronEnergy doesn't actually depend on the particle being a neutron)
  double genThermalNeutronVelocity(double sqrt_temp_kelvin, double rand);//neutron velocity [m/s]
  double genThermalNeutronEnergy(double temp_kelvin, double rand);//neutron energy [eV]
  double genThermalNeutronWavelength(double sqrt_temp_kelvin, double rand);//neutron wavelength [Aa]
  double genThermalY(double rand);//Generate y=v/sqrt(kT/m) [dimensionless]

  //sincos options:
  void sincos(double A,double&cosA, double& sinA);//slow cos(A) and sin(A) for any A
  void sincos_mpipi(double A,double&cosA, double& sinA);//fast cos(A) and sin(A) when -pi<=A<=pi
  void sincos_0pi(double A,double&cosA, double& sinA);//faster cos(A) and sin(A) when 0<=A<=pi
  void sincos_mpi2pi2(double A,double&cosA, double& sinA);//faster cos(A) and sin(A) when -pi/2<=A<=pi/2

  //numerical integration
  void gauleg_10_ord(const double x1, const double x2, std::vector<double>& x, std::vector<double>& w);

  //vector operation
  double trapz(const std::vector<double> &y,const std::vector<double> &x);
  void flip(const std::vector<double> & arr, std::vector<double> & arr_dst, bool opposite_sign=false);
  std::vector<double> logspace(double start, double stop, unsigned num);
  std::vector<double> linspace(double start, double stop, unsigned num);

}


////////////////////////////
// Inline implementations //
////////////////////////////

#if __cplusplus >= 201103L
inline double NCrystal::ncmin(double a, double b) { return std::fmin(a,b); }
inline double NCrystal::ncmax(double a, double b) { return std::fmax(a,b); }
inline double NCrystal::ncabs(double a) { return std::fabs(a); }
inline bool NCrystal::ncisnan(double a) { return std::isnan(a); }
inline bool NCrystal::ncisinf(double a) { return std::isinf(a); }
#else
inline double NCrystal::ncmin(double a, double b) { return a < b ? a : b; }
inline double NCrystal::ncmax(double a, double b) { return a > b ? a : b; }
inline double NCrystal::ncabs(double a) { return a < 0 ? -a : a; }
inline bool NCrystal::ncisnan(double a) {
#  if defined(__GNUC__) or defined(__clang__)
  return __builtin_isnan(a);
#  else
  return a!=a;//in principle only guaranteed for IEEE numbers and might break with -fastmath
#  endif
}
inline bool NCrystal::ncisinf(double a) {
#  if defined(__GNUC__) or defined(__clang__)
  return __builtin_isinf(a);
#  else
  return ncabs(x) > std::numeric_limits<double>::max();
#  endif
}
#endif

#if __cplusplus >= 201103L
  inline double NCrystal::ncerf(double x) { return std::erf(x); }
#elif defined(__GNUC__) or defined(__clang__)
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

inline void NCrystal::sincos_mpipi(double A,double&cosA, double& sinA)
{
  nc_assert(ncabs(A)<=M_PI);
  cosA = std::cos(A);
#if __cplusplus >= 201103L
  sinA = std::copysign(std::sqrt(1.0-cosA*cosA),A);
#else
  sinA = std::sqrt(1.0-cosA*cosA)*(A<0?-1:1);
#endif
}

inline void NCrystal::sincos_0pi(double A,double&cosA, double& sinA)
{
  nc_assert(0<=A&&A<=M_PI);
  cosA = std::cos(A);
  sinA = std::sqrt(1.0-cosA*cosA);
}

inline void NCrystal::sincos_mpi2pi2(double A,double&cosA, double& sinA)
{
  nc_assert(ncabs(A)<=M_PI_2);
  sinA = std::sin(A);
  cosA = std::sqrt(1.0-sinA*sinA);
}

#endif
