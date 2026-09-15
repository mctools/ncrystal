#ifndef NCrystal_StableDbl_hh
#define NCrystal_StableDbl_hh

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

#include "NCrystal/core/NCDefs.hh"

namespace NCRYSTAL_NAMESPACE {

  // StableDbl uses two doubles internally to keep track of a single value,
  // resulting in less loss of precision during addition, subtraction and
  // multiplication. Like StableSum, this class is intended for use cases where
  // precision is more important than efficiency.

  class StableDbl final {
  public:
    StableDbl() = default;//zero initialises
    StableDbl(double x);
    double value() const;

    //Various operators for +-*:
    StableDbl operator+(const StableDbl&) const;
    StableDbl operator-(const StableDbl&) const;
    StableDbl operator*(const StableDbl&) const;
    StableDbl operator*(double) const;
    StableDbl& operator+=(const StableDbl&);
    StableDbl& operator-=(const StableDbl&);
    StableDbl& operator*=(const StableDbl&);
    StableDbl& operator*=(double);
    StableDbl operator-() const;
    friend StableDbl operator*(double, const StableDbl&);

    //Access internal state:
    static StableDbl fromState(PairDD);
    PairDD state() const;
  private:
    friend std::ostream& operator<<(std::ostream&, const StableDbl&);
    double m_h = 0.0;//leading part
    double m_l = 0.0;//correction
    static void twoSum(double a, double b, double& s, double& e);
    static StableDbl twoProd(double a, double b);
    static StableDbl ddAdd(StableDbl a, StableDbl b);
    StableDbl(double hi, double lo);
  };

  inline std::ostream& operator<<(std::ostream& os, const StableDbl& sd)
  {
    return os << fmt(sd.value());
  }
}

////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  inline StableDbl StableDbl::fromState(PairDD st)
  {
    return { st.first, st.second };
  }

  inline PairDD StableDbl::state() const
  {
    return {m_h,m_l};
  }

  inline StableDbl::StableDbl(double x)
    : m_h(x)
  {
    nc_assert(std::isfinite(x));
  }

  inline StableDbl::StableDbl(double hi, double lo)
    : m_h(hi), m_l(lo)
  {
    nc_assert(std::isfinite(m_h));
    nc_assert(std::isfinite(m_l));
  }

  inline void StableDbl::twoSum(double a, double b, double& s, double& e)
  {
    // Error-free transformation of a floating-point sum.
    // Knuth's six-operation TwoSum algorithm.
    nc_assert(std::isfinite(a));
    nc_assert(std::isfinite(b));
    s = a + b;
    const double z = s - a;
    e = (a - (s - z)) + (b - z);
  }

  inline StableDbl StableDbl::twoProd(double a, double b)
  {
    // FMA-based TwoProd error-free transformation: p + e = a * b
    // Dekker (1971); FMA formulation: Ogita, Rump & Oishi (2005).
    nc_assert(std::isfinite(a));
    nc_assert(std::isfinite(b));
    const double p = a * b;
    return StableDbl(p, std::fma(a, b, -p));
  }

  inline StableDbl StableDbl::ddAdd(StableDbl a, StableDbl b)
  {
    // Double-double addition by expansion addition followed by renormalization.
    // Uses TwoSum error-free transformations to combine the four components.
    // Cf. Joldes, Muller & Popescu (2017), "AccurateDWPlusDW".
    double s, e, t, f, u, v;
    twoSum(a.m_h, b.m_h, s, e);
    twoSum(a.m_l, b.m_l, t, f);
    twoSum(s, e + t, u, v);
    twoSum(u, v + f, s, e);
    return StableDbl(s, e);
  }

  inline double StableDbl::value() const { return m_h + m_l; }

  inline StableDbl StableDbl::operator+(const StableDbl& o) const
  {
    nc_assert(std::isfinite(m_h));
    nc_assert(std::isfinite(m_l));
    nc_assert(std::isfinite(o.m_h));
    nc_assert(std::isfinite(o.m_l));
    return ddAdd(*this, o);
  }

  inline StableDbl StableDbl::operator-(const StableDbl& o) const
  {
    nc_assert(std::isfinite(m_h));
    nc_assert(std::isfinite(m_l));
    nc_assert(std::isfinite(o.m_h));
    nc_assert(std::isfinite(o.m_l));
    return ddAdd(*this, StableDbl(-o.m_h, -o.m_l));
  }

  inline StableDbl StableDbl::operator*(const StableDbl& o) const
  {
    nc_assert(std::isfinite(m_h));
    nc_assert(std::isfinite(m_l));
    nc_assert(std::isfinite(o.m_h));
    nc_assert(std::isfinite(o.m_l));
    StableDbl r = twoProd(m_h, o.m_h);
    r = ddAdd(r, twoProd(m_h, o.m_l));
    r = ddAdd(r, twoProd(m_l, o.m_h));
    r = ddAdd(r, twoProd(m_l, o.m_l));
    return r;
  }

  inline StableDbl StableDbl::operator*(double x) const
  {
    nc_assert(std::isfinite(m_h));
    nc_assert(std::isfinite(m_l));
    nc_assert(std::isfinite(x));
    StableDbl r = twoProd(m_h, x);
    return ddAdd(r, twoProd(m_l, x));
  }

  inline StableDbl& StableDbl::operator+=(const StableDbl& o)
  {
    *this = *this + o;
    return *this;
  }

  inline StableDbl& StableDbl::operator-=(const StableDbl& o)
  {
    *this = *this - o;
    return *this;
  }

  inline StableDbl& StableDbl::operator*=(const StableDbl& o)
  {
    *this = *this * o;
    return *this;
  }

  inline StableDbl& StableDbl::operator*=(double x)
  {
    *this = *this * x;
    return *this;
  }

  inline StableDbl StableDbl::operator-() const
  {
    nc_assert(std::isfinite(m_h));
    nc_assert(std::isfinite(m_l));
    return StableDbl(-m_h, -m_l);
  }

  inline StableDbl operator*(double x, const StableDbl& y)
  {
    return y * x;
  }
}

#endif
