#ifndef NCrystal_Vector_hh
#define NCrystal_Vector_hh

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

#include "NCrystal/core/NCDefs.hh"

namespace NCRYSTAL_NAMESPACE {

  class Vector final : public StronglyTypedFixedVector<Vector,double,3> {
  public:

    // Simple Vector class which adds various mathematical operations on top of
    // the StronglyTypedFixedVector functionality.

    using StronglyTypedFixedVector::StronglyTypedFixedVector;

    ncconstexpr17 Vector operator*(const Vector&) const;
    ncconstexpr17 Vector operator* (double ) const;
    ncconstexpr17 Vector operator/ (double ) const;
    ncconstexpr17 Vector& operator*= (double );
    ncconstexpr17 Vector& operator+= (const Vector&);
    ncconstexpr17 Vector& operator-= (const Vector&);
    ncconstexpr17 Vector& operator/= (double );

    ncconstexpr17 Vector operator-() const;
    ncconstexpr17 Vector operator+( const Vector&) const;
    ncconstexpr17 Vector operator-( const Vector&) const;

    Vector unit() const;//slow
    void normalise();//better
    ncconstexpr17 Vector cross(const Vector&) const;
    ncconstexpr17 void cross_inplace(const Vector&);
    constexpr double dot(const Vector&) const;
    double angle(const Vector&) const;//slow
    double angle_highres(const Vector&) const;//very slow, but precise even for small angles
    double mag() const;//slow
    constexpr double mag2() const;//better
    void setMag( double );//slow
    ncconstexpr17 bool isParallel(const Vector&, double epsilon = 1e-10) const;
    ncconstexpr17 bool isOrthogonal(const Vector&, double epsilon = 1e-10) const;
    bool isUnitVector(double tolerance = 1e-10) const;
    constexpr bool isStrictNullVector() const noexcept;

    constexpr const double& x() const noexcept { return m_data[0]; }
    constexpr const double& y() const noexcept { return m_data[1]; }
    constexpr const double& z() const noexcept { return m_data[2]; }
    ncconstexpr17 double& x() noexcept { return m_data[0]; }
    ncconstexpr17 double& y() noexcept { return m_data[1]; }
    ncconstexpr17 double& z() noexcept { return m_data[2]; }

    //Sort x, then y, the z (assumes no NaN):
    constexpr int lexCmp( const Vector& ) const noexcept;

  };

}


////////////////////////////
// Inline implementations //
////////////////////////////

inline double NCrystal::Vector::mag() const
{
  return std::sqrt( mag2() );
}

inline constexpr double NCrystal::Vector::mag2() const
{
  return m_data[0]*m_data[0]+m_data[1]*m_data[1]+m_data[2]*m_data[2];
}

inline constexpr double NCrystal::Vector::dot(const NCrystal::Vector& o) const
{
  return m_data[0]*o.m_data[0]+m_data[1]*o.m_data[1]+m_data[2]*o.m_data[2];
}

inline ncconstexpr17 NCrystal::Vector& NCrystal::Vector::operator+=( const NCrystal::Vector& o)
{
  m_data[0] += o.m_data[0];
  m_data[1] += o.m_data[1];
  m_data[2] += o.m_data[2];
  return *this;
}

inline ncconstexpr17 NCrystal::Vector& NCrystal::Vector::operator-=( const NCrystal::Vector& o)
{
  m_data[0] -= o.m_data[0];
  m_data[1] -= o.m_data[1];
  m_data[2] -= o.m_data[2];
  return *this;
}

inline ncconstexpr17 NCrystal::Vector& NCrystal::Vector::operator*= (double f)
{
  m_data[0] *= f;
  m_data[1] *= f;
  m_data[2] *= f;
  return *this;
}

inline ncconstexpr17 NCrystal::Vector& NCrystal::Vector::operator/= (double f)
{
  const double ff(1.0/f);
  m_data[0] *= ff;
  m_data[1] *= ff;
  m_data[2] *= ff;
  return *this;
}

inline ncconstexpr17 NCrystal::Vector NCrystal::Vector::operator/ (double f) const
{
  const double ff(1.0/f);
  return Vector( m_data[0]*ff, m_data[1]*ff, m_data[2]*ff );
}

inline ncconstexpr17 NCrystal::Vector NCrystal::Vector::operator* (double f) const
{
  return Vector( m_data[0]*f, m_data[1]*f, m_data[2]*f );
}

inline ncconstexpr17 NCrystal::Vector NCrystal::Vector::operator*(const NCrystal::Vector& o) const
{
  return Vector( m_data[1]*o.m_data[2] - m_data[2]*o.m_data[1],
                 m_data[2]*o.m_data[0] - m_data[0]*o.m_data[2],
                 m_data[0]*o.m_data[1] - m_data[1]*o.m_data[0] );
}

inline ncconstexpr17 void NCrystal::Vector::cross_inplace(const Vector&o)
{
  double xx = m_data[1]*o.m_data[2] - m_data[2]*o.m_data[1];
  double yy = m_data[2]*o.m_data[0] - m_data[0]*o.m_data[2];
  m_data[2] = m_data[0]*o.m_data[1] - m_data[1]*o.m_data[0];
  m_data[0] = xx;
  m_data[1] = yy;
}

inline ncconstexpr17 NCrystal::Vector NCrystal::Vector::operator-() const
{
  return Vector( -m_data[0], -m_data[1], -m_data[2] );
}

inline ncconstexpr17 NCrystal::Vector NCrystal::Vector::operator+( const NCrystal::Vector& o ) const
{
  return Vector( m_data[0]+o.m_data[0], m_data[1]+o.m_data[1], m_data[2]+o.m_data[2] );
}

inline ncconstexpr17 NCrystal::Vector NCrystal::Vector::operator-( const NCrystal::Vector& o ) const
{
  return Vector( m_data[0]-o.m_data[0], m_data[1]-o.m_data[1], m_data[2]-o.m_data[2] );
}

inline NCrystal::Vector NCrystal::Vector::unit() const
{
  double themag2 = mag2();
  if (themag2==1.0)
    return *this;
  if (!themag2)
    NCRYSTAL_THROW(CalcError,"NCVector::unit(): Can't scale null-vector.");
  double factor = 1.0 / std::sqrt(themag2);
  return { m_data[0]*factor, m_data[1]*factor, m_data[2]*factor };
}

inline constexpr int NCrystal::Vector::lexCmp( const Vector& o ) const noexcept
{
  return ( x() != o.x()
           ? ( x() < o.x() ? -1 : 1 )
           : ( y() != o.y()
               ? ( y() < o.y() ? -1 : 1 )
               : ( z() != o.z()
                   ? ( z() < o.z() ? -1 : 1 )
                   : 0
                   )
               )
           );
}

inline void NCrystal::Vector::normalise()
{
  double themag2 = mag2();
  constexpr double one_low = 1.0 - 2.0 * std::numeric_limits<double>::epsilon();
  constexpr double one_high = 1.0 + 2.0 * std::numeric_limits<double>::epsilon();
  if ( themag2 >= one_low && themag2 <= one_high )
    return;//already normalised (normalising it again might change the values slightly)
  if (!themag2)
    NCRYSTAL_THROW(CalcError,"NCVector::normalise(): Can't scale null-vector.");
  if (std::isinf(themag2))
    NCRYSTAL_THROW(CalcError,"NCVector::normalise(): Can't scale vector with infinite length.");
  *this /= std::sqrt(themag2);
}

inline ncconstexpr17 NCrystal::Vector NCrystal::Vector::cross(const NCrystal::Vector& o) const
{
  return *this * o;
}

inline ncconstexpr17 bool NCrystal::Vector::isParallel(const NCrystal::Vector& vec2, double epsilon) const
{
  //NB: using '>' rather than '>=' to have null-vectors never be parallel to
  //anything (including themselves, which we could of course check against).
  double dp = dot(vec2);
  return dp*dp > mag2() * vec2.mag2() * ( 1.0 - epsilon);
}

inline ncconstexpr17 bool NCrystal::Vector::isOrthogonal(const Vector& vec2, double epsilon) const
{
  //NB: using '<' rather than '<=' to have null-vectors never be orthogonal to
  //anything.
  double dp = dot(vec2);
  return dp*dp < mag2() * vec2.mag2() * epsilon;
}

inline double NCrystal::Vector::angle(const NCrystal::Vector& vec2) const
{
  double norm = std::sqrt( mag2()*vec2.mag2() );
  if (!norm)
    NCRYSTAL_THROW(CalcError,"NCVector::angle(): Can't find angle to/from null-vector.");
  double result = dot(vec2) / norm;
  return std::acos( std::min(1.,std::max(-1.,result)) );
}

inline double NCrystal::Vector::angle_highres(const NCrystal::Vector& vec2) const
{
  //Based on formula on page 47 of
  //https://people.eecs.berkeley.edu/~wkahan/Mindless.pdf "How Futile are
  //Mindless Assessments of Roundoff in Floating-Point Computation?" by W. Kahan
  //(Jan 11, 2006):
  NCrystal::Vector a(*this);
  NCrystal::Vector b(vec2);
  double mag2_a = a.mag2();
  double mag2_b = b.mag2();
  if (!mag2_a||!mag2_b)
    NCRYSTAL_THROW(CalcError,"NCVector::angle_highres(): Can't find angle to/from null-vector.");
  a *= 1.0/std::sqrt(mag2_a);
  b *= 1.0/std::sqrt(mag2_b);
  return 2*std::atan2((a-b).mag(),(a+b).mag());
}

inline constexpr bool NCrystal::Vector::isStrictNullVector() const noexcept
{
  return m_data[0]==0.0 && m_data[1]==0.0 && m_data[2]==0.0;
}

inline bool NCrystal::Vector::isUnitVector(double tolerance) const
{
  return std::abs( mag2() - 1.0 ) < tolerance;
}

#endif
