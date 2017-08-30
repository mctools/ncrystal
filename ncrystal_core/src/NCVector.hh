#ifndef NCrystal_Vector_hh
#define NCrystal_Vector_hh

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
#include "NCMath.hh"
#include <ostream>

//Simple Vector class

namespace NCrystal {

  class Vector
  {
  public:

    Vector();
    Vector(double x, double y, double z);
    ~Vector(){};

    Vector& operator=( const Vector&);
    Vector operator*(const Vector&) const;
    Vector operator* (double ) const;
    Vector operator/ (double ) const;
    Vector& operator*= (double );
    Vector& operator+= (const Vector&);
    Vector& operator-= (const Vector&);
    Vector& operator/= (double );

    Vector operator-() const;
    Vector operator+( const Vector&) const;
    Vector operator-( const Vector&) const;

    friend std::ostream& operator << (std::ostream &, const Vector&);

    bool operator==( const Vector&) const;
    bool operator!=( const Vector&) const;

    void print() const;
    Vector unit() const;//slow
    void normalise();//better
    Vector cross(const Vector&) const;
    void cross_inplace(const Vector&);
    double dot(const Vector&) const;
    double angle(const Vector&) const;//slow
    double mag() const;//slow
    double mag2() const;//better
    void set(double x, double y, double z);
    void setMag(double );//slow
    bool isParallel(const Vector&, double epsilon = 1e-10) const;

    inline const double& x() const {return m_x;}
    inline const double& y() const {return m_y;}
    inline const double& z() const {return m_z;}

  protected:
    //Keep data members exactly like this, so Vector objects can reliably be
    //reinterpreted as double[3] arrays and vice versa:
    double m_x;
    double m_y;
    double m_z;
  };

  //For interpreting double[3] arrays as Vector:
  static inline Vector& asVect( double (&v)[3] ) { return *reinterpret_cast<Vector*>(&v); }
  static inline const Vector& asVect( const double (&v)[3] ) { return *reinterpret_cast<const Vector*>(&v); }

  std::ostream& operator << (std::ostream &, const Vector&);

}


////////////////////////////
// Inline implementations //
////////////////////////////

inline NCrystal::Vector::Vector()
  : m_x(0.), m_y(0.), m_z(0.)
{
}

inline NCrystal::Vector::Vector(double vx, double vy, double vz)
  : m_x(vx), m_y(vy), m_z(vz)
{
}

inline void NCrystal::Vector::set(double xx, double yy, double zz)
{
  m_x = xx;
  m_y = yy;
  m_z = zz;
}

inline double NCrystal::Vector::mag() const
{
  return std::sqrt( m_x*m_x + m_y*m_y + m_z*m_z );
}

inline double NCrystal::Vector::mag2() const
{
  return m_x*m_x + m_y*m_y + m_z*m_z;
}

inline double NCrystal::Vector::dot(const NCrystal::Vector& o) const
{
  return m_x*o.m_x + m_y*o.m_y + m_z*o.m_z;
}

inline NCrystal::Vector& NCrystal::Vector::operator=( const NCrystal::Vector& o)
{
  m_x = o.m_x;
  m_y = o.m_y;
  m_z = o.m_z;
  return *this;
}

inline NCrystal::Vector& NCrystal::Vector::operator+=( const NCrystal::Vector& o)
{
  m_x += o.m_x;
  m_y += o.m_y;
  m_z += o.m_z;
  return *this;
}

inline NCrystal::Vector& NCrystal::Vector::operator-=( const NCrystal::Vector& o)
{
  m_x -= o.m_x;
  m_y -= o.m_y;
  m_z -= o.m_z;
  return *this;
}

inline NCrystal::Vector& NCrystal::Vector::operator*= (double f)
{
  m_x *= f;
  m_y *= f;
  m_z *= f;
  return *this;
}

inline NCrystal::Vector& NCrystal::Vector::operator/= (double f)
{
  double ff(1.0/f);
  m_x *= ff;
  m_y *= ff;
  m_z *= ff;
  return *this;
}

inline bool NCrystal::Vector::operator==( const NCrystal::Vector& o) const
{
  return ( m_x==o.m_x && m_y==o.m_y && m_z==o.m_z );
}

inline bool NCrystal::Vector::operator!=( const NCrystal::Vector& o) const
{
  return !( (*this) == o );
}

inline NCrystal::Vector NCrystal::Vector::operator/ (double f) const
{
  return Vector( m_x/f, m_y/f, m_z/f );
}

inline NCrystal::Vector NCrystal::Vector::operator* (double f) const
{
  return Vector( m_x*f, m_y*f, m_z*f );
}

inline NCrystal::Vector NCrystal::Vector::operator*(const NCrystal::Vector& o) const
{
  return Vector( m_y*o.m_z - m_z*o.m_y,
                 m_z*o.m_x - m_x*o.m_z,
                 m_x*o.m_y - m_y*o.m_x );
}

inline void NCrystal::Vector::cross_inplace(const Vector&o)
{
  double xx = m_y*o.m_z - m_z*o.m_y;
  double yy = m_z*o.m_x - m_x*o.m_z;
  m_z = m_x*o.m_y - m_y*o.m_x;
  m_x = xx;
  m_y = yy;
}

inline NCrystal::Vector NCrystal::Vector::operator-() const
{
  return NCrystal::Vector( -m_x, -m_y, -m_z );
}

inline NCrystal::Vector NCrystal::Vector::operator+( const NCrystal::Vector& o ) const
{
  return NCrystal::Vector( m_x+o.m_x, m_y+o.m_y, m_z+o.m_z );
}

inline NCrystal::Vector NCrystal::Vector::operator-( const NCrystal::Vector& o ) const
{
  return NCrystal::Vector( m_x-o.m_x, m_y-o.m_y, m_z-o.m_z );
}

inline NCrystal::Vector NCrystal::Vector::unit() const
{
  double themag2 = mag2();
  if (themag2==1.0)
    return *this;
  if (!themag2)
    NCRYSTAL_THROW(CalcError,"NCVector::unit(): Can't scale null-vector.");
  double factor = 1.0/std::sqrt(themag2);
  return NCrystal::Vector(m_x*factor, m_y*factor, m_z*factor);
}

inline void NCrystal::Vector::normalise()
{
  double themag2 = mag2();
  if (themag2==1.0)
    return;
  if (!themag2)
    NCRYSTAL_THROW(CalcError,"NCVector::normalise(): Can't scale null-vector.");
  double f = 1.0/std::sqrt(themag2);
  m_x *= f;
  m_y *= f;
  m_z *= f;
}

inline NCrystal::Vector NCrystal::Vector::cross(const NCrystal::Vector& o) const
{
  return *this * o;
}

inline bool NCrystal::Vector::isParallel(const NCrystal::Vector& vec2, double epsilon) const
{
  //NB: using '>' rather than '>=' to have null-vectors never be parallel to
  //anything (including themselves, which we could of course check against).
  double dp = dot(vec2);
  return dp*dp > mag2() * vec2.mag2() * ( 1.0 - epsilon);
}

inline double NCrystal::Vector::angle(const NCrystal::Vector& vec2) const
{
  double result = dot(vec2) / std::sqrt( mag2()*vec2.mag2() );
  return std::acos( result );
}

#endif
