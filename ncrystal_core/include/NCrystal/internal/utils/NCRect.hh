#ifndef NCrystal_Rect_hh
#define NCrystal_Rect_hh

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
#include "NCrystal/internal/utils/NCMath.hh"

namespace NCRYSTAL_NAMESPACE {

  class Rectangle final {
  public:
    // Immutable class representing an axis-aligned rectangle, i.e. the set of
    // (x,y) values with x0<=x<=x1 and y0<=y<=y1.
    // It is not allowed to access coordinates of empty rectangles.
    Rectangle();//Empty
    Rectangle(const PairDD& xRange, const PairDD& yRange);
    Rectangle(double x0, double x1, double y0, double y1);

    bool operator==(const Rectangle& r) const;
    bool operator<(const Rectangle& r) const;

    double x0() const;
    double x1() const;
    double y0() const;
    double y1() const;

    const PairDD& xRange() const;
    const PairDD& yRange() const;

    bool isEmpty() const;

    Rectangle getUnion(const Rectangle& r) const;
    Rectangle getIntersection(const Rectangle& r) const;

    bool ptIsInside(double x, double y) const;
    bool ptIsInside(const PairDD&) const;

  private:
    friend std::ostream& operator<<(std::ostream&, const Rectangle&);
    PairDD m_xRange;
    PairDD m_yRange;
  };

  std::ostream& operator<<(std::ostream&, const Rectangle&);
}


////////////////////////////
// Inline implementations //
////////////////////////////


namespace NCRYSTAL_NAMESPACE {
  inline Rectangle::Rectangle()
    : m_xRange(0.0, 0.0), m_yRange(0.0, 0.0)
  {
  }

  inline Rectangle::Rectangle(const PairDD& xRange, const PairDD& yRange)
    : m_xRange(xRange), m_yRange(yRange)
  {
    nc_assert(std::isfinite(m_xRange.first));
    nc_assert(std::isfinite(m_xRange.second));
    nc_assert(std::isfinite(m_yRange.first));
    nc_assert(std::isfinite(m_yRange.second));
    nc_assert(m_xRange.first < m_xRange.second);
    nc_assert(m_yRange.first < m_yRange.second);
  }

  inline Rectangle::Rectangle(double x0, double x1, double y0, double y1)
    : Rectangle(PairDD(x0, x1), PairDD(y0, y1))
  {
  }

  inline double Rectangle::x0() const
  {
    nc_assert(!isEmpty());
    return m_xRange.first;
  }

  inline double Rectangle::x1() const
  {
    nc_assert(!isEmpty());
    return m_xRange.second;
  }

  inline double Rectangle::y0() const
  {
    nc_assert(!isEmpty());
    return m_yRange.first;
  }

  inline double Rectangle::y1() const
  {
    nc_assert(!isEmpty());
    return m_yRange.second;
  }

  inline const PairDD& Rectangle::xRange() const
  {
    nc_assert(!isEmpty());
    return m_xRange;
  }

  inline const PairDD& Rectangle::yRange() const
  {
    nc_assert(!isEmpty());
    return m_yRange;
  }

  inline bool Rectangle::isEmpty() const
  {
    //All empty rectangles have data members: ((0,0),(0,0))
    nc_assert( bool(!(m_xRange.first < m_xRange.second))
               == bool( m_xRange.first == 0.0
                        && m_xRange.second == 0.0
                        && m_yRange.first == 0.0
                        && m_yRange.second == 0.0 ) );
    return !(m_xRange.first < m_xRange.second);
  }

  inline Rectangle Rectangle::getUnion(const Rectangle& r) const
  {
    nc_assert(std::isfinite(r.m_xRange.first));
    nc_assert(std::isfinite(r.m_xRange.second));
    nc_assert(std::isfinite(r.m_yRange.first));
    nc_assert(std::isfinite(r.m_yRange.second));
    if (isEmpty())
      return r;
    if (r.isEmpty())
      return *this;
    return Rectangle(PairDD(ncmin(m_xRange.first, r.m_xRange.first),
                            ncmax(m_xRange.second, r.m_xRange.second)),
                     PairDD(ncmin(m_yRange.first, r.m_yRange.first),
                            ncmax(m_yRange.second, r.m_yRange.second)));
  }

  inline Rectangle Rectangle::getIntersection(const Rectangle& r) const
  {
    nc_assert(std::isfinite(r.m_xRange.first));
    nc_assert(std::isfinite(r.m_xRange.second));
    nc_assert(std::isfinite(r.m_yRange.first));
    nc_assert(std::isfinite(r.m_yRange.second));
    if (isEmpty() || r.isEmpty())
      return Rectangle();
    const PairDD xRange(ncmax(m_xRange.first, r.m_xRange.first),
                        ncmin(m_xRange.second, r.m_xRange.second));
    const PairDD yRange(ncmax(m_yRange.first, r.m_yRange.first),
                        ncmin(m_yRange.second, r.m_yRange.second));
    if (xRange.first >= xRange.second || yRange.first >= yRange.second)
      return Rectangle();
    return Rectangle(xRange, yRange);
  }

  inline bool Rectangle::ptIsInside(const PairDD& xy) const
  {
    return ptIsInside( xy.first, xy.second );
  }

  inline bool Rectangle::ptIsInside(double x, double y) const
  {
    nc_assert(std::isfinite(x));
    nc_assert(std::isfinite(y));
    return ( !isEmpty() &&
             m_xRange.first <= x && x <= m_xRange.second &&
             m_yRange.first <= y && y <= m_yRange.second );
  }

  inline bool Rectangle::operator<(const Rectangle& r) const
  {
    if (m_xRange < r.m_xRange)
      return true;
    if (r.m_xRange < m_xRange)
      return false;
    return m_yRange < r.m_yRange;
  }

  inline bool Rectangle::operator==(const Rectangle& r) const
  {
    return m_xRange == r.m_xRange && m_yRange == r.m_yRange;
  }
}

#endif
