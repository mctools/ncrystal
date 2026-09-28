#ifndef NCrystal_SymOp_hh
#define NCrystal_SymOp_hh

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

#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCrystal/internal/utils/NCVector.hh"

namespace NCRYSTAL_NAMESPACE {

  //Immutable class representing a space group symmetry operation (R,t),
  //acting on fractional coordinates as x' = R*x + t. The rotation part R must
  //be one of the 64 rotations occurring in the 530 space group settings of
  //ITVB: the 48 signed permutation matrices (point group m-3m with
  //orthogonal-type axes), or the 24 rotations of point group 6/mmm with
  //hexagonal axes (8 of which are also signed permutation matrices). These
  //all have entries in {-1,0,1}, and each of the two sets is closed under
  //composition and inversion. The translation t is stored exactly in units
  //of 1/24 (sufficient for all operations of the 530 settings, as well as
  //origin shifts between them), and always reduced modulo 1 (i.e. to
  //0..23/24), since operations are only defined modulo lattice translations.
  //
  //String representations are like "x,-y,z+1/2" or "x-y,x,z+1/6" (as in CIF
  //_space_group_symop_operation_xyz): lower case, no spaces, and translations
  //last. When parsing, spaces, upper case letters, translations before
  //variables ("1/2+x"), any order of terms, and leading "+" signs are also
  //accepted. Translations must be integers or fractions (not decimals) which
  //are exact multiples of 1/24.

  class SymOp final {
  public:
    static constexpr int tdenom = 24;

    //Identity operation:
    constexpr SymOp() noexcept : m_data{ { 1,0,0, 0,1,0, 0,0,1, 0,0,0 } } {}

    //From rotation matrix (row-major) and translation (in units of 1/24, any
    //integers which are reduced modulo 24). Throws BadInput if invalid:
    SymOp( const std::array<int,9>& rot, const std::array<int,3>& trans );

    //From strings like "x,-y,z+1/2". Throws BadInput if invalid:
    explicit SymOp( StrView );
    explicit SymOp( const std::string& s ) : SymOp( StrView( s ) ) {}
    explicit SymOp( const char * s ) : SymOp( StrView( s ) ) {}

    int rot( unsigned i, unsigned j ) const ncnoexceptndebug;//row i, column j
    int trans( unsigned i ) const ncnoexceptndebug;//in units of 1/24 (0..23)
    bool isIdentity() const noexcept;
    bool hasTranslation() const noexcept;
    int determinant() const noexcept;//+1 or -1
    int rotationType() const ncnoexceptndebug;//1,2,3,4,6,-1,-2(=m),-3,-4,-6
    unsigned rotationOrder() const ncnoexceptndebug;//smallest n>0 with R^n=I

    //Composition, a*b meaning "b followed by a", and inverse. Composing
    //operations with hexagonal-only and orthogonal-only rotation parts (which
    //can not belong to the same space group) would not give a valid operation,
    //and throws LogicError:
    SymOp operator*( const SymOp& ) const;
    SymOp inverse() const ncnoexceptndebug;

    //Whether a matrix (row-major) is one of the 64 allowed rotation parts:
    static bool isAllowedRotation( const std::array<int,9>& ) ncnoexceptndebug;

    //Apply to fractional coordinates (no wrapping into the unit cell):
    Vector apply( const Vector& ) const noexcept;

    std::string toString() const;//e.g. "x-y,x,z+1/6"

    bool operator==( const SymOp& o ) const noexcept;
    bool operator!=( const SymOp& o ) const noexcept;
    bool operator<( const SymOp& o ) const noexcept;//arbitrary but fixed

  private:
    //Rotation (row-major), followed by translation:
    std::array<std::int8_t,12> m_data;
    struct NoCheck {};
    SymOp( NoCheck, const std::array<std::int8_t,12>& d ) noexcept
      : m_data( d ) {}
    static std::int8_t reduceTrans( int ) noexcept;
    [[noreturn]] static void throwInvalidProduct( const SymOp&, const SymOp& );
  };

  std::ostream& operator<<( std::ostream&, const SymOp& );

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  inline int SymOp::rot( unsigned i, unsigned j ) const ncnoexceptndebug
  {
    nc_assert( i < 3 && j < 3 );
    return m_data[ 3 * i + j ];
  }
  inline int SymOp::trans( unsigned i ) const ncnoexceptndebug
  {
    nc_assert( i < 3 );
    return m_data[ 9 + i ];
  }
  inline bool SymOp::operator==( const SymOp& o ) const noexcept
  {
    return m_data == o.m_data;
  }
  inline bool SymOp::operator!=( const SymOp& o ) const noexcept
  {
    return m_data != o.m_data;
  }
  inline bool SymOp::operator<( const SymOp& o ) const noexcept
  {
    return m_data < o.m_data;
  }
  inline std::int8_t SymOp::reduceTrans( int t ) noexcept
  {
    return static_cast<std::int8_t>( ( ( t % tdenom ) + tdenom ) % tdenom );
  }
  inline SymOp SymOp::operator*( const SymOp& o ) const
  {
    //(R1,t1)*(R2,t2) = (R1*R2, R1*t2+t1):
    std::array<int,9> r;
    std::array<std::int8_t,12> d;
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j )
        r[3*i+j] = ( m_data[3*i] * o.m_data[j]
                     + m_data[3*i+1] * o.m_data[3+j]
                     + m_data[3*i+2] * o.m_data[6+j] );
      d[9+i] = reduceTrans( m_data[3*i] * o.m_data[9]
                            + m_data[3*i+1] * o.m_data[10]
                            + m_data[3*i+2] * o.m_data[11]
                            + m_data[9+i] );
    }
    if ( !isAllowedRotation( r ) )
      throwInvalidProduct( *this, o );
    for ( int i = 0; i < 9; ++i )
      d[i] = static_cast<std::int8_t>( r[i] );
    return SymOp( NoCheck(), d );
  }
  inline Vector SymOp::apply( const Vector& v ) const noexcept
  {
    constexpr double f = 1.0 / tdenom;
    return Vector( m_data[0]*v.x() + m_data[1]*v.y() + m_data[2]*v.z()
                   + m_data[9] * f,
                   m_data[3]*v.x() + m_data[4]*v.y() + m_data[5]*v.z()
                   + m_data[10] * f,
                   m_data[6]*v.x() + m_data[7]*v.y() + m_data[8]*v.z()
                   + m_data[11] * f );
  }
}

#endif
