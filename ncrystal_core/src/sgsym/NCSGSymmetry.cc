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

#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"
#include <set>
#include <mutex>

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {

  namespace detail {
    //Raw generation of all operations (incl. centring, in canonical order)
    //from a Hall symbol. Not intended for general usage (use SpaceGroup
    //instead), but made available for testing purposes. Throws BadInput in
    //case of invalid symbols:
    std::vector<SymOp> rawSymOpsFromHallSymbol( StrView );
  }

  namespace {

    using Mat = std::array<int,9>;//row-major

    [[noreturn]] void throwBadHall( StrView h, const std::string& reason )
    {
      NCRYSTAL_THROW2(BadInput,"Invalid Hall symbol \""<<h<<"\": "<<reason);
    }

    constexpr Mat identityMat() { return { { 1,0,0, 0,1,0, 0,0,1 } }; }

    //Proper rotation matrices of Hall symbols (cf. ITVB Table A1.4.2.4), for
    //rotation order n along an axis ('x','y','z', face diagonals '\'' and
    //'"' which depend on the preceding axis, or body diagonal '*'). Returns
    //false if the combination is not valid:
    bool hallRotation( int n, char axis, char preceding_axis, Mat& m )
    {
      if ( n == 1 ) {
        m = identityMat();
        return true;
      }
      if ( axis == 'x' || axis == 'y' || axis == 'z' ) {
        static const Mat tbl[3][4] = {
          //x: 2, 3, 4, 6
          { { { 1,0,0, 0,-1,0, 0,0,-1 } }, { { 1,0,0, 0,0,-1, 0,1,-1 } },
            { { 1,0,0, 0,0,-1, 0,1,0 } }, { { 1,0,0, 0,1,-1, 0,1,0 } } },
          //y: 2, 3, 4, 6
          { { { -1,0,0, 0,1,0, 0,0,-1 } }, { { -1,0,1, 0,1,0, -1,0,0 } },
            { { 0,0,1, 0,1,0, -1,0,0 } }, { { 0,0,1, 0,1,0, -1,0,1 } } },
          //z: 2, 3, 4, 6
          { { { -1,0,0, 0,-1,0, 0,0,1 } }, { { 0,-1,0, 1,-1,0, 0,0,1 } },
            { { 0,-1,0, 1,0,0, 0,0,1 } }, { { 1,-1,0, 1,0,0, 0,0,1 } } } };
        const int ni = ( n == 2 ? 0 : n == 3 ? 1 : n == 4 ? 2 : n == 6 ? 3 : -1 );
        if ( ni < 0 )
          return false;
        m = tbl[axis - 'x'][ni];
        return true;
      }
      if ( axis == '\'' && preceding_axis == '*' )
        preceding_axis = 'z';//a-b axis (perpendicular to a+b+c and c)
      if ( axis == '\'' || axis == '"' ) {
        if ( n != 2 || !( preceding_axis == 'x' || preceding_axis == 'y'
                          || preceding_axis == 'z' ) )
          return false;
        static const Mat tbl[2][3] = {
          //' after x, y, z:
          { { { -1,0,0, 0,0,-1, 0,-1,0 } }, { { 0,0,-1, 0,-1,0, -1,0,0 } },
            { { 0,-1,0, -1,0,0, 0,0,-1 } } },
          //" after x, y, z:
          { { { -1,0,0, 0,0,1, 0,1,0 } }, { { 0,0,1, 0,-1,0, 1,0,0 } },
            { { 0,1,0, 1,0,0, 0,0,-1 } } } };
        m = tbl[ axis == '"' ? 1 : 0 ][ preceding_axis - 'x' ];
        return true;
      }
      if ( axis == '*' && n == 3 ) {
        m = { { 0,0,1, 1,0,0, 0,1,0 } };
        return true;
      }
      return false;
    }

    SymOp makeOp( const Mat& m, const std::array<int,3>& t )
    {
      return SymOp( m, t );
    }

    struct CanonicalOps {
      std::vector<SymOp> ops, reps, centring;
    };

    //Canonical ordering of a complete group of operations (see
    //NCSGSymmetry.hh):
    CanonicalOps canonicalOrder( const std::set<SymOp>& group )
    {
      CanonicalOps res;
      std::map<Mat,SymOp> rot2rep;//(std::set is sorted, so first is smallest)
      for ( auto& op : group ) {
        Mat m;
        for ( auto i : ncrange( 9u ) )
          m[i] = op.rot( i / 3, i % 3 );
        if ( m == identityMat() )
          res.centring.push_back( op );
        rot2rep.emplace( m, op );//no-op if already present
      }
      for ( auto& e : rot2rep )
        res.reps.push_back( e.second );
      std::sort( res.reps.begin(), res.reps.end() );
      std::sort( res.centring.begin(), res.centring.end() );
      //Identity first (the representative with identity rotation is the
      //identity itself, since that has the smallest translation):
      auto it = std::find( res.reps.begin(), res.reps.end(), SymOp() );
      nc_assert_always( it != res.reps.end() );
      std::rotate( res.reps.begin(), it, it + 1 );
      nc_assert_always( !res.centring.empty()
                        && res.centring.front() == SymOp() );
      for ( auto& c : res.centring )
        for ( auto& r : res.reps )
          res.ops.push_back( c * r );
      nc_assert_always( res.ops.size() == group.size() );
      nc_assert_always( std::set<SymOp>( res.ops.begin(), res.ops.end() )
                        == group );
      return res;
    }

    //Parse Hall symbol and generate complete group:
    std::set<SymOp> groupFromHallSymbol( StrView hall, char* lattice_symbol )
    {
      //Separate the change-of-basis vector "(vx vy vz)" (units of 1/12):
      std::string main = hall.to_string();
      Optional<std::array<int,3>> shift;
      const auto ipar = main.find( '(' );
      if ( ipar != std::string::npos ) {
        const auto iend = main.find( ')', ipar );
        if ( iend == std::string::npos || iend + 1 != main.size() )
          throwBadHall( hall, "invalid change-of-basis vector" );
        const auto vparts = StrView( main ).substr( ipar + 1, iend - ipar - 1 )
          .split();
        if ( vparts.size() != 3 )
          throwBadHall( hall, "change-of-basis vector must have three"
                        " components" );
        std::array<int,3> v;
        for ( auto i : ncrange( 3u ) ) {
          auto iv = vparts.at( i ).toInt32();
          if ( !iv.has_value() || iv.value() < -24 || iv.value() > 24 )
            throwBadHall( hall, "invalid change-of-basis vector component" );
          v[i] = 2 * iv.value();//1/12 -> 1/24 units
        }
        shift = v;
        main.resize( ipar );
      }
      const auto tokens = StrView( main ).split();
      if ( tokens.empty() )
        throwBadHall( hall, "empty symbol" );

      //Lattice symbol:
      std::vector<SymOp> generators;
      StrView lat = tokens.front();
      if ( lat.startswith( '-' ) ) {
        generators.push_back( makeOp( { { -1,0,0, 0,-1,0, 0,0,-1 } },
                                      { { 0,0,0 } } ) );
        lat = lat.substr( 1 );
      }
      if ( lat.size() != 1 )
        throwBadHall( hall, "invalid lattice symbol" );
      switch ( lat[0] ) {
      case 'P': break;
      case 'A': generators.push_back( makeOp( identityMat(), { { 0,12,12 } } ) );
        break;
      case 'B': generators.push_back( makeOp( identityMat(), { { 12,0,12 } } ) );
        break;
      case 'C': generators.push_back( makeOp( identityMat(), { { 12,12,0 } } ) );
        break;
      case 'I': generators.push_back( makeOp( identityMat(), { { 12,12,12 } } ) );
        break;
      case 'R': generators.push_back( makeOp( identityMat(), { { 16,8,8 } } ) );
        break;
      case 'F':
        generators.push_back( makeOp( identityMat(), { { 0,12,12 } } ) );
        generators.push_back( makeOp( identityMat(), { { 12,0,12 } } ) );
        break;
      default:
        throwBadHall( hall, "invalid lattice symbol" );
      }
      *lattice_symbol = lat[0];

      //Rotation symbols:
      int prev_n = 0;
      char prev_axis = 0;
      for ( auto itok : ncrange( std::size_t(1), tokens.size() ) ) {
        const StrView tok = tokens.at( itok );
        const unsigned irot = static_cast<unsigned>( itok - 1 );
        std::size_t i = 0;
        const bool improper = tok.startswith( '-' );
        if ( improper )
          ++i;
        if ( !( i < tok.size() && ( tok[i] == '1' || tok[i] == '2'
                                    || tok[i] == '3' || tok[i] == '4'
                                    || tok[i] == '6' ) ) )
          throwBadHall( hall, "invalid rotation symbol \"" + tok.to_string()
                        + "\"" );
        const int n = tok[i++] - '0';
        //Screw translation digit (e.g. "31", "61"):
        int screw = 0;
        if ( i < tok.size() && tok[i] >= '1' && tok[i] <= '5' ) {
          screw = tok[i++] - '0';
          if ( !( n >= 3 && screw < n ) )
            throwBadHall( hall, "invalid screw translation in \""
                          + tok.to_string() + "\"" );
        }
        //Axis:
        char axis = 0;
        if ( i < tok.size() && ( tok[i] == 'x' || tok[i] == 'y'
                                 || tok[i] == 'z' || tok[i] == '\''
                                 || tok[i] == '"' || tok[i] == '*' ) )
          axis = tok[i++];
        if ( !axis && n != 1 ) {
          //Default axis rules:
          if ( irot == 0 )
            axis = 'z';
          else if ( irot == 1 && n == 2 && ( prev_n == 2 || prev_n == 4 ) )
            axis = 'x';
          else if ( irot == 1 && n == 2 && ( prev_n == 3 || prev_n == 6 ) )
            axis = '\'';
          else if ( irot == 2 && n == 3 )
            axis = '*';
          else
            throwBadHall( hall, "can not determine axis of \""
                          + tok.to_string() + "\"" );
        }
        //Translation letters:
        std::array<int,3> t = { { 0, 0, 0 } };
        for ( ; i < tok.size(); ++i ) {
          switch ( tok[i] ) {
          case 'a': t[0] += 12; break;
          case 'b': t[1] += 12; break;
          case 'c': t[2] += 12; break;
          case 'n': t[0] += 12; t[1] += 12; t[2] += 12; break;
          case 'u': t[0] += 6; break;
          case 'v': t[1] += 6; break;
          case 'w': t[2] += 6; break;
          case 'd': t[0] += 6; t[1] += 6; t[2] += 6; break;
          default:
            throwBadHall( hall, "invalid rotation symbol \"" + tok.to_string()
                          + "\"" );
          }
        }
        if ( screw ) {
          if ( !( axis == 'x' || axis == 'y' || axis == 'z' ) )
            throwBadHall( hall, "screw translation requires principal axis" );
          t[axis - 'x'] += 24 * screw / n;
        }
        Mat m;
        if ( !hallRotation( n, axis, prev_axis, m ) )
          throwBadHall( hall, "invalid rotation symbol \"" + tok.to_string()
                        + "\"" );
        if ( improper )
          for ( auto& e : m )
            e = -e;
        try {
          generators.push_back( makeOp( m, t ) );
        } catch ( Error::BadInput& e ) {
          throwBadHall( hall, e.what() );
        }
        if ( n != 1 ) {
          prev_n = n;
          prev_axis = axis;
        }
      }

      //Closure (right-multiplying by generators, starting from identity):
      std::set<SymOp> group{ SymOp() };
      std::vector<SymOp> todo{ SymOp() };
      try {
        while ( !todo.empty() ) {
          const SymOp a = todo.back();
          todo.pop_back();
          for ( auto& g : generators ) {
            const SymOp p = a * g;
            if ( group.insert( p ).second ) {
              if ( group.size() > 192 )
                throwBadHall( hall, "too many operations" );
              todo.push_back( p );
            }
          }
        }
      } catch ( Error::LogicError& ) {
        throwBadHall( hall, "inconsistent rotation symbols (mixing hexagonal"
                      " and orthogonal-type axes)" );
      }

      //Change of basis (origin shift):
      if ( shift.has_value() ) {
        const SymOp s( identityMat(), shift.value() );
        const SymOp sinv = s.inverse();
        std::set<SymOp> shifted;
        for ( auto& op : group )
          shifted.insert( s * op * sinv );
        group.swap( shifted );
      }
      return group;
    }
  }
}

std::vector<NC::SymOp> NC::detail::rawSymOpsFromHallSymbol( StrView hall )
{
  char lattice;
  return canonicalOrder( groupFromHallSymbol( hall, &lattice ) ).ops;
}

NC::SGSymmetry::SGSymmetry( internal_t, SpaceGroup sg )
  : m_sg( sg )
{
  const char * hall = sg.hallSymbol();
  auto canonical = canonicalOrder( groupFromHallSymbol( StrView( hall ),
                                                        &m_lattice ) );
  m_ops = std::move( canonical.ops );
  m_reps = std::move( canonical.reps );
  m_centring = std::move( canonical.centring );
  const SymOp inversion( { { -1,0,0, 0,-1,0, 0,0,-1 } }, { { 0, 0, 0 } } );
  m_centrosymmetric = false;
  for ( auto& op : m_reps ) {
    bool inv = true;
    for ( auto i : ncrange( 3u ) )
      for ( auto j : ncrange( 3u ) )
        if ( op.rot( i, j ) != inversion.rot( i, j ) )
          inv = false;
    if ( inv )
      m_centrosymmetric = true;
  }
}

const NC::SGSymmetry& NC::SGSymmetry::get( SpaceGroup sg )
{
  //Created on first use, and never deleted:
  static std::mutex s_mtx;
  static std::array<const SGSymmetry*,SpaceGroupHallNumber::max_value> s_cache;
  NCRYSTAL_LOCK_GUARD( s_mtx );
  const SGSymmetry*& p = s_cache[ sg.hallNumber().get() - 1 ];
  if ( !p )
    p = new SGSymmetry( internal_t(), sg );
  return *p;
}
