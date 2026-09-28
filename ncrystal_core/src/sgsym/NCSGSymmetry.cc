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
#include "NCrystal/internal/utils/NCMath.hh"
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

    //Cell constraints: The metric tensor G (symmetric, with the six
    //independent entries g=(g11,g22,g33,g12,g13,g23)) must satisfy
    //R^T*G*R=G for all rotation parts R. This gives a homogeneous linear
    //system M*g=0 with small integer coefficients. A linear relation L*g=0
    //holds for all allowed metrics iff L is in the row space of M, i.e. iff
    //appending L to M does not increase the rank (computed exactly, with
    //integer arithmetic):
    using Row = std::array<std::int64_t,6>;

    unsigned exactRank( std::vector<Row> rows )
    {
      unsigned rank = 0;
      for ( unsigned col = 0; col < 6 && rank < rows.size(); ++col ) {
        std::size_t piv = rank;
        while ( piv < rows.size() && rows[piv][col] == 0 )
          ++piv;
        if ( piv == rows.size() )
          continue;
        std::swap( rows[rank], rows[piv] );
        for ( std::size_t r = rank + 1; r < rows.size(); ++r ) {
          if ( !rows[r][col] )
            continue;
          const std::int64_t f1 = rows[rank][col];
          const std::int64_t f2 = rows[r][col];
          std::int64_t g = 0;
          for ( unsigned k = 0; k < 6; ++k ) {
            rows[r][k] = rows[r][k] * f1 - rows[rank][k] * f2;
            std::int64_t v = rows[r][k] < 0 ? -rows[r][k] : rows[r][k];
            while ( v ) {//gcd, to keep numbers small
              const std::int64_t tmp = g % v;
              g = v;
              v = tmp;
            }
          }
          if ( g > 1 )
            for ( auto& e : rows[r] )
              e /= g;
        }
        ++rank;
      }
      return rank;
    }

    SGCellConstraints deriveCellConstraints( const std::vector<SymOp>& reps )
    {
      //Index of G(k,l) in g:
      auto gidx = []( unsigned k, unsigned l ) -> unsigned
      {
        if ( k == l )
          return k;
        const unsigned s = k + l;//(0,1)->1, (0,2)->2, (1,2)->3
        return s == 1 ? 3u : s == 2 ? 4u : 5u;
      };
      std::vector<Row> M;
      for ( auto& op : reps ) {
        for ( unsigned i = 0; i < 3; ++i ) {
          for ( unsigned j = i; j < 3; ++j ) {
            //(R^T G R)_ij - G_ij = sum_kl R_ki G_kl R_lj - G_ij:
            Row row{};
            for ( unsigned k = 0; k < 3; ++k )
              for ( unsigned l = 0; l < 3; ++l )
                row[gidx(k,l)] += op.rot( k, i ) * op.rot( l, j );
            row[gidx(i,j)] -= 1;
            if ( row != Row{} )
              M.push_back( row );
          }
        }
      }
      const unsigned rankM = exactRank( M );
      auto implied = [&M,rankM]( const Row& L )
      {
        std::vector<Row> tmp = M;
        tmp.push_back( L );
        return exactRank( tmp ) == rankM;
      };
      using LC = SGCellConstraints::LengthConstraint;
      using AC = SGCellConstraints::AngleConstraint;
      SGCellConstraints res;
      //g = (g11,g22,g33,g12,g13,g23):
      res.b = implied( { { -1,1,0,0,0,0 } } ) ? LC::EqualToA : LC::Free;
      res.c = implied( { { -1,0,1,0,0,0 } } ) ? LC::EqualToA : LC::Free;
      //alpha (b,c <-> g23), beta (a,c <-> g13) and gamma (a,b <-> g12). An
      //angle of 120 between vectors u and v of equal lengths means
      //2*(u.v)+u.u=0. Angles equal to alpha are only possible for equal
      //lengths a=b=c, in which case cos(angles) are proportional to g:
      auto angle = [&implied]( const Row& is90, const Row& is120,
                               const Row* eqalpha )
      {
        if ( implied( is90 ) )
          return AC::Is90;
        if ( implied( is120 ) )
          return AC::Is120;
        if ( eqalpha && implied( *eqalpha ) )
          return AC::EqualToAlpha;
        return AC::Free;
      };
      const Row beta_eq_alpha{ { 0,0,0,0,1,-1 } };
      const Row gamma_eq_alpha{ { 0,0,0,1,0,-1 } };
      res.alpha = angle( { { 0,0,0,0,0,1 } }, { { 0,1,0,0,0,2 } }, nullptr );
      res.beta = angle( { { 0,0,0,0,1,0 } }, { { 1,0,0,0,2,0 } },
                        &beta_eq_alpha );
      res.gamma = angle( { { 0,0,0,1,0,0 } }, { { 1,0,0,2,0,0 } },
                         &gamma_eq_alpha );
      //Consistency: 120 degree angles and equal angles require equal
      //lengths, and the constraints must account for all degrees of freedom
      //of the metric:
      nc_assert_always( res.alpha != AC::Is120 || res.c == LC::EqualToA );
      nc_assert_always( res.beta != AC::Is120 || res.c == LC::EqualToA );
      nc_assert_always( res.gamma != AC::Is120 || res.b == LC::EqualToA );
      if ( res.beta == AC::EqualToAlpha || res.gamma == AC::EqualToAlpha )
        nc_assert_always( res.b == LC::EqualToA && res.c == LC::EqualToA );
      const unsigned nfree = ( 1 + ( res.b == LC::Free ? 1 : 0 )
                               + ( res.c == LC::Free ? 1 : 0 )
                               + ( res.alpha == AC::Free ? 1 : 0 )
                               + ( res.beta == AC::Free ? 1 : 0 )
                               + ( res.gamma == AC::Free ? 1 : 0 ) );
      nc_assert_always( nfree == 6 - rankM );
      return res;
    }

    const char * angleName( unsigned i )
    {
      return i == 0 ? "alpha" : ( i == 1 ? "beta" : "gamma" );
    }
  }
}

void NC::SGCellConstraints::check( const CellParameters& cp ) const
{
  using LC = LengthConstraint;
  using AC = AngleConstraint;
  const double lengths[3] = { cp.a, cp.b, cp.c };
  const LC lc[3] = { LC::Free, b, c };
  const char * lnames[3] = { "a", "b", "c" };
  for ( unsigned i = 0; i < 3; ++i ) {
    if ( lengths[i] == 0.0 )
      NCRYSTAL_THROW2(BadInput,"Lattice parameter "<<lnames[i]<<" is 0 ("
                      <<( lc[i] == LC::Free
                          ? "it must be specified"
                          : "implied values must be filled in first" )<<")");
    if ( !( lengths[i] > 0.0 ) || !std::isfinite( lengths[i] ) )
      NCRYSTAL_THROW2(BadInput,"Lattice parameter "<<lnames[i]
                      <<" must be a positive number (got "<<lengths[i]<<")");
    if ( lc[i] == LC::EqualToA && lengths[i] != cp.a )
      NCRYSTAL_THROW2(BadInput,"Lattice parameters a and "<<lnames[i]
                      <<" must be equal for this space group (got a="<<cp.a
                      <<" and "<<lnames[i]<<"="<<lengths[i]<<")");
  }
  const double angles[3] = { cp.alpha, cp.beta, cp.gamma };
  const AC ac[3] = { alpha, beta, gamma };
  for ( unsigned i = 0; i < 3; ++i ) {
    if ( angles[i] == 0.0 )
      NCRYSTAL_THROW2(BadInput,"Lattice angle "<<angleName(i)<<" is 0 ("
                      <<( ac[i] == AC::Free
                          ? "it must be specified"
                          : "implied values must be filled in first" )<<")");
    if ( !( angles[i] > 0.0 && angles[i] < 180.0 ) )
      NCRYSTAL_THROW2(BadInput,"Lattice angle "<<angleName(i)<<" must be in"
                      " the range (0,180) degrees (got "<<angles[i]<<")");
    const double expected = ( ac[i] == AC::Is90 ? 90.0
                              : ac[i] == AC::Is120 ? 120.0
                              : ac[i] == AC::EqualToAlpha ? cp.alpha
                              : angles[i] );
    if ( angles[i] != expected )
      NCRYSTAL_THROW2(BadInput,"Lattice angle "<<angleName(i)<<" must be "
                      <<( ac[i] == AC::EqualToAlpha ? "equal to alpha"
                          : ac[i] == AC::Is90 ? "90 degrees" : "120 degrees" )
                      <<" for this space group (got "<<angles[i]<<")");
  }
  if ( beta == AC::EqualToAlpha && gamma == AC::EqualToAlpha
       && !( cp.alpha < 120.0 ) )
    NCRYSTAL_THROW2(BadInput,"Lattice angles alpha=beta=gamma must be less"
                    " than 120 degrees (got "<<cp.alpha<<")");
}

void NC::SGCellConstraints::complete( CellParameters& cp ) const
{
  using LC = LengthConstraint;
  using AC = AngleConstraint;
  if ( b == LC::EqualToA && cp.b == 0.0 )
    cp.b = cp.a;
  if ( c == LC::EqualToA && cp.c == 0.0 )
    cp.c = cp.a;
  double* angles[3] = { &cp.alpha, &cp.beta, &cp.gamma };
  const AC ac[3] = { alpha, beta, gamma };
  for ( unsigned i = 0; i < 3; ++i ) {
    if ( *angles[i] != 0.0 )
      continue;
    if ( ac[i] == AC::Is90 )
      *angles[i] = 90.0;
    else if ( ac[i] == AC::Is120 )
      *angles[i] = 120.0;
    else if ( ac[i] == AC::EqualToAlpha )
      *angles[i] = cp.alpha;//(alpha is never EqualToAlpha itself)
  }
  check( cp );
}

std::ostream& NC::operator<<( std::ostream& os, const SGCellConstraints& cc )
{
  using LC = SGCellConstraints::LengthConstraint;
  using AC = SGCellConstraints::AngleConstraint;
  os << "a, b" << ( cc.b == LC::EqualToA ? "=a" : "" )
     << ", c" << ( cc.c == LC::EqualToA ? "=a" : "" );
  const AC ac[3] = { cc.alpha, cc.beta, cc.gamma };
  for ( unsigned i = 0; i < 3; ++i ) {
    os << ", " << angleName( i );
    switch ( ac[i] ) {
    case AC::Free: break;
    case AC::Is90: os << "=90"; break;
    case AC::Is120: os << "=120"; break;
    case AC::EqualToAlpha: os << "=alpha"; break;
    }
  }
  return os;
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
  m_cellconstraints = deriveCellConstraints( m_reps );
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    constexpr double site_merge_tol = 5e-4;
    constexpr double site_distinct_tol = 1e-2;

    double wrapCoord( double x )
    {
      double r = x - std::floor( x );
      if ( !( r < 1.0 ) )
        r = 0.0;//x-floor(x) can round to exactly 1.0
      return r == 0.0 ? 0.0 : r;//also maps -0 to 0
    }

    Vector wrapPos( const Vector& v )
    {
      return Vector( wrapCoord( v.x() ), wrapCoord( v.y() ),
                     wrapCoord( v.z() ) );
    }

    //Periodic difference of each coordinate, in [-0.5,0.5]:
    Vector periodicDiff( const Vector& a, const Vector& b )
    {
      auto f = []( double d ) { return d - std::round( d ); };
      return Vector( f( a.x() - b.x() ), f( a.y() - b.y() ),
                     f( a.z() - b.z() ) );
    }

    //Largest periodic coordinate difference:
    double periodicDist( const Vector& a, const Vector& b )
    {
      const Vector d = periodicDiff( a, b );
      return ncmax( ncabs( d.x() ), ncabs( d.y() ), ncabs( d.z() ) );
    }

    [[noreturn]] void throwNearSpecial( const Vector& site, double dist )
    {
      NCRYSTAL_THROW2(BadInput,"The site ("<<fmt(site.x())<<", "
                      <<fmt(site.y())<<", "<<fmt(site.z())<<") is close to a"
                      " special position but not on it (symmetry-equivalent"
                      " positions are only "<<fmt(dist,"%.2g")<<" apart in"
                      " fractional coordinates). For sites on special"
                      " positions, specify coordinates with sufficient"
                      " precision, preferably as exact fractions like 1/3.");
    }
  }
}

NC::SGSiteOrbit NC::expandSiteToOrbit( const SGSymmetry& sym,
                                       const Vector& input_site )
{
  if ( !std::isfinite( input_site.x() ) || !std::isfinite( input_site.y() )
       || !std::isfinite( input_site.z() ) )
    NCRYSTAL_THROW2(BadInput,"Invalid site coordinates: ("<<input_site.x()
                    <<", "<<input_site.y()<<", "<<input_site.z()<<")");
  const Vector site = wrapPos( input_site );
  const auto& ops = sym.operations();

  //Symmetrise: average (the nearest periodic copies of) the images of the
  //site under the operations which map it onto itself. These operations
  //(with their translations adjusted by the corresponding lattice vectors)
  //form a group, so the average is exactly invariant under them:
  Vector sum( 0.0, 0.0, 0.0 );
  unsigned nstab = 0;
  for ( auto& op : ops ) {
    const Vector img = op.apply( site );
    const Vector d = periodicDiff( img, site );
    const double dist = ncmax( ncabs( d.x() ), ncabs( d.y() ),
                               ncabs( d.z() ) );
    if ( dist < site_merge_tol ) {
      sum += site + d;
      ++nstab;
    } else if ( dist < site_distinct_tol ) {
      throwNearSpecial( site, dist );
    }
  }
  nc_assert_always( nstab >= 1 );//at least the identity
  const Vector symsite = wrapPos( sum * ( 1.0 / nstab ) );

  //Expand the symmetrised site, keeping the first occurrence of each
  //distinct position (all pairs of distinct positions must be clearly
  //separated):
  SGSiteOrbit res;
  for ( auto& op : ops ) {
    const Vector img = wrapPos( op.apply( symsite ) );
    bool found = false;
    for ( auto& p : res.positions ) {
      const double dist = periodicDist( img, p );
      if ( dist < site_merge_tol ) {
        found = true;
        break;
      }
      if ( dist < site_distinct_tol )
        throwNearSpecial( site, dist );
    }
    if ( !found )
      res.positions.push_back( img );
  }
  const std::size_t npos = res.positions.size();
  nc_assert_always( npos * nstab == ops.size() );
  res.siteSymmetryOrder = nstab;
  return res;
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
