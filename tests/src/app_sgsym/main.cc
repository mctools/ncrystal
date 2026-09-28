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

// Tests of the sgsym component (space groups and their symmetries).

#include "NCrystal/internal/sgsym/NCSpaceGroup.hh"
#include "NCrystal/internal/sgsym/NCSymOp.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"
#include "NCrystal/internal/phys_utils/NCEqRefl.hh"
#include <iostream>
#include <sstream>
#include <set>
#include <map>

namespace NC = NCrystal;
#define REQUIRE(x) nc_assert_always(x)

namespace NCRYSTAL_NAMESPACE {
  namespace detail {
    //Not declared in any header (only for testing):
    std::vector<SymOp> rawSymOpsFromHallSymbol( StrView );
  }
}

////////////////////////////////////////////////////////////////////////////////
// Space group table
////////////////////////////////////////////////////////////////////////////////

//Tests the SpaceGroup and SpaceGroupHallNumber classes. The complete table of
//space group settings is printed, so the reference log also guards against
//accidental changes to it.

namespace test_table {

namespace {
  void expectBad( const char * s )
  {
    try {
      NC::SpaceGroup sg( s );
    } catch ( NC::Error::BadInput& e ) {
      std::cout << "  \"" << s << "\" -> BadInput: " << e.what() << std::endl;
      return;
    }
    nc_assert_always( false );
  }
  void expectBadHall( std::uint16_t hn )
  {
    try {
      NC::SpaceGroup sg{ NC::SpaceGroupHallNumber{ hn } };
    } catch ( NC::Error::BadInput& e ) {
      std::cout << "  Hall number " << hn << " -> BadInput: " << e.what()
                << std::endl;
      return;
    }
    nc_assert_always( false );
  }
}

void run()
{
  static_assert( sizeof( NC::SpaceGroup ) == 2, "" );
  static_assert( std::is_trivially_copyable<NC::SpaceGroup>::value, "" );
  //Constructors are explicit (no implicit conversions):
  static_assert( !std::is_convertible<const char*,NC::SpaceGroup>::value, "" );
  static_assert( !std::is_convertible<std::string,NC::SpaceGroup>::value, "" );
  static_assert( !std::is_convertible<NC::StrView,NC::SpaceGroup>::value, "" );
  static_assert( !std::is_convertible<NC::SpaceGroupHallNumber,
                                      NC::SpaceGroup>::value, "" );
  static_assert( std::is_constructible<NC::SpaceGroup,std::string&&>::value,
                 "" );

  std::cout << "All settings (Hall number, setting, Hall symbol, flags):"
            << std::endl;
  std::set<unsigned> numbers, numbers_multi;
  unsigned last_number = 0;
  for ( std::uint16_t hn = 1; hn <= 530; ++hn ) {
    const NC::SpaceGroup sg{ NC::SpaceGroupHallNumber{ hn } };
    REQUIRE( sg.hallNumber().get() == hn );
    REQUIRE( sg.number() >= 1 && sg.number() <= 230 );
    REQUIRE( sg.number() >= last_number );//sorted
    REQUIRE( sg.isDefaultSetting() == ( sg.number() != last_number ) );
    last_number = sg.number();
    numbers.insert( sg.number() );
    if ( sg.hasMultipleSettings() )
      numbers_multi.insert( sg.number() );
    REQUIRE( *sg.hallSymbol() );
    //Empty choice codes only for default settings:
    REQUIRE( sg.choice() != nullptr );
    REQUIRE( *sg.choice() || sg.isDefaultSetting() );//(*: non-empty?)
    //String round trip and stream output:
    const std::string str = sg.toString();
    REQUIRE( NC::SpaceGroup( str ) == sg );
    REQUIRE( NC::SpaceGroup( std::string( str ) ) == sg );//temporary
    REQUIRE( NC::SpaceGroup( NC::StrView( str ) ) == sg );
    REQUIRE( NC::SpaceGroup( str.c_str() ) == sg );
    //Hall number -> SpaceGroup -> string -> SpaceGroup -> Hall number:
    REQUIRE( NC::SpaceGroup( str ).hallNumber().get() == hn );
    std::ostringstream ss;
    ss << sg;
    REQUIRE( ss.str() == str );
    //Bare number gives default setting:
    const std::string nstr = std::to_string( sg.number() );
    REQUIRE( ( NC::SpaceGroup( nstr ) == sg ) == sg.isDefaultSetting() );
    REQUIRE( ( NC::SpaceGroup( nstr ) < sg ) == !sg.isDefaultSetting() );
    std::cout << "  " << hn << " " << sg << " \"" << sg.hallSymbol() << "\""
              << ( sg.isDefaultSetting() ? " [default]" : "" )
              << ( sg.hasMultipleSettings() ? "" : " [single]" ) << std::endl;
  }
  REQUIRE( numbers.size() == 230 );
  std::cout << "Numbers with multiple settings: " << numbers_multi.size()
            << std::endl;
  REQUIRE( numbers_multi.size() == 90 );

  std::cout << "Examples:" << std::endl;
  for ( auto s : { "227", "227:1", "227:2", "166", "166:H", "166:R", "14",
                   "14:b1", "14:c1", "62", "62:cab", "68:1cab", "70:2", "15:-b2",
                   "1", "225" } ) {
    const NC::SpaceGroup sg( s );
    std::cout << "  \"" << s << "\" -> " << sg << " (Hall number "
              << sg.hallNumber() << ", Hall symbol \"" << sg.hallSymbol()
              << "\")" << std::endl;
  }
  REQUIRE( NC::SpaceGroup( "227" ).hallNumber().get() == 525 );
  REQUIRE( NC::SpaceGroup( "227:2" ).hallNumber().get() == 526 );

  std::cout << "Invalid input:" << std::endl;
  for ( auto s : { "", "0", "231", "999", "1000", "-1", "+1", " 227", "227 ",
                   "abc", "227:", "227:3", "227:1:2", "166:h", "62:abc",
                   "225:1", "14:b", "3:b1", ":2", "22a" } )
    expectBad( s );
  expectBadHall( 0 );
  expectBadHall( 531 );
  std::cout << "All OK" << std::endl;
}

}

////////////////////////////////////////////////////////////////////////////////
// Symmetry operations
////////////////////////////////////////////////////////////////////////////////

//Tests the SymOp class (symmetry operations).

namespace test_symop {

namespace {

  using Mat = std::array<int,9>;

  Mat matMult( const Mat& a, const Mat& b )
  {
    Mat q;
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        q[3*i+j] = a[3*i]*b[j] + a[3*i+1]*b[3+j] + a[3*i+2]*b[6+j];
    return q;
  }

  //Independent construction of the 64 allowed rotations: the 48 signed
  //permutation matrices, and the 24 rotations of 6/mmm in hexagonal axes
  //(generated here by 6-fold rotations about z, 2-fold rotations about the
  //a-axis, and inversion):
  std::set<Mat> allowedRotations()
  {
    std::set<Mat> res;
    const int perms[6][3] = { {0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1},
                              {2,1,0} };
    for ( auto& pm : perms ) {
      for ( int signs = 0; signs < 8; ++signs ) {
        Mat m{};
        for ( int i = 0; i < 3; ++i )
          m[3*i+pm[i]] = ( signs & ( 1 << i ) ) ? -1 : 1;
        res.insert( m );
      }
    }
    std::set<Mat> hex{ Mat{ { 1,0,0, 0,1,0, 0,0,1 } } };
    const Mat six{ { 1,-1,0, 1,0,0, 0,0,1 } };
    const Mat twofold_a{ { 1,-1,0, 0,-1,0, 0,0,-1 } };
    const Mat inv{ { -1,0,0, 0,-1,0, 0,0,-1 } };
    for ( bool changed = true; changed; ) {
      changed = false;
      for ( auto m : std::vector<Mat>( hex.begin(), hex.end() ) )
        for ( auto& g : { six, twofold_a, inv } )
          if ( hex.insert( matMult( m, g ) ).second )
            changed = true;
    }
    nc_assert_always( res.size() == 48 && hex.size() == 24 );
    res.insert( hex.begin(), hex.end() );
    return res;
  }

  bool entriesOK( const NC::SymOp& op )
  {
    for ( unsigned i = 0; i < 3; ++i ) {
      for ( unsigned j = 0; j < 3; ++j )
        if ( !( op.rot(i,j) >= -1 && op.rot(i,j) <= 1 ) )
          return false;
      if ( !( op.trans(i) >= 0 && op.trans(i) < NC::SymOp::tdenom ) )
        return false;
    }
    return true;
  }

  //Group generated by the given operations (only for testing):
  std::vector<NC::SymOp> closure( std::vector<NC::SymOp> gens )
  {
    std::set<NC::SymOp> g{ NC::SymOp() };
    for ( bool changed = true; changed; ) {
      changed = false;
      std::vector<NC::SymOp> cur( g.begin(), g.end() );
      for ( auto& a : cur )
        for ( auto& b : gens )
          if ( g.insert( a * b ).second )
            changed = true;
    }
    return { g.begin(), g.end() };
  }

  void expectBad( const char * s )
  {
    try {
      NC::SymOp op( s );
    } catch ( NC::Error::BadInput& e ) {
      std::cout << "  \"" << s << "\" -> BadInput: " << e.what() << std::endl;
      return;
    }
    nc_assert_always( false );
  }
}

void run()
{
  static_assert( sizeof( NC::SymOp ) == 12, "" );
  static_assert( std::is_trivially_copyable<NC::SymOp>::value, "" );
  static_assert( !std::is_convertible<const char*,NC::SymOp>::value, "" );
  static_assert( !std::is_convertible<std::string,NC::SymOp>::value, "" );
  static_assert( !std::is_convertible<NC::StrView,NC::SymOp>::value, "" );
  static_assert( std::is_constructible<NC::SymOp,std::string&&>::value, "" );

  REQUIRE( NC::SymOp().isIdentity() );
  REQUIRE( NC::SymOp().toString() == "x,y,z" );

  //All 3^9 matrices with entries in {-1,0,1}:
  const std::array<std::array<int,3>,5> translations = { {
      { { 0, 0, 0 } }, { { 12, 0, 0 } }, { { 1, 5, 23 } }, { { 8, 16, 3 } },
      { { -3, 33, 48 } } } };
  const std::set<Mat> allowed = allowedRotations();
  REQUIRE( allowed.size() == 64 );
  unsigned nvalid = 0;
  std::map<int,unsigned> ntypes;
  for ( int code = 0; code < 19683; ++code ) {
    Mat m;
    for ( int i = 0, c = code; i < 9; ++i, c /= 3 )
      m[i] = c % 3 - 1;
    const bool valid = allowed.count( m ) > 0;
    REQUIRE( valid == NC::SymOp::isAllowedRotation( m ) );
    bool accepted = true;
    try {
      NC::SymOp tmp( m, { { 0, 0, 0 } } );
    } catch ( NC::Error::BadInput& ) {
      accepted = false;
    }
    REQUIRE( valid == accepted );
    if ( !valid )
      continue;
    ++nvalid;
    for ( auto& t : translations ) {
      const NC::SymOp op( m, t );
      REQUIRE( entriesOK( op ) );
      for ( unsigned i = 0; i < 3; ++i ) {
        REQUIRE( op.trans(i) == ( ( t[i] % 24 ) + 24 ) % 24 );
        for ( unsigned j = 0; j < 3; ++j )
          REQUIRE( op.rot(i,j) == m[3*i+j] );
      }
      //Rotation type, determinant, order, powers:
      const int rt = op.rotationType();
      REQUIRE( ( rt > 0 ) == ( op.determinant() == 1 ) );
      if ( &t == &translations.front() )
        ++ntypes[rt];
      const unsigned order = op.rotationOrder();
      NC::SymOp p = op;
      for ( unsigned k = 1; k < order; ++k ) {
        REQUIRE( entriesOK( p ) );
        REQUIRE( !( p.isIdentity() ) || op.hasTranslation() );
        p = p * op;
      }
      REQUIRE( entriesOK( p ) );
      //R^order=I (translations need not vanish, e.g. for screw axes):
      REQUIRE( NC::SymOp( { { p.rot(0,0), p.rot(0,1), p.rot(0,2),
                              p.rot(1,0), p.rot(1,1), p.rot(1,2),
                              p.rot(2,0), p.rot(2,1), p.rot(2,2) } },
                          { { 0, 0, 0 } } ).isIdentity() );
      //Inverse:
      const NC::SymOp inv = op.inverse();
      REQUIRE( entriesOK( inv ) );
      REQUIRE( ( op * inv ).isIdentity() );
      REQUIRE( ( inv * op ).isIdentity() );
      REQUIRE( inv.inverse() == op );
      //String round trips and stream output:
      const std::string s = op.toString();
      REQUIRE( NC::SymOp( s ) == op );
      REQUIRE( NC::SymOp( NC::StrView( s ) ) == op );
      REQUIRE( NC::SymOp( s.c_str() ) == op );
      std::ostringstream ss;
      ss << op;
      REQUIRE( ss.str() == s );
    }
  }
  std::cout << "Accepted rotation matrices: " << nvalid << std::endl;
  REQUIRE( nvalid == 64 );

  //All products and inverses within each of the two families are valid, but
  //mixing hexagonal-only and orthogonal-only rotations is not:
  unsigned nmixed_invalid = 0;
  for ( auto& a : allowed ) {
    for ( auto& b : allowed ) {
      const NC::SymOp opa( a, { { 0, 0, 0 } } );
      const NC::SymOp opb( b, { { 0, 0, 0 } } );
      bool ok = true;
      try {
        NC::SymOp tmp = opa * opb;
        REQUIRE( entriesOK( tmp ) );
      } catch ( NC::Error::LogicError& ) {
        ok = false;
      }
      REQUIRE( ok == ( allowed.count( matMult( a, b ) ) > 0 ) );
      if ( !ok )
        ++nmixed_invalid;
    }
  }
  std::cout << "Products of two accepted rotations which are invalid: "
            << nmixed_invalid << " (of " << 64*64 << ")" << std::endl;
  try {
    NC::SymOp( "x-y,x,z" ) * NC::SymOp( "z,x,y" );
    REQUIRE( false );
  } catch ( NC::Error::LogicError& e ) {
    std::cout << "  Example: " << e.what() << std::endl;
  }
  std::cout << "Numbers per rotation type:";
  for ( auto& e : ntypes )
    std::cout << " " << e.first << ":" << e.second;
  std::cout << std::endl;

  //Associativity and group properties within cubic and hexagonal groups
  //(incl. translations):
  for ( auto gens : { std::vector<NC::SymOp>{ NC::SymOp( "z,x,y" ),
                                               NC::SymOp( "-y+1/4,x+3/4,z+1/4" ),
                                               NC::SymOp( "-x,-y,-z" ) },
                      std::vector<NC::SymOp>{ NC::SymOp( "x-y,x,z+1/6" ),
                                               NC::SymOp( "y,x,-z" ),
                                               NC::SymOp( "-x,-y,-z" ) } } ) {
    const auto g = closure( gens );
    std::cout << "Closure of";
    for ( auto& e : gens )
      std::cout << " " << e;
    std::cout << ": " << g.size() << " operations (modulo 1)" << std::endl;
    for ( auto& a : g ) {
      REQUIRE( entriesOK( a ) );
      for ( auto& b : g ) {
        REQUIRE( entriesOK( a * b ) );
        for ( std::size_t k = 0; k < g.size(); k += 7 ) {
          const auto& c = g[k];
          REQUIRE( ( a * b ) * c == a * ( b * c ) );
        }
      }
    }
  }

  //Apply:
  {
    const NC::SymOp op( "-y+1/4,x+3/4,z+1/4" );
    const NC::Vector v = op.apply( NC::Vector( 0.1, 0.2, 0.3 ) );
    REQUIRE( NC::floateq( v.x(), -0.2 + 0.25, 1e-15, 1e-15 ) );
    REQUIRE( NC::floateq( v.y(), 0.1 + 0.75, 1e-15, 1e-15 ) );
    REQUIRE( NC::floateq( v.z(), 0.3 + 0.25, 1e-15, 1e-15 ) );
    //Inverse gives back original position modulo 1 (since translations are
    //reduced modulo 1):
    const NC::Vector w = op.inverse().apply( v );
    auto mod1 = []( double x ) { return x - std::floor( x ); };
    REQUIRE( NC::floateq( mod1( w.x() ), 0.1, 1e-14, 1e-14 ) );
    REQUIRE( NC::floateq( mod1( w.y() ), 0.2, 1e-14, 1e-14 ) );
    REQUIRE( NC::floateq( mod1( w.z() ), 0.3, 1e-14, 1e-14 ) );
  }

  std::cout << "Accepted input:" << std::endl;
  for ( auto s : { "x,y,z", " X , Y , Z ", "1/2+x,1/2-y,z", "-y+x,x,z+1/6",
                   "+x,-y,+z", "x+1/2+1/4,y,z", "x-1/2,y,z", "x+3/2,y,z",
                   "x+2/4,y,z", "x+1/8,y,z", "x,y,z+1", "x,y,z-1/24",
                   "-x+y,-x,z+2/3", "y,x,-z+5/12" } )
    std::cout << "  \"" << s << "\" -> " << NC::SymOp( s ) << std::endl;

  std::cout << "Invalid input:" << std::endl;
  for ( auto s : { "", "x,y", "x,y,z,x", ",y,z", "2x,y,z", "x,x,z",
                   "x+y+z,y,z", "x+1/5,y,z", "x+0.5,y,z", "x,y,z+.5",
                   "a,b,c", "x y,y,z", "x+,y,z", "x+1/0,y,z", "x-x,y,z",
                   "x+1//2,y,z", "x,y,-z-", "--x,y,z", "x,y,1/2",
                   "x+1234567890,y,z" } )
    expectBad( s );
  try {
    NC::SymOp op( { { 2,0,0, 0,1,0, 0,0,1 } }, { { 0, 0, 0 } } );
    REQUIRE( false );
  } catch ( NC::Error::BadInput& e ) {
    std::cout << "  (2,0,0,0,1,0,0,0,1) -> BadInput: " << e.what()
              << std::endl;
  }
  std::cout << "All OK" << std::endl;
}

}

////////////////////////////////////////////////////////////////////////////////
// Symmetry of space groups
////////////////////////////////////////////////////////////////////////////////

//Tests SGSymmetry (symmetry operations derived from the Hall symbols of the
//530 space group settings): group properties, canonical ordering, total
//counts, consistency of the Laue classes with the independently implemented
//EqRefl class, and errors from the Hall symbol parser.

namespace test_symmetry {

namespace {
  using Mat = std::array<int,9>;
  Mat rotOf( const NC::SymOp& op )
  {
    Mat m;
    for ( unsigned i = 0; i < 9; ++i )
      m[i] = op.rot( i / 3, i % 3 );
    return m;
  }
  constexpr Mat identity = { { 1,0,0, 0,1,0, 0,0,1 } };
  constexpr Mat inversion = { { -1,0,0, 0,-1,0, 0,0,-1 } };

  //Does EqRefl support the setting (only default orientations of Laue
  //classes, i.e. unique axis b for monoclinic, and hexagonal or rhombohedral
  //axes for trigonal)?
  bool eqReflApplies( const NC::SpaceGroup& sg )
  {
    if ( sg.number() >= 3 && sg.number() <= 15 ) {
      const std::string c = sg.choice();
      return !c.empty() && ( c[0] == 'b' || c.substr( 0, 2 ) == "-b" );
    }
    return true;
  }

  //Compare Laue orbits of hkl (from rotation parts plus inversion) with
  //EqRefl families. Returns number of hkl points checked:
  unsigned checkEqRefl( const NC::SpaceGroup& sg,
                        const NC::SGSymmetry& sym )
  {
    const NC::EqRefl eqrefl( static_cast<int>( sg.number() ),
                             std::string( sg.choice() ) == "R" );
    unsigned n = 0;
    for ( int h = -4; h <= 4; ++h ) {
      for ( int k = -4; k <= 4; ++k ) {
        for ( int l = -4; l <= 4; ++l ) {
          if ( !h && !k && !l )
            continue;
          //Reflections transform with the transposed rotation matrix:
          std::set<std::array<int,3>> orbit;
          for ( auto& op : sym.representatives() ) {
            std::array<int,3> v;
            for ( unsigned j = 0; j < 3; ++j )
              v[j] = op.rot(0,j) * h + op.rot(1,j) * k + op.rot(2,j) * l;
            orbit.insert( v );
            orbit.insert( { { -v[0], -v[1], -v[2] } } );
          }
          std::set<std::array<int,3>> fam;
          for ( auto& e : eqrefl.getEquivalentReflections( h, k, l ) ) {
            fam.insert( { { e.h, e.k, e.l } } );
            fam.insert( { { -e.h, -e.k, -e.l } } );
          }
          REQUIRE( orbit == fam );
          ++n;
        }
      }
    }
    return n;
  }

  void expectBadHall( const char * s )
  {
    try {
      NC::detail::rawSymOpsFromHallSymbol( s );
    } catch ( NC::Error::BadInput& e ) {
      std::cout << "  \"" << s << "\" -> BadInput: " << e.what() << std::endl;
      return;
    }
    nc_assert_always( false );
  }
}

void run()
{
  std::cout << "Hall number, setting, order, number of representatives,"
            << " lattice symbol, centrosymmetric:" << std::endl;
  unsigned total_ops = 0, total_reps = 0, n_eqrefl = 0, n_eqrefl_hkl = 0;
  std::set<Mat> all_rots;
  for ( std::uint16_t hn = 1; hn <= 530; ++hn ) {
    const NC::SpaceGroup sg{ NC::SpaceGroupHallNumber{ hn } };
    const NC::SGSymmetry& sym = NC::SGSymmetry::get( sg );
    REQUIRE( &sym == &NC::SGSymmetry::get( sg ) );//cached
    REQUIRE( sym.spaceGroup() == sg );
    const auto& ops = sym.operations();
    const auto& reps = sym.representatives();
    const auto& cent = sym.centringVectors();
    REQUIRE( sym.order() == ops.size() );
    REQUIRE( ops.size() == reps.size() * cent.size() );
    total_ops += sym.order();
    total_reps += static_cast<unsigned>( reps.size() );
    //Canonical order: identity first, then blocks per centring vector:
    REQUIRE( ops.front().isIdentity() && reps.front().isIdentity()
             && cent.front().isIdentity() );
    for ( std::size_t ic = 0; ic < cent.size(); ++ic )
      for ( std::size_t ir = 0; ir < reps.size(); ++ir )
        REQUIRE( ops.at( ic * reps.size() + ir ) == cent.at( ic ) * reps.at( ir ) );
    REQUIRE( std::is_sorted( reps.begin() + 1, reps.end() ) );
    REQUIRE( std::is_sorted( cent.begin(), cent.end() ) );
    //Centring vectors are pure translations, representatives have distinct
    //rotations, and each has the smallest translation among operations with
    //that rotation:
    for ( auto& c : cent )
      REQUIRE( rotOf( c ) == identity );
    std::set<Mat> rots;
    for ( auto& r : reps ) {
      REQUIRE( rots.insert( rotOf( r ) ).second );
      for ( auto& op : ops )
        if ( rotOf( op ) == rotOf( r ) )
          REQUIRE( !( op < r ) );
    }
    all_rots.insert( rots.begin(), rots.end() );
    //Group properties (closure, inverses):
    const std::set<NC::SymOp> opset( ops.begin(), ops.end() );
    REQUIRE( opset.size() == ops.size() );
    for ( auto& a : ops ) {
      REQUIRE( opset.count( a.inverse() ) );
      for ( auto& b : ops )
        REQUIRE( opset.count( a * b ) );
    }
    //Lattice symbol and centrosymmetry:
    const char lat = sym.latticeSymbol();
    const std::size_t ncent_expected = ( lat == 'P' ? 1 : lat == 'R' ? 3
                                         : lat == 'F' ? 4 : 2 );
    REQUIRE( std::string( "PABCIRF" ).find( lat ) != std::string::npos );
    REQUIRE( cent.size() == ncent_expected );
    REQUIRE( sym.isCentrosymmetric() == ( rots.count( inversion ) > 0 ) );
    //Laue classes consistent with EqRefl:
    if ( eqReflApplies( sg ) ) {
      ++n_eqrefl;
      n_eqrefl_hkl += checkEqRefl( sg, sym );
    }
    std::cout << "  " << hn << " " << sg << " " << sym.order() << " "
              << reps.size() << " " << lat << " "
              << ( sym.isCentrosymmetric() ? 1 : 0 ) << std::endl;
  }
  std::cout << "Total number of operations: " << total_ops
            << ", representatives: " << total_reps
            << ", distinct rotations: " << all_rots.size() << std::endl;
  REQUIRE( total_ops == 7388 && total_reps == 4462 && all_rots.size() == 64 );
  std::cout << "Laue classes consistent with EqRefl for " << n_eqrefl
            << " settings (" << n_eqrefl_hkl << " hkl points in total)"
            << std::endl;

  std::cout << "Examples:" << std::endl;
  for ( auto s : { "14:b1", "15:b1", "166:H", "166:R", "194", "227:2" } ) {
    const NC::SpaceGroup sg( s );
    const auto& sym = NC::SGSymmetry::get( sg );
    std::cout << "  " << sg << " (Hall symbol \"" << sg.hallSymbol()
              << "\"):" << std::endl << "    representatives:";
    for ( auto& op : sym.representatives() )
      std::cout << " " << op;
    std::cout << std::endl << "    centring vectors:";
    for ( auto& op : sym.centringVectors() )
      std::cout << " " << op;
    std::cout << std::endl;
  }

  std::cout << "Invalid Hall symbols:" << std::endl;
  for ( auto s : { "", "X 1", "-", "P 5", "P 2q", "P 22", "P 2 2 2 2",
                   "P 3* 2\"", "P 6 3*", "P 4 (0 0", "P 1 (0 0 1 1)",
                   "P 2 (a b c)", "P 1 (0 0 1) x", "P 32z*", "P 3 3" } )
    expectBadHall( s );
  std::cout << "All OK" << std::endl;
}

}

int main()
{
  std::cout << "==== Space group table ====" << std::endl;
  test_table::run();
  std::cout << "==== Symmetry operations ====" << std::endl;
  test_symop::run();
  std::cout << "==== Symmetry of space groups ====" << std::endl;
  test_symmetry::run();
  return 0;
}
