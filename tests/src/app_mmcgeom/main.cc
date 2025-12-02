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

#include "NCrystal/internal/minimc/NCMMC_Geom.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include <iostream>

#define REQUIRE nc_assert_always
#define NCTESTPREC 1e-13

#define NCTESTPREC_SPHERE 1e-9 // due to squared numbers, max precision for some
                               //sphere tests are only ~half of normal.
#ifndef NDEBUG
#  define REQUIREFLTEQ(x,y) assert(NC::floateq(x,y,NCTESTPREC,NCTESTPREC))
#else
#  define REQUIREFLTEQ(x,y) nc_assert_always(NC::floateq(x,y,NCTESTPREC,NCTESTPREC))
#endif
#ifndef NDEBUG
#  define REQUIREFLTEQ_SPHERE(x,y) assert(NC::floateq(x,y,NCTESTPREC_SPHERE,NCTESTPREC_SPHERE))
#else
#  define REQUIREFLTEQ_SPHERE(x,y) nc_assert_always(NC::floateq(x,y,NCTESTPREC_SPHERE,NCTESTPREC_SPHERE))
#endif

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

namespace {

  void initBasket( NCMMC::NeutronBasket& b,
                   const NC::Vector& pos,
                   const NC::Vector& dir )
  {
    b.nused = 1;
    b.x[0] = pos.x();
    b.y[0] = pos.y();
    b.z[0] = pos.z();
    auto u = dir.unit();
    b.ux[0] = u.x();
    b.uy[0] = u.y();
    b.uz[0] = u.z();
    b.w[0] = 1.0;
    b.ekin[0] = 0.025;
  }

  double distToVolEntry( const NCMMC::Geometry& geom,
                         NC::Vector pos,
                         NC::Vector dir )
  {
    NCMMC::NeutronBasket b;
    initBasket( b, pos, dir );
    double res[1];
    geom.distToVolumeEntry( b, res );
    return res[0];
  }

  double distToVolExit( const NCMMC::Geometry& geom,
                        NC::Vector pos,
                        NC::Vector dir )
  {
    NCMMC::NeutronBasket b;
    initBasket( b, pos, dir );
    double res[1];
    geom.distToVolumeExit( b, res );
    return res[0];
  }

  class TestSlab final {
    double m_dz;
  public:
    TestSlab( double dz ) : m_dz(dz) {}

    bool pointIsInside( const NC::Vector& pos ) const
    {
      return NC::ncabs(pos.z()) <= m_dz;
    }

    bool pointIsCompletelyInside( const NC::Vector& pos ) const
    {
      return NC::ncabs(pos.z()) < m_dz;
    }

    double distToExit( const NC::Vector& pos, const NC::Vector& dir_raw ) const
    {
      nc_assert_always(pointIsInside(pos));
      const double pz = pos.z();
      const double uz = dir_raw.unit().z();
      if ( uz > 0.0 )
        return (m_dz-pz)/uz;
      if ( uz < 0.0 )
        return -(m_dz+pz)/uz;
      return NC::kInfinity;
    }

    double distToEntry( const NC::Vector& pos, const NC::Vector& dir_raw ) const
    {
      const double pz = pos.z();
      if ( NC::ncabs(pz)<m_dz )
        return 0.0;
      const double uz = dir_raw.unit().z();
      if ( NC::ncabs(pz)==m_dz ) {
        //edge
        return pz*uz > 0.0 ? -1.0 : -0.0;
      }
      //outside:
      if ( pz > 0.0 )
        return ( uz >= 0.0 ? -1.0 : (pz-m_dz)/(-uz) );
      return ( uz <= 0.0 ? -1.0 : (-m_dz-pz)/uz );
    }
  };

  class TestBox final {
    double m_dx, m_dy, m_dz;
    struct Face {
      NC::Vector center;
      NC::Vector normal;
      NC::Vector uspan1, uspan2;
      const char* name = nullptr;
    };
    void addFace( NC::Vector C, NC::Vector s1, NC::Vector s2, const char * nn )
    {
      m_faces.emplace_back();
      auto& f = m_faces.back();
      f.center = C;
      f.normal = C.unit();
      f.uspan1 = s1.unit()/s1.mag();
      f.uspan2 = s2.unit()/s2.mag();
      f.name = nn;
    }
    std::vector<Face> m_faces;
    NC::Optional<double> dist_to_face( const Face& f,
                                       const NC::Vector& pos,
                                       NC::Vector dir ) const
    {
      //We assume dir is normalised.
      const double speed = -f.normal.dot(dir);
      const NC::Vector posmC = pos-f.center;
      if ( speed == 0.0 ) {
        //Moving orthogonal to the face normal. We count this as a miss, unless
        //we are moving inside the face.
        if ( f.normal.dot(posmC) == 0.0 ) {
          const double s1 = f.uspan1.dot(posmC);
          const double s2 = f.uspan2.dot(posmC);
          if ( NC::ncmax(NC::ncabs(s1),NC::ncabs(s2)) < 1.0 )
            return 0.0;
        }
        return NC::NullOpt;//miss
      }
      const double t = f.normal.dot(posmC) / speed;
      if ( t < 0.0 )
        return NC::NullOpt;//miss
      const NC::Vector local_p = posmC+(dir*t);
      const double s1 = f.uspan1.dot(local_p);
      const double s2 = f.uspan2.dot(local_p);
      if ( NC::ncabs(s1)>1.0 || NC::ncabs(s2) > 1.0 )
        return NC::NullOpt;//miss
      return t;
    }
  public:
    TestBox( double dx, double dy, double dz )
      : m_dx(dx), m_dy(dy), m_dz(dz)
    {
      //center, normal, span1, span2.
      addFace( { -dx, 0, 0 }, {0,0,dz}, {0,dy,0}, "-dx" );
      addFace( { +dx, 0, 0 }, {0,0,dz}, {0,dy,0}, "+dx" );
      addFace( { 0, -dy, 0 }, {dx,0,0}, {0,0,dz}, "-dy" );
      addFace( { 0, +dy, 0 }, {dx,0,0}, {0,0,dz}, "+dy" );
      addFace( { 0, 0, -dz }, {dx,0,0}, {0,dy,0}, "-dz" );
      addFace( { 0, 0, +dz }, {dx,0,0}, {0,dy,0}, "+dz" );
    }

    bool pointIsInside( const NC::Vector& pos ) const
    {
      return ( NC::ncabs(pos.x()) <= m_dx
               && NC::ncabs(pos.y()) <= m_dy
               && NC::ncabs(pos.z()) <= m_dz );
    }

    bool pointIsCompletelyInside( const NC::Vector& pos ) const
    {
      return ( NC::ncabs(pos.x()) < m_dx
               && NC::ncabs(pos.y()) < m_dy
               && NC::ncabs(pos.z()) < m_dz );
    }

    double distToExit( const NC::Vector& pos, const NC::Vector& dir_raw ) const
    {
      nc_assert_always(pointIsInside(pos));
      NC::Vector dir = dir_raw.unit();
      NC::Optional<double> dist;
      for ( auto& f : m_faces ) {
        auto d = dist_to_face( f, pos, dir );
        if ( !d.has_value() )
          continue;
        if ( d.value() == 0.0 ) {
          //We are ON the face. Result will depend on whether or not we are
          //heading into or out of the box.
          if ( dir.dot(f.normal) > 0.0 ) {
            return 0.0;//clearly heading out
          } else {
            //Heading in (or on an edge and heading out). Ignore this one.
            continue;
          }
        }
        if ( !dist.has_value() || d.value() < dist.value() )
          dist = d;
      }
      nc_assert_always(dist.has_value());
      return dist.value();

    }

    double distToEntry( const NC::Vector& pos, const NC::Vector& dir_raw ) const
    {
      //First a simple check if we are completely inside.
      if ( pointIsCompletelyInside( pos ) )
        return 0.0;
      //Ok, pos is outside or on the faces, check the intersection with the
      //faces, and find the shortest. If on the face, check direction vs.  face
      //normal.
      NC::Vector dir = dir_raw.unit();
      NC::Optional<double> dist;
      for ( auto& f : m_faces ) {
        auto d = dist_to_face( f, pos, dir );
        if ( !d.has_value() )
          continue;
        if ( d.value() == 0.0 ) {
          //We are ON the face. Result will depend on whether or not we are
          //heading into or out of the box.
          if ( dir.dot(f.normal) > 0.0 )
            return -1;//clearly heading out
          //Heading in (unless also on some other face and heading out, which
          //can happen on an edge of the box - that will be caught by the check
          //for that face.):
          dist = 0.0;
        }

        if ( !dist.has_value() || d.value() < dist.value() )
          dist = d;
      }
      return dist.value_or(-1);
    }
  };
  struct BoxTestCase final {
    NC::Vector pos;
    NC::Vector dir;
    NC::Optional<double> expected_disttoentry;
    NC::Optional<double> expected_disttoexit;
    BoxTestCase( NC::Vector p, NC::Vector d,
                 NC::Optional<double> d2entry = NC::NullOpt,
                 NC::Optional<double> d2exit = NC::NullOpt )
      : pos(p), dir(d.unit()),
        expected_disttoentry(d2entry),
        expected_disttoexit(d2exit) {}
  };
  using SlabTestCase = BoxTestCase;
  using SphereTestCase = BoxTestCase;

  class TestSphere final {
    double m_r;
    double m_rsq;
  public:
    TestSphere( double r ) : m_r(r), m_rsq(r*r) {}

    bool pointIsInside( const NC::Vector& pos ) const
    {
      return NC::ncabs(pos.mag2()) <= m_rsq;
    }

    bool pointIsCompletelyInside( const NC::Vector& pos ) const
    {
      return NC::ncabs(pos.mag2()) < m_rsq;
    }

    double distToExit( const NC::Vector& pos, const NC::Vector& dir_raw ) const
    {
      nc_assert_always(pointIsInside(pos));

      const NC::Vector dir = dir_raw.unit();
      const double B = 2.0 * pos.dot(dir);
      const double C = pos.mag2()-m_rsq;
      const double D = B*B - 4*C;//A=1
      REQUIRE( D>=-1e-15 );
      const double sqrtD = std::sqrt(NC::ncmax(0.0,D));
      const double t1 = ( -B - sqrtD )*0.5;
      const double t2 = ( -B + sqrtD )*0.5;
      const double tmin = NC::ncmin(t1,t2);
      const double tmax = NC::ncmax(t1,t2);
      if ( tmax <= 0.0 )
        return 0.0;//near edge, heading out
      nc_assert_always(tmin<=0.0);
      nc_assert_always(tmax>=0.0);
      return tmax;
    }

    double distToEntry( const NC::Vector& pos, const NC::Vector& dir_raw ) const
    {
      if ( pointIsCompletelyInside(pos) )
        return 0.0;
      const NC::Vector dir = dir_raw.unit();
      const double posmag2 = pos.mag2();
      if ( posmag2 > m_rsq*(1.0-1e15) && posmag2<m_rsq*(1.0+1e-15) ) {
        //edge, decide based on normal if headed in or out.
        return pos.dot(dir) < 0.0 ? 0.0 : -1.0 ;
      }
      const double B = 2.0 * pos.dot(dir);
      const double C = pos.mag2()-m_rsq;
      const double D = B*B - 4*C;//A=1
      if ( D <= 0.0 )
        return -1.0;//miss
      const double sqrtD = std::sqrt(D);
      const double t1 = ( -B - sqrtD )*0.5;
      const double t2 = ( -B + sqrtD )*0.5;
      const double tmin = NC::ncmin(t1,t2);
      const double tmax = NC::ncmax(t1,t2);
      if ( tmax < 0.0 )
        return -1.0;//intersection in backwards direction, miss
      if ( tmin > 0.0 )
        return tmin;//intersects twice

      //Presumably we are on the edge if we can get here. Repeat edge trick from
      //above:
      nc_assert_always( pos.mag() > m_r*0.99999 );
      return pos.dot(dir) < 0.0 ? 0.0 : -1.0 ;
    }
  };

}

void testSlabCases()
{
  using V = NC::Vector;
  const double dz(3.0);
  auto geom = NCMMC::createGeometry( "slab;dz=3.0" );
  TestSlab testslab( dz );
  auto rng  = NC::getRNG();
  std::vector<SlabTestCase> test_cases;
  //Manually added (and calculated) test cases:
  test_cases.emplace_back( V( 0.0, 0.0, -10 ), V( 0,0,1 ), 7 );
  test_cases.emplace_back( V( 0.0, 0.0, -2 ), V( 0,0,1 ), 0, 5 );
  test_cases.emplace_back( V( 0.0, 0.0, -2 ), V( 0,1,1 ), 0, 5*std::sqrt(2) );

  test_cases.emplace_back( V( 0.0, 0.0, dz ), V( 0,1,0 ), 0, NC::kInfinity );
  test_cases.emplace_back( V( 0.0, 0.0, dz ), V( 0,1,1 ), -1, 0 );
  test_cases.emplace_back( V( 0.0, 0.0, dz ), V( 0,1,-1 ), 0, 6*std::sqrt(2) );

  //And a bunch of generated ones:
  for ( std::size_t i = 0; i < 100000; ++i ) {
    double use_dz;
    if ( i%3 == 0 )
      use_dz = 2 * dz * ( -1.0 + 2.0 * rng->generate() );
    else if ( i%3 == 1 )
      use_dz = dz;
    else
      use_dz = -dz;
    double xyscale = 10.0;
    if ( rng->generate() < 0.1 )
      xyscale = 1e99;
    V pos = V( xyscale * ( -1.0 + 2.0 * rng->generate() ),
               xyscale * ( -1.0 + 2.0 * rng->generate() ),
               use_dz );
    V dir;
    if ( rng->generate() < 0.7 ) {
      dir = NC::randIsotropicDirection( rng );
    } else {
      auto puc = randPointOnUnitCircle( rng );
      dir.set( puc.first,puc.second,0.0 );
    }
    test_cases.emplace_back(pos,dir);
  }

  const bool verbose = false;

  for ( auto& tcase : test_cases ) {
    const V& pos = tcase.pos;
    const V& dir = tcase.dir;
    const bool testexit = testslab.pointIsInside( pos );

    const double refdist_2entry = testslab.distToEntry( pos, dir );
    const double dist_2entry = distToVolEntry( geom, pos, dir );
    const bool has_manualref_d2entry = tcase.expected_disttoentry.has_value();
    const bool has_manualref_d2exit = tcase.expected_disttoexit.has_value();

    const double refdist_2exit = ( testexit
                                   ? testslab.distToExit( pos, dir )
                                   : -17.0 );
    const double dist_2exit = ( testexit
                                ? distToVolExit( geom, pos, dir )
                                : -17.0 );
    const bool has_manualref = has_manualref_d2entry || has_manualref_d2exit;

    if ( verbose || has_manualref ) {
      std::cout<<"Slab case: pos = "<<pos<<"  dir="<<dir<<std::endl;
      std::cout<<"   DistToEntry:"<<std::endl;
      std::cout<<"               refdist = "<<NC::fmtg(refdist_2entry)<<std::endl;
      std::cout<<"                  dist = "<<NC::fmtg(dist_2entry)<<std::endl;
      if ( has_manualref_d2entry )
        std::cout<<"     expected (manual) = "
                 <<NC::fmtg(tcase.expected_disttoentry.value())<<std::endl;
      if ( testexit ) {
        std::cout<<"   DistToExit:"<<std::endl;
        std::cout<<"               refdist = "<<NC::fmtg(refdist_2exit)<<std::endl;
        std::cout<<"                  dist = "<<NC::fmtg(dist_2exit)<<std::endl;
        if ( has_manualref_d2exit )
          std::cout<<"     expected (manual) = "
                   <<NC::fmtg(tcase.expected_disttoexit.value())<<std::endl;
      }
      std::cout<<std::endl;
    }
    REQUIREFLTEQ( refdist_2entry,
                  tcase.expected_disttoentry.value_or(refdist_2entry));
    REQUIREFLTEQ( refdist_2entry, dist_2entry );
    if ( testexit ) {
      REQUIREFLTEQ( refdist_2exit,
                  tcase.expected_disttoexit.value_or(refdist_2exit));
      REQUIREFLTEQ( refdist_2exit, dist_2exit );
    }
  }
  std::cout<<"Tested a total of "<<test_cases.size()
           <<" Slab test cases!"<<std::endl;
}


void testBoxCases()
{
  using V = NC::Vector;
  const double dx(2.0), dy(3.0), dz(5.0);
  auto geom = NCMMC::createGeometry( "box;dx=2.0;dy=3.0;dz=5.0" );
  TestBox testbox( dx,dy,dz );
  auto rng  = NC::getRNG();
  std::vector<BoxTestCase> test_cases;
  //Manually added (and calculated) test cases:
  test_cases.emplace_back( V( -dx, 99, 0 ), V( 0,-1,0 ), 99-dy );
  test_cases.emplace_back( V( -dx, 99, 0 ), V( 0,1,0 ), -1.0 );

  test_cases.emplace_back( V( -dx, 0, 0 ), V( -1,0,0 ), -1.0, 0.0 );
  test_cases.emplace_back( V( -dx, 0, 0 ), V( -1,1,0 ), -1.0, 0.0 );
  test_cases.emplace_back( V( -dx, 0, 0 ), V( 1,1,0 ), 0.0, std::sqrt(18) );
  test_cases.emplace_back( V( -dx-2.5, 0, 0 ), V( 1,0,0 ), 2.5 );
  test_cases.emplace_back( V( -dx-2.5, 0, 0 ), V( -1,0,0 ), -1 );

  test_cases.emplace_back( V( dx, 0, 0 ), V( -1,1,0 ), 0, std::sqrt(18) );
  test_cases.emplace_back( V( dx, 0, 0 ), V( 1,1,0 ), -1, 0.0 );
  test_cases.emplace_back( V( dx+2.5, 0, 0 ), V( -1,0,0 ), 2.5 );
  test_cases.emplace_back( V( dx+2.5, 0, 0 ), V( 1,0,0 ), -1 );

  test_cases.emplace_back( V( 0, -dy, 0 ), V( 1,-1,0 ), -1.0, 0.0 );
  test_cases.emplace_back( V( 0, -dy, 0 ), V( 1,1,0 ), 0.0, std::sqrt(8) );
  test_cases.emplace_back( V( 0, -dy-2.5, 0 ), V( 0,1,0 ), 2.5 );
  test_cases.emplace_back( V( 0, -dy-2.5, 0 ), V( 0,-1,0 ), -1 );

  test_cases.emplace_back( V( 0, dy, 0 ), V( 1,-1,0 ), 0.0, std::sqrt(8) );
  test_cases.emplace_back( V( 0, dy, 0 ), V( 1,1,0 ), -1.0, 0.0 );
  test_cases.emplace_back( V( 0, dy+2.5, 0 ), V( 0,-1,0 ), 2.5 );
  test_cases.emplace_back( V( 0, dy+2.5, 0 ), V( 0,1,0 ), -1 );

  test_cases.emplace_back( V( 0, 0, -dz ), V( 1,0,-1 ), -1.0, 0.0 );
  test_cases.emplace_back( V( 0, 0, -dz ), V( 1,0,1 ), 0.0, std::sqrt(8) );
  test_cases.emplace_back( V( 0, 0, -dz-2.5 ), V( 0,0,1 ), 2.5 );
  test_cases.emplace_back( V( 0, 0, -dz-2.5 ), V( 0,0,-1 ), -1 );

  test_cases.emplace_back( V( 0, 0, dz ), V( 1,0,-1 ), 0.0, std::sqrt(8) );
  test_cases.emplace_back( V( 0, 0, dz ), V( 1,0,1 ), -1.0, 0.0 );
  test_cases.emplace_back( V( 0, 0, dz+2.5 ), V( 0,0,-1 ), 2.5 );
  test_cases.emplace_back( V( 0, 0, dz+2.5 ), V( 0,0,1 ), -1 );


  test_cases.emplace_back( V( -dx*1.01, -dy*1.01, -dz*1.01 ), V( 1,1,1 ) );

  test_cases.emplace_back( V( -dx, -dy, -dz ), V( 1,1,1 ), 0, std::sqrt(3*16) );
  test_cases.emplace_back( V( -dx, -dy, -dz ), V( 0,0,1 ), 0, 10 );
  test_cases.emplace_back( V( -dx, -dy, -dz ), V( 0,1,0 ), 0, 6 );
  test_cases.emplace_back( V( -dx, -dy, -dz ), V( 1,0,0 ), 0, 4 );
  test_cases.emplace_back( V( -dx, -dy, -dz ), V( 1,1,0 ), 0, std::sqrt(32));
  test_cases.emplace_back( V( -dx, -dy, -dz ), V( 1,0,1 ), 0, std::sqrt(2*16) );
  test_cases.emplace_back( V( -dx, -dy, -dz ), V( 0,1,1 ), 0, std::sqrt(2*36) );
  test_cases.emplace_back( V( 0, 0, 0 ), V( 0, 0, 1 ), 0, 5 );
  test_cases.emplace_back( V( -1, -1, 3 ), V( -1,-2, 1 ), 0, std::sqrt(6) );

  test_cases.emplace_back( V( 0, 0, 0 ), V( 2,3,5 ), 0, std::sqrt(38) );
  test_cases.emplace_back( V( 0, 0, 0 ), V( -2,3,5 ), 0, std::sqrt(38) );
  test_cases.emplace_back( V( 0, 0, 0 ), V( 2,-3,5 ), 0, std::sqrt(38) );
  test_cases.emplace_back( V( 0, 0, 0 ), V( 2,3,-5 ), 0, std::sqrt(38) );
  test_cases.emplace_back( V( 0, 0, 0 ), V( -2,-3,5 ), 0, std::sqrt(38) );
  test_cases.emplace_back( V( 0, 0, 0 ), V( -2,3,-5 ), 0, std::sqrt(38) );
  test_cases.emplace_back( V( 0, 0, 0 ), V( 2,-3,-5 ), 0, std::sqrt(38) );

  test_cases.emplace_back( V( 2,3,5 ), V(1,1,1), -1, 0 );
  test_cases.emplace_back( V( -2,3,5 ), V(-1,1,1), -1, 0 );
  test_cases.emplace_back( V( 2,-3,5 ), V(1,-1,1), -1, 0 );
  test_cases.emplace_back( V( 2,3,-5 ), V(1,1,-1), -1, 0 );
  test_cases.emplace_back( V( -2,-3,5 ), V(-1,-1,1), -1, 0 );
  test_cases.emplace_back( V( -2,3,-5 ), V(-1,1,-1), -1, 0 );
  test_cases.emplace_back( V( 2,-3,-5 ), V(1,-1,-1), -1, 0 );


  //sitting in the center of a face, direction in the plane of the face.:
  test_cases.emplace_back( V( -dx, 0, 0 ), V( 0,1,0 ), 0.0 );
  test_cases.emplace_back( V( dx, 0, 0 ), V( 0,1,1 ), 0.0 );
  test_cases.emplace_back( V( 0, -dy, 0 ), V( 1,0,1 ), 0.0 );
  test_cases.emplace_back( V( 0, dy, 0 ), V( 1,0,2 ), 0.0 );
  test_cases.emplace_back( V( 0, 0, -dz ), V( 0,1,0 ), 0.0 );
  test_cases.emplace_back( V( 0, 0, dz ), V( -1,1, 0 ), 0.0 );

  //And a bunch of generated ones:
  for ( std::size_t i = 0; i < 100000; ++i ) {
    V pos = V( 2 * dx * ( -1.0 + 2.0 * rng->generate() ),
               2 * dy * ( -1.0 + 2.0 * rng->generate() ),
               2 * dz * ( -1.0 + 2.0 * rng->generate() ) );
    V dir = NC::randIsotropicDirection( rng );
    test_cases.emplace_back(pos,dir);
  }

  const bool verbose = false;

  for ( auto& tcase : test_cases ) {
    const V& pos = tcase.pos;
    const V& dir = tcase.dir;
    const bool testexit = testbox.pointIsInside( pos );


    const double refdist_2entry = testbox.distToEntry( pos, dir );
    const double dist_2entry = distToVolEntry( geom, pos, dir );
    const bool has_manualref_d2entry = tcase.expected_disttoentry.has_value();
    const bool has_manualref_d2exit = tcase.expected_disttoexit.has_value();

    const double refdist_2exit = ( testexit
                                   ? testbox.distToExit( pos, dir )
                                   : -17.0 );
    const double dist_2exit = ( testexit
                                ? distToVolExit( geom, pos, dir )
                                : -17.0 );
    const bool has_manualref = has_manualref_d2entry || has_manualref_d2exit;

    if ( verbose || has_manualref ) {
      std::cout<<"Box case: pos = "<<pos<<"  dir="<<dir<<std::endl;
      std::cout<<"   DistToEntry:"<<std::endl;
      std::cout<<"               refdist = "
               <<NC::fmtg(refdist_2entry)<<std::endl;
      std::cout<<"                  dist = "
               <<NC::fmtg(dist_2entry)<<std::endl;
      if ( has_manualref_d2entry )
        std::cout<<"     expected (manual) = "
                 <<NC::fmtg(tcase.expected_disttoentry.value())<<std::endl;
      if ( testexit ) {
        std::cout<<"   DistToExit:"<<std::endl;
        std::cout<<"               refdist = "
                 <<NC::fmtg(refdist_2exit)<<std::endl;
        std::cout<<"                  dist = "
                 <<NC::fmtg(dist_2exit)<<std::endl;
        if ( has_manualref_d2exit )
          std::cout<<"     expected (manual) = "
                   <<NC::fmtg(tcase.expected_disttoexit.value())<<std::endl;
      }
      std::cout<<std::endl;
    }
    REQUIREFLTEQ( refdist_2entry,
                  tcase.expected_disttoentry.value_or(refdist_2entry));
    REQUIREFLTEQ( refdist_2entry, dist_2entry );
    if ( testexit ) {
      REQUIREFLTEQ( refdist_2exit,
                  tcase.expected_disttoexit.value_or(refdist_2exit));
      REQUIREFLTEQ( refdist_2exit, dist_2exit );
    }
  }
  std::cout<<"Tested a total of "<<test_cases.size()
           <<" Box test cases!"<<std::endl;
}

void testSphereCases()
{
  using V = NC::Vector;
  const double r(3.0);
  auto geom = NCMMC::createGeometry( "sphere;r=3.0" );
  TestSphere testsphere( r );
  auto rng  = NC::getRNG();
  std::vector<SphereTestCase> test_cases;
  //Manually added (and calculated) test cases:
  test_cases.emplace_back( V( 0.0, 0.0, -10 ), V( 0,0,1 ), 7 );
  test_cases.emplace_back( V( 0.0, 0.0, -2 ), V( 0,0,1 ), 0, 5 );
  test_cases.emplace_back( V( 0.0, 0.0, -3 ), V( 0,1,1 ), 0, std::sqrt(18) );


  test_cases.emplace_back( V( 0, 0, -3 ), V( 0,1, 0 ), -1, 0 );
  test_cases.emplace_back( V( 0, -3, 0 ), V( 1,0, 0 ), -1, 0 );
  test_cases.emplace_back( V( -3, 0, 0 ), V( 0,0, 1 ), -1, 0 );

  test_cases.emplace_back( V( 0, 0, -3 ), V( 0,0,1 ), 0, 6 );
  test_cases.emplace_back( V( 0, -3, 0 ), V( 0,1,0 ), 0, 6 );
  test_cases.emplace_back( V( -3, 0, 0 ), V( 1,0,0 ), 0, 6 );

  //And a bunch of generated ones (but not exactly on the edge, those we handle
  //manually above):
  for ( std::size_t i = 0; i < 100000; ++i ) {
    V pos = V( 2 * r * ( -1.0 + 2.0 * rng->generate() ),
               2 * r * ( -1.0 + 2.0 * rng->generate() ),
               2 * r * ( -1.0 + 2.0 * rng->generate() ) );
    if ( rng->generate() < 0.1 )
      pos = pos.unit() * r;
    V dir = NC::randIsotropicDirection( rng );
    test_cases.emplace_back(pos,dir);
  }

  const bool verbose = false;

  for ( auto& tcase : test_cases ) {
    const V& pos = tcase.pos;
    const V& dir = tcase.dir;
    const bool testexit = testsphere.pointIsInside( pos );
    const double refdist_2entry = testsphere.distToEntry( pos, dir );
    const double dist_2entry = distToVolEntry( geom, pos, dir );
    const bool has_manualref_d2entry = tcase.expected_disttoentry.has_value();
    const bool has_manualref_d2exit = tcase.expected_disttoexit.has_value();

    const double refdist_2exit = ( testexit
                                   ? testsphere.distToExit( pos, dir )
                                   : -17.0 );
    const double dist_2exit = ( testexit
                                ? distToVolExit( geom, pos, dir )
                                : -17.0 );
    const bool has_manualref = has_manualref_d2entry || has_manualref_d2exit;

    if ( verbose || has_manualref ) {
      std::cout<<"Sphere case: pos = "<<pos<<"  dir="<<dir<<std::endl;
      std::cout<<"   DistToEntry:"<<std::endl;
      std::cout<<"               refdist = "
               <<NC::fmtg(refdist_2entry)<<std::endl;
      std::cout<<"                  dist = "
               <<NC::fmtg(dist_2entry)<<std::endl;
      if ( has_manualref_d2entry )
        std::cout<<"     expected (manual) = "
                 <<NC::fmtg(tcase.expected_disttoentry.value())<<std::endl;
      if ( testexit ) {
        std::cout<<"   DistToExit:"<<std::endl;
        std::cout<<"               refdist = "
                 <<NC::fmtg(refdist_2exit)<<std::endl;
        std::cout<<"                  dist = "
                 <<NC::fmtg(dist_2exit)<<std::endl;
        if ( has_manualref_d2exit )
          std::cout<<"     expected (manual) = "
                   <<NC::fmtg(tcase.expected_disttoexit.value())<<std::endl;
      }
      std::cout<<std::endl;
    }
    REQUIREFLTEQ_SPHERE( refdist_2entry,
                  tcase.expected_disttoentry.value_or(refdist_2entry));
    REQUIREFLTEQ_SPHERE( refdist_2entry, dist_2entry );
    if ( testexit ) {
      REQUIREFLTEQ_SPHERE( refdist_2exit,
                           tcase.expected_disttoexit.value_or(refdist_2exit));
      REQUIREFLTEQ_SPHERE( refdist_2exit, dist_2exit );
    }
  }
  std::cout<<"Tested a total of "<<test_cases.size()
           <<" Sphere test cases!"<<std::endl;
}

void testSphereCases2()
{
  {
    NCMMC::NeutronBasket b;
    const double x[]  = { -30.0, -30.0, 30.0, 0.0, 0.0, 10.0 };
    const double y[]  = {   0.0,   0.0,  0.0, 0.0, 0.0, 0.0 };
    const double z[]  = {   0.0,   0.0,  0.0, 0.0, 10.0*(1.0-1.0e-14), 0.0 };
    const double ux[] = {   1.0,   0.0,  1.0, 0.0, 0.0, 0.0 };
    const double uy[] = {   0.0,  -1.0,  0.0, 1.0, 0.0, 1.0 };
    const double uz[] = {   0.0,   0.0,  0.0, 0.0, 1.0, 0.0 };
    const double dist_to_entry[] = { 20.0, -1.0, -1.0, 0.0, 0.0, -1.0 };
    constexpr std::size_t n = sizeof(x) / sizeof(*x);
    static_assert(n <= NCMMC::basket_N, "");
    b.nused = 0;
    for( std::size_t i = 0; i < n; ++i ) {
      ++b.nused;
      b.x[i] = x[i];
      b.y[i] = y[i];
      b.z[i] = z[i];
      b.ux[i] = ux[i];
      b.uy[i] = uy[i];
      b.uz[i] = uz[i];
      b.w[i] = 1.0;
      b.ekin[i] = 0.025;
    }
    double buf[NCMMC::basket_N];
    auto geom = NCMMC::createGeometry( "sphere;r=10.0" );
    geom->distToVolumeEntry( b, buf );
    for ( std::size_t i = 0; i < n; ++i )
      REQUIREFLTEQ( buf[i], dist_to_entry[i]);
  }

  {
    NCMMC::NeutronBasket b;
    const double x[]  = { -9.999,  0.0,  5.0, 9.999, 0.0, -10.0, -10.0, -10.0,
                          -10.0 };
    const double y[]  = {   0.0,   0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    const double z[]  = {   0.0,   0.0,  0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0 };
    const double ux[] = {   1.0,   0.0,  1.0, 1.0, 0.0, 0.0, -1.0, 1.0,
                            NC::kInvSqrt2 };
    const double uy[] = {   0.0,  -1.0,  0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                            NC::kInvSqrt2 };
    const double uz[] = {   0.0,   0.0,  0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
    const double dist_to_exit[] = { 19.999, 10.0, 5.0, 0.001, 0.0, 0.0, 0.0,
                                    20.0, 10.0*NC::kSqrt2 };
    constexpr std::size_t n = sizeof(x) / sizeof(*x);
    static_assert(n <= NCMMC::basket_N, "");
    b.nused = 0;
    for( std::size_t i = 0; i < n; ++i ) {
      ++b.nused;
      b.x[i] = x[i];
      b.y[i] = y[i];
      b.z[i] = z[i];
      b.ux[i] = ux[i];
      b.uy[i] = uy[i];
      b.uz[i] = uz[i];
      b.w[i] = 1.0;
      b.ekin[i] = 0.025;
    }
    double buf[NCMMC::basket_N];
    auto geom = NCMMC::createGeometry( "sphere;r=10.0" );
    geom->distToVolumeExit( b, buf );
    for ( std::size_t i = 0; i < n; ++i )
      REQUIREFLTEQ( buf[i], dist_to_exit[i]);
  }
}

int main(int,char**) {

  testSphereCases();
  testSphereCases2();
  testSlabCases();
  testBoxCases();
  std::vector<double> miscvals
    = { -1e199, -1e5, -10, -2.0, -1.0, -0.1, -1e-20, -1e-199, 0.0,
         1e-199, 1e5, 10, 2.0, 1.0, 0.1, 1e20, 1e199 };

  std::vector<double> vals_abslteq2
    = { -2.0, (-2.0+1e-15), (-2.0+1e-14), -1.5, -1.0, -0.5, -1e-20, -1e-199, 0.0,
        2.0, (2.0-1e-15), (2.0-1e-14), 1.5, 1.0, 0.5, 1e-20, 1e-199 };
  {
    auto geom = NCMMC::createGeometry( "slab;dz=2.0" );
    std::size_t ii(0);
    for ( auto x : miscvals ) {
      for ( auto y : miscvals ) {
        for ( auto z : vals_abslteq2 ) {
          if ( (ii++)%50 == 0 )
            std::cout << "x="<<NC::fmt(x)<<", y="<<NC::fmt(y)
                      <<", z="<<NC::fmt(z) <<std::endl;
          REQUIRE( geom->pointIsInside( { 0, 0, z } ) );
          REQUIREFLTEQ( distToVolExit( geom, { x, y, z }, { 0, 0, 1 } ),
                        2.0-z);
          REQUIREFLTEQ( distToVolExit( geom, { x, y, z }, { 0, 0, -1 } ),
                        z + 2.0 );
          REQUIREFLTEQ( distToVolExit( geom, { x, y, z }, { 1, 1, 1 } ),
                        std::sqrt(3)*(2.0-z));
          REQUIREFLTEQ( distToVolExit( geom, { x, y, z }, { 1, 1, -1 } ),
                        std::sqrt(3)*(z+2.0));
          REQUIRE( distToVolEntry( geom, { x, y, z }, { 0, 0, 1 } )
                   == ( z==2.0?-1.0:0.0) );
          REQUIRE( distToVolEntry( geom, { x, y, z }, { -10, 0, 1 } )
                   == ( z==2.0?-1.0:0.0) );
          REQUIRE( distToVolEntry( geom, { x, y, z }, { 0, 0.9, 1 } )
                   == ( z==2.0?-1.0:0.0) );
          REQUIRE( distToVolEntry( geom, { x, y, z }, { 0, 0, -1 } )
                   == ( z==-2.0?-1.0:0.0) );
          REQUIRE( distToVolEntry( geom, { x, y, z }, { 9999, 0, -1 } )
                   == ( z==-2.0?-1.0:0.0) );
          REQUIRE( distToVolEntry( geom, { x, y, z }, { -99, -99999, -1 } )
                   == ( z==-2.0?-1.0:0.0) );
        }
      }
    }
  }

  {
    auto geom = NCMMC::createGeometry( "box;dx=2.0;dy=3.0;dz=5.0" );
    REQUIREFLTEQ( distToVolExit( geom, { 2.0, 3.0, 5.0 }, { 0, 0, 1 } ), 0.0);
    REQUIREFLTEQ( distToVolExit( geom, { 2.0, 3.0, 5.0 }, { 0, 0, -1 } ), 10.0);
    REQUIREFLTEQ( distToVolExit( geom, { 2.0, 3.0, 5.0 }, { 0, 1, 0 } ), 0.0);
    REQUIREFLTEQ( distToVolExit( geom, { 2.0, 3.0, 5.0 }, { 0, -1, 0 } ), 6.0);
    REQUIREFLTEQ( distToVolExit( geom, { 2.0, 3.0, 5.0 }, { 1, 0, 0 } ), 0.0);
    REQUIREFLTEQ( distToVolExit( geom, { 2.0, 3.0, 5.0 }, { -1, 0, 0 } ), 4.0);
    REQUIREFLTEQ( distToVolExit( geom, { 1.0, 1.0, 0.0 }, { 1, 1, 0 } ),
                  std::sqrt(2));

    REQUIREFLTEQ( distToVolEntry( geom, { 2.0+1, 3.0+1, 5.0+1 },
                                  { -1, -1, -1 } ), std::sqrt(3) );
    REQUIREFLTEQ( distToVolEntry( geom, { 2.0+1, 3.0+1, 5.0+1 },
                                  { 1, 1, 1 } ), -1.0 );

    {
      auto d = distToVolEntry( geom, { -4.0, -2, 0 }, { 1, 1, 0 } );
      REQUIREFLTEQ( d, 2*std::sqrt(2) );
    }
    {
      auto d = distToVolEntry( geom, { -3.0, -6, 0 }, { 1, 1, 0 } );
      REQUIREFLTEQ( d, 3*std::sqrt(2) );
    }
    {
      auto d = distToVolEntry( geom, { 0.0, -999, 0 }, { 0, 1, 0 } );
      REQUIREFLTEQ( d, 999.0-3.0 );
    }
    {
      auto d = distToVolEntry( geom, { 2.0, -999, 0 }, { 0, 1, 0 } );
      REQUIREFLTEQ( d, 999.0-3.0 );
    }
  }

  {
    auto geom = NCMMC::createGeometry( "sphere;r=2.0" );
    REQUIREFLTEQ( distToVolExit( geom, { 0, 0, 0 }, { 0, 0, 1 } ), 2.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, 0, 0 }, { 1, 1, 1 } ), 2.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, 0, 0 }, { 1, 1, 0 } ), 2.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, 0, 0 }, { 1, 0, 1 } ), 2.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, 0, 0 }, { 0, 1, 1 } ), 2.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, 0, 0 }, { 1, 0, 0 } ), 2.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, 0, 0 }, { 0, 1, 0 } ), 2.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, 0, -2 }, { 0, 0, 1 } ), 4.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, 0, 2 }, { 0, 0, -1 } ), 4.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, -2, 0 }, { 0, 1, 0 } ), 4.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, 2, 0 }, { 0, -1, 0 } ), 4.0);
    REQUIREFLTEQ( distToVolExit( geom, { -2, 0, 0 }, { 1, 0, 0 } ), 4.0);
    REQUIREFLTEQ( distToVolExit( geom, { 2, 0, 0 }, { -1, 0, 0 } ), 4.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, 0, -2 }, { 0, 0, -1 } ), 0.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, 0, 2 }, { 0, 0, 1 } ), 0.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, -2, 0 }, { 0, -1, 0 } ), 0.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, 2, 0 }, { 0, 1, 0 } ), 0.0);
    REQUIREFLTEQ( distToVolExit( geom, { -2, 0, 0 }, { -1, 0, 0 } ), 0.0);
    REQUIREFLTEQ( distToVolExit( geom, { 2, 0, 0 }, { 1, 0, 0 } ), 0.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, 0, -1 }, { 0, 0, 1 } ), 3.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, 0, 1 }, { 0, 0, 1 } ), 1.0);
    REQUIREFLTEQ( distToVolExit( geom, { 0, 0, -1.9 }, { 0, 0, 1 } ), 3.9);
    REQUIREFLTEQ( distToVolExit( geom, { 0, 0, 1.9 }, { 0, 0, 1 } ), 0.1);
    const double sq3 = std::sqrt(3);
    REQUIREFLTEQ( distToVolExit( geom, { 1/sq3, 1/sq3 , 1/sq3 }, { 1, 1, 1 } ),
                  1.0);

    REQUIREFLTEQ( distToVolEntry( geom, { 0, 0, 0 }, { 0, 0, 1 } ), 0.0);
    REQUIREFLTEQ( distToVolEntry( geom, { 0, 0, 2 }, { 0, 0, 1 } ), -1.0);
    REQUIREFLTEQ( distToVolEntry( geom, { 0, 0, 2 }, { 0, 0, -1 } ), 0.0);
    REQUIREFLTEQ( distToVolEntry( geom, { 0, 2, 0 }, { 0, 0, 1 } ), -1.0);
    REQUIREFLTEQ( distToVolEntry( geom, { 0, 2, 0 }, { 0, 0, -1 } ), -1.0);
    REQUIREFLTEQ( distToVolEntry( geom, { 0, 2, 0 }, { 0, -1, 0 } ), 0.0);
    REQUIREFLTEQ( distToVolEntry( geom, { 0, 2, 0 }, { 0, 1, 0 } ), -1.0);
    REQUIREFLTEQ( distToVolEntry( geom, { 2, 0, 0 }, { 0, 0, 1 } ), -1.0);
    REQUIREFLTEQ( distToVolEntry( geom, { 2, 0, 0 }, { 0, 0, -1 } ), -1.0);
    REQUIREFLTEQ( distToVolEntry( geom, { 2, 0, 0 }, { 1, 0, 0 } ), -1.0);
    REQUIREFLTEQ( distToVolEntry( geom, { 2, 0, 0 }, { -1, 0, 0 } ), 0.0);

    REQUIREFLTEQ( distToVolEntry( geom, { 1/sq3, 1/sq3 , 1/sq3 }, { 1, 1, 1 } ),
                  0.0);
    REQUIREFLTEQ( distToVolEntry( geom, { 3/sq3, 3/sq3 , 3/sq3 },
                                  { -1, -1, -1 } ), 1.0);
    REQUIREFLTEQ( distToVolEntry( geom, { 0, 0, 2 }, { 0, 0, 1 } ), -1.0);
    REQUIREFLTEQ( distToVolEntry( geom, { 0, 0, 2 }, { 0, 0, -1 } ), 0.0);
    REQUIREFLTEQ( distToVolEntry( geom, { 0, 2, 0 }, { 0, 0, 1 } ), -1.0);
    REQUIREFLTEQ( distToVolEntry( geom, { 2, 0, 0 }, { 0, 0, 1 } ), -1.0);

    REQUIREFLTEQ( distToVolEntry( geom, { 2/sq3, 2/sq3 , 2/sq3 }, { 1, 1, 1 } ),
                  -1.0);
    REQUIREFLTEQ( distToVolEntry( geom, { 2/sq3, 2/sq3 , 2/sq3 },
                                  { -1, -1, -1 } ),0.0);

    REQUIREFLTEQ( distToVolEntry( geom, { 2-1e-15, -17.0, 0 },
                                  { 0, 1, 0 } ), 17.0);
    //NB: 1e-15 in the next line fails (returns ~17 not -1.0):
    REQUIREFLTEQ( distToVolEntry( geom, { 2+1e-14, -17.0, 0 },
                                  { 0, 1, 0 } ), -1.0);
  }

  std::cout<<"All tests passed!"<<std::endl;
  return 0;
}
