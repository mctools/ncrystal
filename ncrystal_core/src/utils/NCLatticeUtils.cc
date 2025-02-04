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

#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCString.hh"
namespace NC = NCrystal;

#define NCRYSTAL_EXACT_LATTICEROTS_FOR_SPECIAL_CASES

NC::RotMatrix NC::getLatticeRot( double a, double b, double c,
                                 double alpha, double beta, double gamma )
{
  //sanity check angles in (0,pi) - also confirms that angles are not in degrees by mistake.
  nc_assert_always(alpha<kPi&&alpha>0);
  nc_assert_always(beta<kPi&&beta>0);
  nc_assert_always(gamma<kPi&&gamma>0);
  nc_assert_always(a>0);
  nc_assert_always(b>0);
  nc_assert_always(c>0);
  double cg(2), sg(2), ca(2), cb(2), sb(2);

#ifdef NCRYSTAL_EXACT_LATTICEROTS_FOR_SPECIAL_CASES
  //Explicitly look for angles 90 and 120 that are very common in crystals, so
  //we can avoid some numerical imprecisions in those cases and instead use
  //exact values for the sines and cosines. This is needed because
  //e.g. std::cos(kPiHalf)=6.1e-17 rather than 0.0 exactly, resulting in
  //incorrect tiny non-zero matrix elements where there should be zero. This
  //propagates further, for instance resulting in funny crystal normals like {
  //0, 1, -6.12323e-17 } instead of the correct { 0, 1, 0 }. This is mostly a
  //cosmetic issue of course.
  const double tol = 1e-14;
  bool alpha90 = ncabs(alpha-kPiHalf)<tol;
  bool beta90 = ncabs(beta-kPiHalf)<tol;
  bool gamma90 = ncabs(gamma-kPiHalf)<tol;
  bool gamma120 = ncabs(gamma-kPi*(2/3.))<tol;

  if (gamma90) { cg = 0.0; sg = 1.0; }
  else if (gamma120) { cg = -0.5; sg = 0.86602540378443864676372317075293618347140262690519; }
  if (alpha90) { ca = 0.0; }
  if (beta90) { cb = 0.0; sb = 1.0; }
  if (!ca&&!cb&&!cg) {
    double m[9] = { a, 0., 0.,
                    0., b, 0.,
                    0., 0., c };
    return RotMatrix(m);
  }
#endif

  if (cg==2) cg = std::cos(gamma);
  if (sg==2) sg = std::sin(gamma);
  if (ca==2) ca = std::cos(alpha);
  if (cb==2) cb = std::cos(beta);
  if (sb==2) sb = std::sin(beta);

  nc_assert(sg>0);
  const double m57 = c*(ca-cb*cg)/sg;
  double m[9] = { a, 0., 0.,
                  b*cg, b*sg, m57,
                  c*cb, m57, 0. };//last entry changed below
  if ( !cb && !m57 ) {
    //Avoid the sqrt and potential introduction of numerical imprecision:
    m[8] = c;
  } else if ( cb && !m57 ) {
    //Avoid the sqrt and potential introduction of numerical imprecision:
    m[8] = c*sb;
  } else {
    //most general case
    double tmp = sb*sb - ncsquare( (ca-cb*cg)/sg );
    nc_assert( tmp >=0.0 );
    m[8]  = c * std::sqrt(tmp);
  }
  return RotMatrix(m);
}


NC::RotMatrix NC::getReciprocalLatticeRot( double a, double b, double c,
                                           double alpha, double beta, double gamma )
{
  //sanity check angles in (0,pi) - also confirms that angles are not in degrees by mistake.
  nc_assert_always(alpha<kPi&&alpha>0);
  nc_assert_always(beta<kPi&&beta>0);
  nc_assert_always(gamma<kPi&&gamma>0);
  nc_assert_always(a>0);
  nc_assert_always(b>0);
  nc_assert_always(c>0);
#ifdef NCRYSTAL_EXACT_LATTICEROTS_FOR_SPECIAL_CASES
  //Like in getLatticeRot we explicitly treat angles 90 and 120.
  const double tol = 1e-14;
  bool alpha90 = ncabs(alpha-kPiHalf)<tol;
  bool beta90 = ncabs(beta-kPiHalf)<tol;
  bool gamma90 = ncabs(gamma-kPiHalf)<tol;
  bool gamma120 = ncabs(gamma-kPi*(2/3.))<tol;
  if ( alpha90 && beta90 && gamma90 ) {
    //Usually this is cubic, tetragonal, orthorombic
    std::array<double,9> m{ k2Pi/a, 0., 0.,
                            0., k2Pi/b, 0.,
                            0., 0., k2Pi/c };
    return RotMatrix(m);
  } else if (alpha90 && beta90 && gamma120 ) {
    //Usually this is trigonal/hexagonal
    const double twopidivsqrt3 = 3.6275987284684357011881565152843114645681324961855;
    const double fourpidivsqrt3 = 7.255197456936871402376313030568622929136264992371;
    std::array<double,9> m{ k2Pi/a, 0., 0.,
                            twopidivsqrt3/a, fourpidivsqrt3/b, 0.,
                            0., 0., k2Pi/c };
    return RotMatrix(m);
  } else if ( alpha90 && gamma90 ) {
    //Usually this is monoclinic (with the b-axis unique setting)

    //Todo: alpha=gamma=90 for monoclinic corresponds to the b-axis unique
    //      setting, which is the most widespread setting. We might consider
    //      also supporting the a-axis unique setting (beta=gamma=90) and the
    //      c-axis unique setting (alpha=beta=90).

    const double sb = std::sin(beta);
    const double tanb = std::tan(beta);
    nc_assert(tanb!=0.0);
    const double cotb = 1.0 / tanb;//slightly more accurate than cos(b)/sin(b)
    std::array<double,9> m{ k2Pi/a, 0., 0.,
                             0., k2Pi/b, 0.,
                            -cotb*k2Pi/a, 0., k2Pi/(c*sb) };
    return RotMatrix(m);
  }
#endif
  //general case (usually triclinic and rhombohedral ends here)
  RotMatrix m = getLatticeRot(a,b,c,alpha,beta,gamma);
  m.inv();
  m *= k2Pi;
  return m;
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    static constexpr std::size_t npts_sphere = 1000;
    const std::array<Vector,npts_sphere>& manyPtsAllOverUnitSphere() {
      //Fibonacci method of approximately distributing N points over a sphere.
      static std::array<Vector,npts_sphere> s_v;
      static std::mutex s_mtx;
      NCRYSTAL_LOCK_GUARD(s_mtx);
      if (s_v[0].mag2() != 0.0)
        return s_v;
      constexpr double golden_angle = (3.0-constexpr_sqrt(5.0))*kPi;
      const double dz = (2.0-2.0/npts_sphere)/(npts_sphere-1.0);
      const double z0 = -1.0 + 1.0/npts_sphere;
      for ( auto i : ncrange(npts_sphere) ) {
        const double z = z0+i*dz;
        const double theta = golden_angle * i;
        double cosTh, sinTh;
        sincos( theta, cosTh, sinTh );
        double r = std::sqrt( 1.0 - ncsquare(z) );
        s_v[i] = Vector( r*cosTh, r*sinTh, z );
      }
      return s_v;
    }
  }
}

NC::MaxHKL NC::estimateHKLRange( double dcutoff,
                                 double a, double b, double c,
                                 double alpha, double beta, double gamma )
{
  nc_assert(dcutoff>0);
  const double invdcutoff = 1.0/dcutoff;
  nc_assert_always(alpha<kPi&&alpha>0);
  nc_assert_always(beta<kPi&&beta>0);
  nc_assert_always(gamma<kPi&&gamma>0);
  nc_assert_always(a>0);
  nc_assert_always(b>0);
  nc_assert_always(c>0);
  const double tol = 1e-14;
  bool alpha90 = ncabs(alpha-kPiHalf)<tol;
  bool beta90 = ncabs(beta-kPiHalf)<tol;
  bool gamma90 = ncabs(gamma-kPiHalf)<tol;

  auto floorhkl = []( double x ) -> int
  {
    double xc = std::floor( x * (1.0+1e-12) );
    if ( xc <= 1.0 )
      return 1;
    if ( xc >= static_cast<double>(std::numeric_limits<int>::max()) )
      return std::numeric_limits<int>::max();
    return static_cast<int>(xc+0.5);
  };

  if ( alpha90 && beta90 && gamma90 ) {
    //Handle cubic, tetragonal, and orthorombic easily (could in principle also
    //deal with some highly unusual degenerate mono/tri-clinic cases):
    return { floorhkl( a*invdcutoff ), floorhkl( b*invdcutoff ), floorhkl( c*invdcutoff ) };
  }

  //Continue with generic code. TODO: Add more special cases with analytical
  //code for certain additional symmetries (monoclinic and/or trigonal/hexagonal
  //would be nice).

  //Brute force check of the extreme reach along the (h,k,l) axes of mapping a
  //unit sphere in wavevector space onto (h,k,l) space.

  //The relevant mapping matrix is the following multiplied by 1/2pi (we leave
  //the 1/2pi for later):
  auto trf = getLatticeRot( a, b, c, alpha, beta, gamma );

  double max_reach_h(0.0), max_reach_k(0.0), max_reach_l(0.0);
  for ( auto& pt : manyPtsAllOverUnitSphere() ) {
    auto hkl_pt = trf*pt;
    max_reach_h = ncmax( max_reach_h, ncabs(hkl_pt.x()) );
    max_reach_k = ncmax( max_reach_k, ncabs(hkl_pt.y()) );
    max_reach_l = ncmax( max_reach_l, ncabs(hkl_pt.z()) );
  }

  //Add safety factor since manyPtsAllOverUnitSphere() does not cover ALL pts on
  //sphere:
  max_reach_h *= 1.05;
  max_reach_k *= 1.05;
  max_reach_l *= 1.05;

  //Now, we were not actually interested in the extreme reach of a unit sphere,
  //but rather a sphere of radius 2pi/dcutoff (which is the longest wavevectors
  //allowed by the dcut). Since the transformations are linear, we can simply
  //scale the extreme reach. However, we must remember to multiply with 1/2pi
  //which we left out above, so in the end we just have to multiply with
  //1/dcutoff:
  return { floorhkl( max_reach_h * invdcutoff ),
           floorhkl( max_reach_k * invdcutoff ),
           floorhkl( max_reach_l * invdcutoff ) };
}

void NC::checkAndCompleteLatticeAngles( unsigned sg, double& alpha, double& beta, double& gamma )
{
  if (sg>230)
    NCRYSTAL_THROW2(BadInput,"invalid spacegroup number ("<<sg<<")");
  if ( sg < 1 )
    return;
  auto cs = crystalSystem( sg );
  switch (cs) {
  case Orthorhombic:
  case Tetragonal:
  case Cubic:
    //all 90
    if ( ( alpha>0 && alpha!=90 ) || ( beta>0 && beta!=90 ) || ( gamma>0 && gamma!=90 ) )
      NCRYSTAL_THROW2(BadInput,"Spacegroup ("<<sg<<") requires alpha=beta=gamma=90");
    alpha = beta = gamma = 90;
    return;
  case Trigonal:
  case Hexagonal:
    if ( ( alpha>0 && alpha!=90 ) || ( beta>0 && beta!=90 ) || ( gamma>120 && gamma!=120 ) )
      NCRYSTAL_THROW2(BadInput,"Spacegroup ("<<sg<<") requires alpha=beta=90 and gamma=120");
    alpha = beta = 90;
    gamma = 120;
    return;
  case Monoclinic:
  case Triclinic:
    //Although we might be able to do something more specific, we will for now.
    //simply require all three angles to be set with values <180
    if ( !( alpha>0 && alpha<180 ) || !( beta>0 && beta<180 ) || !( gamma>0 && gamma<180 ) )
      NCRYSTAL_THROW2(BadInput,"Spacegroup ("<<sg<<") requires all three angles to be set (and to values < 180).");
    return;
  };
}


void NC::checkAndCompleteLattice( unsigned sg, double a, double& b, double & c )
{
  if (sg>230)
    NCRYSTAL_THROW2(BadInput,"invalid spacegroup number ("<<sg<<")");
  if (sg>0) {
    //lattice_b and lattice_c does not always need to be provided if identical
    //to lattice_a by spacegroup symmetry:
    bool b_is_a = false;
    bool c_is_a = false;
    if (sg >= 75 && sg <= 194 )
      b_is_a = true;
    else if ( sg >= 195 && sg <= 230 )
      c_is_a = b_is_a = true;

    if (b_is_a) {
      if (!b) {
        b = a;
      } else if (!(b==a)) {
        NCRYSTAL_THROW2(BadInput,"lattice lengths a and b must be equal for spacegroup "<<sg);
      }
    }
    if (c_is_a) {
      if (!c) {
        c = a;
      } else if (!(c==a)) {
        NCRYSTAL_THROW2(BadInput,"lattice lengths a and c must be equal for spacegroup "<<sg);
      }
    }
  }
  if ( !(a>0.0) || !(b>0.0) || !(c>0.0) )
    NCRYSTAL_THROW(BadInput,"lattice parameters must be positive numbers");
}

double NC::dspacingFromHKL( int h, int k, int l, const RotMatrix& rec_lat )
{
  if (h==0&&k==0&&l==0)
    NCRYSTAL_THROW(BadInput,"Can not calculate d-spacing for hkl=000");
  double ksq = ( rec_lat * Vector(h,k,l) ).mag2();
  if (!(ksq>0.0))
    NCRYSTAL_THROW(CalcError,"Created invalid k-vector in d-spacing calculations (bad lattice rotation provided?)");
  return k2Pi / std::sqrt(ksq);
}

  // Validate orientation parameters (check against null-vectors, parallel
  // vectors, consistency of angles in the two frames given dirtol,
  // etc.). Throws BadInput exception in case of issues.
namespace NCRYSTAL_NAMESPACE {
  namespace  {
    constexpr const char * vectName( const CrystalAxis& ) { return "CrystalAxis"; }
    constexpr const char * vectName( const HKLPoint& ) { return "HKLPoint"; }
    constexpr const char * vectName( const LabAxis& ) { return "LabAxis"; }
    static constexpr double orient_paralleltol = 1e-6;
    template<class TVect>
    void checkNotNull( const TVect& v, const char * designator = nullptr)
    {
      if (!v.template as<Vector>().mag2())
        NCRYSTAL_THROW2(BadInput,"Specified "
                        <<(designator?designator:"")
                        <<(designator?" ":"")
                        <<vectName(v)<<" is a null-vector.");
    }
    void baseCheckOD( const OrientDir& odir, const char * designator = nullptr )
    {
      checkNotNull(odir.lab,designator);
      const auto& cc = odir.crystal;
      if ( cc.has_value<HKLPoint>() )
        checkNotNull(cc.get<HKLPoint>(),designator);
      else if ( cc.has_value<CrystalAxis>() )
        checkNotNull(cc.get<CrystalAxis>(),designator);
      else
        NCRYSTAL_THROW2(BadInput,"Invalid crystal direction object (empty crystal direction)");
    }
  }

}

void NC::precheckLatticeOrientTol( double dirtol )
{
  if ( !( dirtol > 0.0 ) || dirtol > kPi )
    NCRYSTAL_THROW(BadInput, "Directional tolerance must be in interval (0.0,pi]");
}

void NC::precheckLatticeOrientDir( const OrientDir& odir )
{
  baseCheckOD(odir);
}

void NC::verifyLatticeOrientDef( const LabAxis& l1, const CrystalAxis& c1,
                                 const LabAxis& l2, const CrystalAxis& c2,
                                 double dirtol )
{
  checkNotNull(l1,"primary");
  checkNotNull(c1,"primary");
  checkNotNull(l2,"secondary");
  checkNotNull(c2,"secondary");
  precheckLatticeOrientTol(dirtol);

  if ( l1.as<Vector>().isParallel( l2.as<Vector>(), orient_paralleltol ) )
    NCRYSTAL_THROW(BadInput, "Specified primary and secondary lab directions are parallel");

  if ( c1.as<Vector>().isParallel( c2.as<Vector>(), orient_paralleltol ) )
    NCRYSTAL_THROW(BadInput, "Specified primary and secondary crystal directions are parallel");

  const double anglec = c1.as<Vector>().angle(c2.as<Vector>());
  const double anglel = l1.as<Vector>().angle(l2.as<Vector>());
  if ( ncabs( anglec - anglel ) > dirtol ) {
    NCRYSTAL_THROW2(BadInput,"Chosen orientation defining directions in the lab frame are "<<dbl2shortstr(anglel*kToDeg)
                    <<" deg apart, while the chosen directions in the crystal frame"
                    " are "<<dbl2shortstr(anglec*kToDeg)<<" deg apart. This is not within the specified"
                    " tolerance of "<<dbl2shortstr(dirtol)<<" rad. = "<<dbl2shortstr(dirtol*kToDeg)<<" deg.");
  }

}

void NC::precheckLatticeOrientDef( const OrientDir& dir1, const OrientDir& dir2, double dirtol )
{
  if ( dir1.crystal.has_value<CrystalAxis>() && dir2.crystal.has_value<CrystalAxis>() ) {
    //We can immediately do a full check:
    verifyLatticeOrientDef( dir1.lab, dir1.crystal.get<CrystalAxis>(),
                            dir2.lab, dir2.crystal.get<CrystalAxis>(),
                            dirtol );
    return;
  }
  //Partially check as much as we can:
  baseCheckOD(dir1,"primary");
  baseCheckOD(dir2,"secondary");
  precheckLatticeOrientTol(dirtol);
  if ( dir1.lab.as<Vector>().isParallel( dir2.lab.as<Vector>(), orient_paralleltol ) )
    NCRYSTAL_THROW(BadInput, "Specified primary and secondary lab directions are parallel");
  if ( dir1.crystal.has_value<HKLPoint>() && dir2.crystal.has_value<HKLPoint>() ) {
    if ( dir1.crystal.get<HKLPoint>().as<Vector>().isParallel( dir2.crystal.get<HKLPoint>().as<Vector>(), orient_paralleltol ) )
      NCRYSTAL_THROW(BadInput, "Specified primary and secondary crystal directions (hkl points) are parallel");
  }
}

NC::RotMatrix NC::verifyLatticeOrientDefAndConstructCrystalRotation( const OrientDir& dir1,
                                                                     const OrientDir& dir2,
                                                                     double dirtol,
                                                                     const RotMatrix& reci_lattice )
{
  precheckLatticeOrientDef( dir1, dir2, dirtol );
  nc_assert(!dir1.crystal.empty());
  nc_assert(!dir2.crystal.empty());

  Vector c1 = ( dir1.crystal.has_value<CrystalAxis>()
                ? dir1.crystal.get<CrystalAxis>().as<Vector>()
                : ( reci_lattice * dir1.crystal.get<HKLPoint>().as<Vector>().unit() ) ).unit();
  Vector c2 = ( dir2.crystal.has_value<CrystalAxis>()
                ? dir2.crystal.get<CrystalAxis>().as<Vector>()
                : ( reci_lattice * dir2.crystal.get<HKLPoint>().as<Vector>().unit() ) ).unit();
  Vector l1 = dir1.lab.as<Vector>().unit();
  Vector l2 = dir2.lab.as<Vector>().unit();

  verifyLatticeOrientDef( l1.as<LabAxis>(), c1.as<CrystalAxis>(),
                          l2.as<LabAxis>(), c2.as<CrystalAxis>(),
                          dirtol );

  //We are within the tolerance and all is OK (nothing null, parallel, etc.), so
  //we now ensure exact angle(c1,c2)==angle(l1,l2) by removing components of the
  //secondary direction parallel to the primary direction:

  c2 -= c1 * c2.dot(c1);
  l2 -= l1 * l2.dot(l1);
  c2.normalise();
  l2.normalise();
  return RotMatrix(l1,c1,l2,c2);
}

NC::CrystalSystem NC::crystalSystem( int sg )
{
  if ( sg < 1 || sg > 230 )
    NCRYSTAL_THROW(BadInput,"Space group number is not in the range 1 to 230");
  if ( sg >= 143 ) {
    if ( sg >= 195 )
      return CrystalSystem::Cubic;
    return sg >= 168 ? CrystalSystem::Hexagonal : CrystalSystem::Trigonal;
  }
  if ( sg >= 75 )
    return CrystalSystem::Tetragonal;
  if ( sg >= 16 )
    return CrystalSystem::Orthorhombic;
  return sg >= 3 ? CrystalSystem::Monoclinic : CrystalSystem::Triclinic;
}

std::tuple<int,int,int> NC::normalAndDSpacingToHKLIndex( const RotMatrix& lattice_rot,
                                                         double dspacing,
                                                         const Vector& normal )
{
  nc_assert( dspacing > 0.0 );
  auto hkl = lattice_rot * normal;
  hkl /= dspacing;
  if ( hkl < -hkl )
    hkl = -hkl;
  Vector hkl_rounded( std::round(hkl[0]), std::round(hkl[1]), std::round(hkl[2]) );
  if ( ( hkl - hkl_rounded ).mag2() > 1e-10 )
    NCRYSTAL_THROW(CalcError,"HKL point estimated from dspacing+normal is not approximately integral.");
  return { int(hkl_rounded[0]), int(hkl_rounded[1]), int(hkl_rounded[2]) };
}
