////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
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

#include "NCLatticeUtils.hh"
#include "NCMath.hh"
#include "NCVector.hh"

#define NCRYSTAL_EXACT_LATTICEROTS_FOR_SPECIAL_CASES

NCrystal::RotMatrix NCrystal::getLatticeRot( double a, double b, double c,
                                             double alpha, double beta, double gamma )
{
  //sanity check angles in (0,pi) - also confirms that angles are not in degrees by mistake.
  nc_assert_always(alpha<kPi&&alpha>0);
  nc_assert_always(beta<kPi&&beta>0);
  nc_assert_always(gamma<kPi&&gamma>0);
  nc_assert_always(a>0);
  nc_assert_always(b>0);
  nc_assert_always(c>0);
  double cg(2), sg(2), ca(2), cb(2);

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
  if (beta90) { cb = 0.0; }
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

  nc_assert(sg>0);
  double m57 = c*(ca-cb*cg)/sg;
  double m[9] = { a, 0., 0.,
                  b*cg, b*sg, m57,
                  c*cb, m57, 0. };
  m[7]=m[5];
  if (cb||m57) {
    double tmp = c*c-m[6]*m[6]-m57*m57;
    nc_assert(tmp>=0.0);
    m[8]  = std::sqrt(tmp);
  } else {
    //Avoid the sqrt and potential introduction of numerical imprecision:
    m[8] = c;
  }
  return RotMatrix(m);

  //TODO for NC2: Our data-library should include files which does not have
  //alpha=beta=90deg, so we actually validate more of the code above.
}


NCrystal::RotMatrix NCrystal::getReciprocalLatticeRot( double a, double b, double c,
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
    double m[9] = { k2Pi/a, 0., 0.,
                    0., k2Pi/b, 0.,
                    0., 0., k2Pi/c };
    return RotMatrix(m);
  } else if (alpha90 && beta90 && gamma120 ) {
    const double twopidivsqrt3 = 3.6275987284684357011881565152843114645681324961855;
    const double fourpidivsqrt3 = 7.255197456936871402376313030568622929136264992371;
    double m[9] = { k2Pi/a, 0., 0.,
                    twopidivsqrt3/a, fourpidivsqrt3/b, 0.,
                    0., 0., k2Pi/c };
    return RotMatrix(m);
  }
#endif
  RotMatrix m = getLatticeRot(a,b,c,alpha,beta,gamma);
  m.inv();
  m *= k2Pi;
  return m;
}

void NCrystal::estimateHKLRange( double dcutoff, const NCrystal::RotMatrix& rec_lat,
                                 int& max_h, int& max_k, int& max_l )
{
  //Comment from original location: not missing a single reflection, but the
  //subsequent hkl loop is going to be slightly larger (say a small factor) than
  //it could be.
  nc_assert(dcutoff>0);
  double max_wavelength = k2Pi / dcutoff;
  double kh = ( rec_lat * Vector(1,0,0) ).mag();
  double kk = ( rec_lat * Vector(0,1,0) ).mag();
  double kl = ( rec_lat * Vector(0,0,1) ).mag();
  max_h = ceil( max_wavelength / kh);
  max_k = ceil( max_wavelength / kk);
  max_l = ceil( max_wavelength / kl);
}

double NCrystal::estimateDCutoff( int max_hkl, const NCrystal::RotMatrix& rec_lat )
{
  nc_assert(max_hkl>0);
  double kh = ( rec_lat * Vector(1,0,0) ).mag();
  double kk = ( rec_lat * Vector(0,1,0) ).mag();
  double kl = ( rec_lat * Vector(0,0,1) ).mag();
  double kmin = ncmin( kh, ncmin(kk,kl) );
  return k2Pi / ( max_hkl * kmin );
}


void NCrystal::checkAndCompleteLattice( unsigned sg, double a, double& b, double & c )
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
  if ( !(a>0.0) || !(a>0.0) || !(a>0.0) )
    NCRYSTAL_THROW(BadInput,"lattice parameters must be positive numbers");
}

double NCrystal::dspacingFromHKL( int h, int k, int l, const NCrystal::RotMatrix& rec_lat )
{
  if (h==0&&k==0&&l==0)
    NCRYSTAL_THROW(BadInput,"Can not calculate d-spacing for hkl=000");
  double ksq = ( rec_lat * Vector(h,k,l) ).mag2();
  if (!(ksq>0.0))
    NCRYSTAL_THROW(CalcError,"Created invalid k-vector in d-spacing calculations (bad lattice rotation provided?)");
  return k2Pi / std::sqrt(ksq);
}
