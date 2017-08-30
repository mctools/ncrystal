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

#include "NCLatticeUtils.hh"
#include "NCMath.hh"
#include "NCVector.hh"

NCrystal::RotMatrix NCrystal::getLatticeRot( double a, double b, double c,
                                             double alpha, double beta, double gamma )
{
  double cg = std::cos(gamma);
  double sg = std::sin(gamma);
  double ca = std::cos(alpha);
  double cb = std::cos(beta);
  nc_assert(sg!=0);
  double m[9] = { a, 0., 0.,
                  b*cg, b*sg, -c*(cb*cg-cos(alpha))/sg,
                  c*cb, c*((ca-cb*cg)/sg), 0. };
  double tmp = c*c-m[6]*m[6]-m[7]*m[7];
  nc_assert(tmp>=0.0);
  m[8]  = sqrt(tmp);
  return RotMatrix(m);
}

NCrystal::RotMatrix NCrystal::getReciprocalLatticeRot( double a, double b, double c,
                                                       double alpha, double beta, double gamma )
{
  RotMatrix m = getLatticeRot(a,b,c,alpha,beta,gamma);
  m.inv();
  m *= (2.0*M_PI);
  return m;
}

void NCrystal::estimateHKLRange( double dcutoff, const NCrystal::RotMatrix& rec_lat,
                                 int& max_h, int& max_k, int& max_l )
{
  //Comment from original location: not missing a single reflection, but the
  //subsequent hkl loop is going to be slightly larger (say a small factor) than
  //it could be.
  nc_assert(dcutoff>0);
  double max_wavelength = 2*M_PI / dcutoff;
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
  return 2*M_PI / ( max_hkl * kmin );
}

