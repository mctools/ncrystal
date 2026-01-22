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

#include "NCrystal/internal/extn_utils/NCExtnBC1974.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NCBC1974 = NCrystal::Extn::BC1974;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    inline double bc1974_calc_cos2th( double sintheta )
    {
      return ncmax( -1.0, 1.0 - 2.0 * ncsquare(sintheta) );
    }

    inline double bc1974_eval_final_formula( double x, double A,
                                             double B, double factor = 2.0 )
    {
      //The original recipe, but including a few safeguards against zero
      //division and returning values above 1.0:
      double t = 1. + B * x;
      //guard against zero division:
      if ( std::fabs( t ) < 1e-20 )
        t = ( t > 0 ? 1e-20 : -1e-20 );
      const double u = 1. + factor * x + ncsquare( x ) * A / t;
      //Simulaneously guard against sqrt of negative numbers, and ensure values
      //are at most 1, by clamping the argument of the sqrt to be >=1:
      return 1. / std::sqrt(std::max<double>(1.0,u));
    }
  }
}

double NCBC1974::y_primary( double x, double sintheta )
{
  const double c = bc1974_calc_cos2th(sintheta);
  const double A = 0.20 + 0.45 * c;
  const double B = 0.22 - 0.12 * ncsquare(0.5 - c);
  return bc1974_eval_final_formula( x, A, B );
}

double NCBC1974::y_scndgauss( double x, double sintheta )
{
  const double c = bc1974_calc_cos2th(sintheta);
  const double A = 0.58 + 0.48  * c + 0.24 * ncsquare(c);
  const double B = 0.02 - 0.025 * c;
  return bc1974_eval_final_formula( x, A, B, 2.12 );
}

double NCBC1974::y_scndlorentz( double x, double sintheta )
{
  const double c = bc1974_calc_cos2th(sintheta);
  const double A = 0.025 + 0.285*c;
  const double B = ( c >= 0.0
                     ? 0.15 - 0.2 * ncsquare( 0.75-c )
                     : -0.45*c );
  return bc1974_eval_final_formula( x, A, B );
}

double NCBC1974::y_scndfresnel( double x, double sintheta )
{
  const double c = bc1974_calc_cos2th(sintheta);
  const double A = 0.48 + 0.6 * c;
  const double B = 0.20 - 0.06 * ncsquare( 0.2 - c );
  return bc1974_eval_final_formula( x, A, B );
}
