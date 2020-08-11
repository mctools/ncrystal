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

#include "NCDebyeMSD.hh"
#include "NCMath.hh"
#include "NCRomberg.hh"

double NCrystal::debyeIsotropicMSD(double debye_temperature, double temperature, double atomic_mass)
{
  nc_assert_always(debye_temperature>0.0&&debye_temperature<1e5);
  nc_assert_always(temperature>=0.0&&temperature<1e5);
  nc_assert_always(atomic_mass>=1.007&&atomic_mass<500);
  return calcDebyeMSDScale( debye_temperature, atomic_mass )*calcDebyeMSDShape(temperature/debye_temperature);
}

double NCrystal::calcDebyeMSDScale( double debye_temperature, double atomic_mass )
{
  nc_assert_always(debye_temperature>0.0);
  nc_assert_always(atomic_mass>=1.007&&atomic_mass<500);
  const double kk = 3.0*constant_hbar*constant_hbar*constant_c*constant_c / ( constant_dalton2eVc2*constant_boltzmann );
  return kk/(atomic_mass*debye_temperature);
}

namespace NCrystal {
  class DebyeMSDShapeIntegral : public Romberg {
  public:
    //Integrate function f(x)=x/(exp(x)-1)
    DebyeMSDShapeIntegral() {}
    virtual ~DebyeMSDShapeIntegral() {}
    virtual double evalFunc(double x) const {
      if (ncabs(x)<1e-4) {
        //evaluate via Taylor expansion for numerical stability
        double x2=x*x;
        return 1-x*0.5+x2*0.08333333333333333333333333333333-x2*x2*0.00138888888888888888888888888888888889;
      }
#if __cplusplus >= 201103L
      return x / std::expm1(x);
#else
      return x / (std::exp(x)-1.0);
#endif
    }
    virtual bool accept(unsigned lvl, double prev_estimate, double estimate,double,double) const
    {
      return lvl>11 || (lvl > 7 && ncabs(estimate-prev_estimate)<1e-10);
    }
  };
}

double NCrystal::calcDebyeMSDShape( double x )
{
  nc_assert_always(x>=0.0);
  if (x<1e-50)
    return 0.25;
  DebyeMSDShapeIntegral integral;
  return 0.25 + x * x * integral.integrate( 0.0, 1.0 / x );
}

double NCrystal::debyeTempFromIsotropicMSD(double msd, double temperature, double atomic_mass)
{
  return findRoot2([msd,temperature,atomic_mass](double dt){ return debyeIsotropicMSD(dt,temperature,atomic_mass)-msd;},
                   0.1,0.999e5,1e-7);
}
