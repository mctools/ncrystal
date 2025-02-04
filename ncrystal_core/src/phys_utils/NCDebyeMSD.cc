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

#include "NCrystal/internal/phys_utils/NCDebyeMSD.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRomberg.hh"

namespace NC = NCrystal;

double NC::debyeIsotropicMSD( DebyeTemperature dt, Temperature t, AtomMass am )
{
  dt.validate();
  //Don't do this since we want to allow t=0: t.validate();
  am.validate();
  nc_assert_always(dt.get()>0.0&&dt.get()<1e5);
  nc_assert_always(t.get()>=0.0&&t.get()<=Temperature::allowed_range.second);
  nc_assert_always(am.get()>=1.007&&am.get()<500);
  return calcDebyeMSDScale( dt, am )*calcDebyeMSDShape(t.get()/dt.get());
}

double NC::calcDebyeMSDScale( DebyeTemperature dt, AtomMass am )
{
  dt.validate();
  am.validate();
  nc_assert_always(dt.get()>0.0);
  nc_assert_always(am.get()>=1.007&&am.get()<500);
  const double kk = 3.0*constant_hbar*constant_hbar*constant_c*constant_c / ( constant_dalton2eVc2*constant_boltzmann );
  return kk/(am.get()*dt.get());
}

namespace NCRYSTAL_NAMESPACE {
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
      return x / std::expm1(x);
    }
    virtual bool accept(unsigned lvl, double prev_estimate, double estimate,double,double) const
    {
      return lvl>11 || (lvl > 7 && ncabs(estimate-prev_estimate)<1e-10);
    }
  };
}

double NC::calcDebyeMSDShape( double x )
{
  nc_assert_always(x>=0.0);
  if (x<1e-50)
    return 0.25;
  DebyeMSDShapeIntegral integral;
  return 0.25 + x * x * integral.integrate( 0.0, 1.0 / x );
}

NC::DebyeTemperature NC::debyeTempFromIsotropicMSD( double msd, Temperature t, AtomMass am )
{
  //For stability and efficiency, first perform a brute-force search of limits
  //before unleashing the generic root finding algorithm.

  auto calcMSD = [&t,&am](double debye_temp) { return debyeIsotropicMSD(DebyeTemperature{debye_temp},t,am); };
  double dt_low(200.0), dt_high(300.0);
  while ( calcMSD(dt_low) <= msd ) {
    dt_high = dt_low;
    dt_low /= 1.5;
    if (dt_low < 1e-6)
      NCRYSTAL_THROW(CalcError,"Can not determine Debye temperature from isotropic MSD (too loosely bound atoms?)");
  }
  while ( calcMSD(dt_high) >= msd ) {
    dt_low = dt_high;
    dt_high *= 1.5;
    if (dt_low > 0.999e6)
      NCRYSTAL_THROW(CalcError,"Can not determine Debye temperature from isotropic MSD (too tightly bound atoms?)");
  }

  return DebyeTemperature { findRoot2( [&calcMSD,msd](double dt) { return calcMSD(dt)-msd; }, dt_low,dt_high,1e-7) };
}
