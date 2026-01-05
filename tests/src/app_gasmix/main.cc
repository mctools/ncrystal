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

#include "NCrystal/internal/gasmix/NCGasMixUtils.hh"
#include "NCrystal/internal/atomdb/NCAtomDB.hh"
#include <iostream>
namespace NC = NCrystal;
namespace NCGM = NCrystal::GasMix;

void test_satvappressure_water()
{

  // Values < 0C from https://journals.ametsoc.org/view/journals/apme/57/6/jamc-d-17-0334.1.xml
  // Values between 0 and 100C are from wikipedia
  // Values above 100C from:
  // https://www.engineeringtoolbox.com/water-vapor-saturation-pressure-d_599.html

  const double ref_T_celcius[] = { -100.0,-80.0,-60.0,-40.0,-20.0,
                                   0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100.,
                                   110.,120.,130.,140.,150.,160.,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,370.
  };
  const double ref_P_kPa[] = { 0.0014049e-3,0.054773e-3,1.0813e-3,12.8412e-3,103.239e-3,
                                 0.6113,0.8726,1.2281,1.7056,2.3388,3.1690,4.2455,5.6267,7.3814,9.5898,12.3440,
                                 15.7520,19.9320,25.0220,31.1760,38.5630,47.3730,57.8150,70.1170,84.5290,101.3200,
                                 143.38,198.67,270.28,361.54,476.16,618.23,1002.8,1554.9,2319.6,3346.9,4692.3,6416.6,8587.9,11284,14601,18666,21044
  };

  unsigned iref(0);
  for ( auto tempval_C : ref_T_celcius ) {
    auto tempval = NC::Temperature{ tempval_C + 273.15 };
    auto ref_pressureval = NC::Pressure{ ref_P_kPa[iref++]*1e3 };
    auto calc_pressureval = NCGM::saturatedVapourPressureOfWater( tempval );
    double reldiff = std::abs(calc_pressureval.dbl()/ref_pressureval.dbl() - 1.0);
    std::cout<<"Water vapour P @ "<<tempval<<" : "<<calc_pressureval<<" (ref. "<<ref_pressureval<<"). Reldiff = "<<reldiff<<std::endl;

    if ( tempval_C >= -100.0 && tempval_C <= 100.0 ) {
      nc_assert_always(reldiff<1e-3);//0.1% precision in core interval of -100C..100C
    }
    if ( tempval_C >= -100.0 && tempval_C <= 200.0 ) {//1.12% precision in -200C..200C
      nc_assert_always(reldiff<1.12e-2);
    }
    if ( tempval_C >= 200.0 && tempval_C <= 300.0 ) {//7% precision in 200C..300C
      nc_assert_always(reldiff<0.07);
    }
    if ( tempval_C >= 300.0 && tempval_C <= 400.0 ) {//15% precision in 300C..400C (well, last datapoint is 370C)
      nc_assert_always(reldiff<0.15);
    }
  }
}

void test_parse()
{
  auto test = []( const char * str )
  {
    std::cout<<"Parsing: >>>"<<str<<"<<<"<<std::endl;
    auto req = NCGM::requestFromString(str);
    auto str2 = requestToString(req);
    std::cout << "  \\-> reencoded: >>>"<<str2<<"<<<"<<std::endl;
    std::cout << "  \\-> gives: "<<NCGM::analyseGasMixRequest(req) << std::endl;
  };
  auto testBad = []( const char * str )
  {
    std::cout<<"Parsing (expects BadInput failure): >>>"<<str<<"<<<"<<std::endl;
    bool goterr = false;
    try {
      auto req = NCGM::requestFromString(str);
      auto str2 = requestToString(req);
      std::cout << "  \\-> reencoded: >>>"<<str2<<"<<<"<<std::endl;
      auto res = NCGM::analyseGasMixRequest(req);
      std::cout << "  \\-> gives: "<< res << std::endl;
    } catch ( NC::Error::BadInput& err ) {
      std::cout<<"  \\-> got expected BadInput error: >>>"<<err.what()<<"<<<"<<std::endl;
      goterr = true;
    }
    if (!goterr)
      NCRYSTAL_THROW(LogicError,"Did not end in error as expected!");
  };


  testBad("");
  testBad("1.0*CO2");
  test("1.0xCO2");
  test("CO2");
  test("OCO");
  test("1.0xCH4");
  test("1.0xCO2////");
  testBad("1.0xCO2\"");
  testBad("1.0xCO2'");
  testBad("1.0xCO2$");
  testBad("1.0xCO2@");
  testBad("1.0xCO2\\");

  testBad("1.0xCO2\t");
  testBad("1.0xCO2\n");
  test("1.0xCO2/10bar");
  testBad("1.0xCO2/1e20bar");
  testBad("1.0xCO2/0.0bar");
  testBad("1.0xCO2/-1.0bar");
  testBad("1.0xCO2/1.0bar/1000Pa");
  testBad("1.0xCO2/10t-4bar");
  test("0.3xCO2+0.7xAr/2atm");
  test("0.3xCO2+0.7xH2O/H_is_D");
  test("air");
  test("0.3xCO2+0.7xAr/0.02relhumidity");
  testBad("0.3xCO2+0.8xAr");
  test("0.3xCO2+0.7000000001xAr");
  test("massfractions/0.3xCO2+0.7xAr");
  testBad("massfractions/0.3xCO2+0.7xAr/molarfractions");
  testBad("massfractions/0.3xCO2+0.7xAr/massfractions");
  testBad("1.0xCO2/1.0bar/2gcm3");
  testBad("1.0xCO2/1kgm3/2gcm3");
  testBad("1.0xCO2/1perAa3");
  test("1.0xCO2/1kgm3");
  testBad("1.0xCO2/0.0kgm3");
  testBad("1.0xCO2/1e20kgm3");
  test("1.0xCO2/0.0relhumidity");
  test("1.0xCO2/1.0relhumidity");
  test("1.0xCO2/0.01relhumidity");
  testBad("1.0xCO2/1.0000000000001relhumidity");
  testBad("1.0xCO2/-1e-20relhumidity");
  testBad("0.3xCO2+0.7xAr/0.02relhumidity/0.02relhumidity");
  testBad("1.0xCO2/1kgm3/1.0xCO2");
  testBad("1.0xCO2/1kgm3/1.0xAr");
  test("1.0xCO2/300K");
  test("1.0xCO2/293.15K");
  test("1.0xCO2/20C");
  test("1.0xCO2/-270C");
  testBad("1.0xCO2/-270CK");
  test("1.0xCO2/0.001K");
  testBad("1.0xCO2/0.0K");
  testBad("1.0xCO2/1e20K");
  test("1.0xCO2/1000000C");
  testBad("1.0xCO2/wuhu");
  testBad("1.0xCO2/270C/10K");
  testBad("1.0xSiH4/Si 28.09u 4.1491fm 0.004b 0.171b");
  testBad("1.0xSiH4/Si;28.09u;4.1491fm;0.004b;0.171b");
  test("1.0xSiH4/Si_28.09u_4.1491fm_0.004b_0.171b/H_is_D");
  testBad("1.0xSiH4/Si_28.09u_4.1491fm_0.004b__0.171b/H_is_D");
  test("1.0xSiH4/ Si_28.09u_4.1491fm_0.004b_0.171b  /  /  H_is_D");
  testBad("1.0xSiH4/ Si_28.09u_4.1491fm_0.004b_0.171b  /  /  H_isss_D");
  testBad("1.0xSiH4/ Si_28.09u_4.1491fm_0.004b_0.171b  /  /  H_is_G");
  testBad("1kgm3");
  test("1.0xHe/1kgm3");
  test("//1.0xHe/1kgm3");
  test("air/1kgm3");
  test("air/2atm");
  test("air/2atm/0.01relhumidity");
  test("air/2atm/0.01relhumidity/molarfractions");
  testBad("air/2atm/0.01relhumidity/massfractions");
  testBad("0.3xCO2+1.01xAr");
  testBad("0.3xCO2+-0.11xAr");
  test("1.0xCO2+0.xAr");
  test("1.0xCO2/10bar");
  test("1.0xCO2/2atm");
  test("1.0xCO2/20F");
  test("1.0xCO2/0F");
  test("1.0xCO2/100F");
  test("1.0xCO2/-1.234F");
  testBad("0.2xCO2+0.8xAr/0.5relhumidity/H_is_D");
}

void test_results()
{
  auto test = []( const char * str )
  {
    std::cout<<"Testing: >>>"<<str<<"<<<"<<std::endl;
    auto req = NCGM::requestFromString(str);
    auto str2 = requestToString(req);
    std::cout << "  \\-> reencoded: >>>"<<str2<<"<<<"<<std::endl;
    auto res = NCGM::analyseGasMixRequest(req);
    std::cout << "  \\-> result: "<<res << std::endl;
    std::cout << "  \\-> flat composition: "<< flattenComponentList( res ) << std::endl;
    std::cout << "  \\-> atomdb: "<<std::endl;
    for ( auto& e : res.atomDB ) {
      nc_assert_always(e.first.isElement());
      auto opt_stdnatelem = NC::AtomDB::getNaturalElement( e.first.Z() );
      bool isStdNaturalElem = opt_stdnatelem!=nullptr && opt_stdnatelem->getUniqueID() == e.second->getUniqueID();
      auto eName = NC::elementZToName(e.first.Z());
      std::cout<<"       | "<< eName<<" : "<<*e.second<<(isStdNaturalElem?" <<STANDARD-ELEM-DATA>>":"")<<std::endl;
    }
  };
  test("CO2");
  test("1.0xCO2");
  test("CO2/C_is_H");
  test("1.0xCO2/C_is_C13");
  test("air");
  test("air/C_is_C13");
  test("air/0.0relhumidity");
  test("air/0.50relhumidity");
  test("air/1.0relhumidity");

  //compare air density(relhumid and compare pts with: https://www.engineeringtoolbox.com/density-air-d_680.html
  test("air/0.0relhumidity/20C");//1.2041520963509085kgm3   rel 1.0
  test("air/0.8relhumidity/20C");//1.195743954182023kgm3    rel 0.9930
  test("air/0.0relhumidity/80C");//0.99956728598405442kgm3  rel 0.8301
  test("air/0.8relhumidity/80C");//0.85809965173152525kgm3  rel 0.7126
  //Checks out!!

  test("0.5xHe+0.5xAr");
  test("0.5xHe+0.5xAr/massfractions");
  test("0.5xAr+0.5xHe");
  test("0.5xAr+0.5xHe/massfractions");

  test("0.4xAr+0.5xHe+0.1xAr");
  test("0.4xHe+0.5xAr+0.1xHe");

  test("0.2xCO2+0.8xAr");
  test("0.2xCO2+0.8xAr/0.5relhumidity");

}

int main()
{
  test_parse();
  test_satvappressure_water();
  test_results();
  return 0;
}
