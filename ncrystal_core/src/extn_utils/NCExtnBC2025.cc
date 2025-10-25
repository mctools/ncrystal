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

#include "NCrystal/internal/extn_utils/NCExtnBC2025.hh"

namespace NCBC2025 = NCrystal::Extn::BC2025;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///                         BC2025 Standard Recipes                          ///
///     Precision guarantee for x<1000: Error less than 1e-3*min(y,1-y)      ///
///        Reference: T. Kittelmann 2025 (publication in preparation)        ///
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double NCBC2025::y_primary( double x, double sintheta )
{
  double y0;
  if ( x < 0.1 ) {
    y0 = 1.0-x*(0.94285714-x*(0.8204-x*(0.593-x*(0.364-0.19*x))));
  } else {
    if ( x > 1e3 )
      return y_primary(1e3,sintheta)*std::sqrt(1e3/x);
    const double xs = std::sqrt(x);
    const double xp = (xs-1.0)/(xs+1.0);
    y0 = ( 0.518212-xp*(0.93036-xp*(0.182006+xp*(1.10097-xp*(0.62625
           +xp*(1.73562-xp*(1.08506+xp*(2.19459-xp*(1.40451+xp*(1.62083
           -xp*(1.1031+xp*(0.49125-0.357611*xp))))))))))) );
  }
  const double s = std::sqrt(sintheta);
  const double u = x*s;
  double ydelta;
  if ( u < 0.1 ) {
    ydelta = u*u*(0.41021645-u*(1.187-u*(2.37-4.18*u)));
  } else {
    const double us = std::sqrt(u);
    const double up = (us-1.0)/(us+1.0);
    ydelta = ( 0.05508+up*(0.1166-up*(0.2099+up*(0.5482-up*(0.5248+up*(1.402
               -up*(1.168+up*(2.096-up*(2.116+up*(1.155-up*(1.952
               -0.6046*up)))))))))) );
  }
  return y0 + sintheta*s*ydelta;
}

double NCBC2025::y_scndgauss( double x, double sintheta )
{
  double y0;
  if ( x < 0.1 ) {
    y0 = 1.0-x*(1.0606602-x*(0.9238-x*(0.667-x*(0.409-0.22*x))));
  } else {
    if ( x > 1e3 )
      return y_scndgauss(1e3,sintheta)*std::pow(x*1e-3,-0.933);
    const double xs = std::sqrt(x);
    const double xp = (xs-1.0)/(xs+1.0);
    y0 = ( 0.4588909-xp*(1.038687-xp*(0.2401003+xp*(1.288282-xp*(0.7641972
           +xp*(1.880246-xp*(1.886916+xp*(2.171852-xp*(3.273034
           +xp*(0.9771599-xp*(2.988445-xp*(0.4993548+xp*(1.037121
           -0.4353142*xp)))))))))))) );
  }
  const double s = std::sqrt(sintheta);
  const double u = x*s;
  double ydelta;
  if ( u < 0.1 ) {
    ydelta = u*u*(0.46188022-u*(1.333-u*(2.66-4.68*u)));
  } else {
    const double us = std::sqrt(u);
    const double up = (us-1.0)/(us+1.0);
    ydelta = ( 0.062289443+up*(0.13177896-up*(0.240705+up*(0.61857545
               -up*(0.61744404+up*(1.4812474-up*(1.5419561+up*(1.9976424
               -up*(2.8090858+up*(0.74297172-up*(2.3120683
               -0.8661981*up)))))))))) );
  }
  return y0 + sintheta*s*ydelta;
}

double NCBC2025::y_scndlorentz( double x, double sintheta )
{
  double y0;
  if ( x < 0.1 ) {
    y0 = 1.0-x*(1.0-x*(1.0667-x*(0.988-x*(0.79-0.55*x))));
  } else {
    if ( x > 1e3 )
      return y_scndlorentz(1e3,sintheta)*std::sqrt(1e3/x);
    const double xs = std::sqrt(x);
    const double xp = (xs-1.0)/(xs+1.0);
    y0 = ( 0.53379-xp*(0.84182-xp*(0.16806+xp*(0.65124-xp*(0.67623
           +xp*(0.47199-xp*(1.0872+xp*(0.030142-xp*(0.91361-xp*(0.28313
           +xp*(0.30078-0.1507*xp)))))))))) );
  }
  const double s = std::sqrt(sintheta);
  const double u = x*s;
  double ydelta;
  if ( u < 0.1 ) {
    ydelta = u*u*(0.53333333-u*(1.9753-u*(5.14-11.9*u)));
  } else {
    const double us = std::sqrt(u);
    const double up = (us-1.0)/(us+1.0);
    ydelta = ( 0.0514714+up*(0.0863117-up*(0.191581+up*(0.266342
               -up*(0.504516+up*(0.32195-up*(0.894662-up*(0.0162501
               +up*(0.708855-0.33707*up)))))))) );
  }
  return y0 + sintheta*s*ydelta;
}

double NCBC2025::y_scndfresnel( double x, double sintheta )
{
  double y0;
  if ( x < 0.1 ) {
    y0 = 1.0-x*(1.0-x*(0.88-x*(0.639-x*(0.394-0.21*x))));
  } else {
    if ( x > 1e3 )
      return y_scndfresnel(1e3,sintheta)*std::sqrt(1e3/x);
    const double xs = std::sqrt(x);
    const double xp = (xs-1.0)/(xs+1.0);
    y0 = ( 0.493354-xp*(0.963692-xp*(0.235067+xp*(1.18222-xp*(0.672931
           +xp*(1.78522-xp*(1.09976+xp*(2.10882-xp*(1.34721+xp*(1.46841
           -xp*(1.0054+xp*(0.426279-0.313436*xp))))))))))) );
  }
  const double s = std::sqrt(sintheta);
  const double u = x*s;
  double ydelta;
  if ( u < 0.1 ) {
    ydelta = u*u*(0.44-u*(1.278-u*(2.56-4.52*u)));
  } else {
    const double us = std::sqrt(u);
    const double up = (us-1.0)/(us+1.0);
    ydelta = ( 0.05839+up*(0.12063-up*(0.233343+up*(0.578753-up*(0.584531
               +up*(1.42753-up*(1.28278+up*(1.95436-up*(2.18561+up*(0.877761
               -up*(1.80505-0.599956*up)))))))))) );
  }
  return y0 + sintheta*s*ydelta;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///                          BC2025 Luxury Recipes                           ///
///     Precision guarantee for x<1000: Error less than 1e-6*min(y,1-y)      ///
///        Reference: T. Kittelmann 2025 (publication in preparation)        ///
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double NCBC2025::y_primary_lux( double x, double sintheta )
{
  double y0;
  if ( x < 0.02 ) {
    y0 = ( 1.0-x*(0.94285714285714-x*(0.8204329004329-x*(0.59332136724729
           -x*(0.36445371604802-0.19428532960037*x)))) );
  } else {
    if ( x > 1e3 )
      return y_primary_lux(1e3,sintheta)*std::sqrt(1e3/x);
    const double xs = std::sqrt(x);
    const double xp = (xs-1.0)/(xs+1.0);
    y0 = ( 0.5181739142759-xp*(0.9312354258708-xp*(0.1838441456426
           +xp*(1.1362742358-xp*(0.6439499617823+xp*(2.206331176315
           -xp*(1.140630599466+xp*(5.318659095481-xp*(1.818212320215
           +xp*(14.55391219715-xp*(10.87410353778+xp*(45.03938805729
           -xp*(110.6857924368+xp*(174.6473769273-xp*(720.5710673832
           +xp*(736.1717109119-xp*(3106.364394659+xp*(2617.91188764
           -xp*(9418.786974417+xp*(6988.557590101-xp*(20670.35083956
           +xp*(13504.21392907-xp*(33107.13743136+xp*(18520.07264798
           -xp*(38415.15440831+xp*(17532.51527802-xp*(31492.96870494
           +xp*(10877.45934053-xp*(17311.15064191+xp*(3975.170684325
           -xp*(5724.961892345+xp*(647.8950238727
           -860.4745834869*xp))))))))))))))))))))))))))))))) );
  }
  const double s = std::sqrt(sintheta);
  const double u = x*s;
  double ydelta;
  if ( u < 0.02 ) {
    ydelta = ( u*u*(0.41021645021645-u*(1.1866427344946-u*(2.3689491543122
               -4.1771345864079*u))) );
  } else {
    const double us = std::sqrt(u);
    const double up = (us-1.0)/(us+1.0);
    ydelta = ( 0.05508580666571+up*(0.1171543796275-up*(0.2094959893224
               +up*(0.5753610003144-up*(0.5029926535087+up*(1.827528591995
               -up*(0.9343233427577+up*(5.560154682465-up*(1.67846357379
               +up*(21.02643609634-up*(13.57802643233+up*(107.3363063748
               -up*(151.8892352219+up*(576.5632688968-up*(1071.234226007
               +up*(2533.754063934-up*(5066.801335663+up*(8223.815227808
               -up*(16955.1507836+up*(18939.97089727-up*(40926.64492568
               +up*(29716.04901001-up*(71039.43814205+up*(29047.97910648
               -up*(86612.37915478+up*(12774.94362795-up*(70438.23189688
               -up*(4826.668786858+up*(34308.55921383-up*(8466.342245127
               +up*(7568.031380324
               -2978.804894556*up)))))))))))))))))))))))))))))) );
  }
  return y0 + sintheta*s*ydelta;
}

double NCBC2025::y_scndgauss_lux( double x, double sintheta )
{
  double y0;
  if ( x < 0.02 ) {
    y0 = ( 1.0-x*(1.0606601717798-x*(0.9237604307034-x*(0.66666666666667
           -x*(0.40888100159996-0.21773242158073*x)))) );
  } else {
    if ( x > 1e3 )
      return y_scndgauss_lux(1e3,sintheta)*std::pow(x*1e-3,-0.933);
    const double xs = std::sqrt(x);
    const double xp = (xs-1.0)/(xs+1.0);
    y0 = ( 0.45898433395-xp*(1.0398885672-xp*(0.23347432787+xp*(1.3320564936
           -xp*(0.65964335645+xp*(2.4573554347-xp*(1.1858452452
           +xp*(6.0101020414-xp*(0.42409952962+xp*(16.109953924
           +xp*(7.6979330841-xp*(40.844896412+xp*(38.282141999
           -xp*(90.529504717+xp*(109.93053034-xp*(168.29068816
           +xp*(216.74542662-xp*(254.04407166+xp*(301.23376255
           -xp*(299.83692304+xp*(290.35481016-xp*(262.9383083
           +xp*(184.61968604-xp*(159.14627416+xp*(69.52062555
           -xp*(58.698294202+xp*(11.715653334
           -9.8564365005*xp)))))))))))))))))))))))))) );
  }
  const double s = std::sqrt(sintheta);
  const double u = x*s;
  double ydelta;
  if ( u < 0.02 ) {
    ydelta = ( u*u*(0.4618802153517-u*(1.3333333333333-u*(2.6577265103997
               -4.6812470639856*u))) );
  } else {
    const double us = std::sqrt(u);
    const double up = (us-1.0)/(us+1.0);
    ydelta = ( 0.06224953019339+up*(0.1328158690775-up*(0.236741722263
               +up*(0.6589627718965-up*(0.5420212226156+up*(2.036463256552
               -up*(0.97288478567+up*(5.846171066704-up*(0.4786858467243
               +up*(16.52467895771+up*(5.177873225504-up*(43.26855122074
               +up*(23.12266235959-up*(97.49162419575+up*(52.38175330609
               -up*(176.313696784+up*(68.31348886279-up*(238.6100261469
               +up*(41.13149813131-up*(222.4411069882-up*(10.56142327328
               +up*(125.8054728166-up*(30.66189892098+up*(32.28090591464
               -13.18049019718*up))))))))))))))))))))))) );
  }
  return y0 + sintheta*s*ydelta;
}

double NCBC2025::y_scndlorentz_lux( double x, double sintheta )
{
  double y0;
  if ( x < 0.02 ) {
    y0 = ( 1.0-x*(1.0-x*(1.0666666666667-x*(0.98765432098765
           -x*(0.79012345679012-0.55308641975309*x)))) );
  } else {
    if ( x > 1e3 )
      return y_scndlorentz_lux(1e3,sintheta)*std::sqrt(1e3/x);
    const double xs = std::sqrt(x);
    const double xp = (xs-1.0)/(xs+1.0);
    y0 = ( 0.5336441259-xp*(0.8429114501-xp*(0.1770013802+xp*(0.6797417272
           -xp*(0.818932607+xp*(0.7066833263-xp*(2.153895064
           +xp*(0.942873933-xp*(5.438402926+xp*(2.034470992
           -xp*(12.91149187+xp*(5.399270166-xp*(26.43349543
           +xp*(11.99710274-xp*(43.38620265+xp*(19.31504654
           -xp*(54.09291652+xp*(21.30361545-xp*(48.71558791
           +xp*(15.19814842-xp*(29.65815201+xp*(6.299367185
           -xp*(10.86478452+xp*(1.14907036
           -1.800708025*xp))))))))))))))))))))))) );
  }
  const double s = std::sqrt(sintheta);
  const double u = x*s;
  double ydelta;
  if ( u < 0.02 ) {
    ydelta = ( u*u*(0.53333333333333-u*(1.9753086419753-u*(5.1358024691358
               -11.891358024691*u))) );
  } else {
    const double us = std::sqrt(u);
    const double up = (us-1.0)/(us+1.0);
    ydelta = ( 0.05157982383+up*(0.08702046572-up*(0.199353768
               +up*(0.2899285322-up*(0.6378148804+up*(0.5502529572
               -up*(1.956981119+up*(0.9631382051-up*(5.54608499
               +up*(2.04851528-up*(14.45528159+up*(4.719780698
               -up*(32.73916116+up*(8.241956394-up*(59.63989821
               +up*(7.43970856-up*(81.0118243-up*(2.53170678+up*(75.52806893
               -up*(14.4309282+up*(42.64462669-up*(14.48078074
               +up*(10.92280777-5.137685796*up)))))))))))))))))))))) );
  }
  return y0 + sintheta*s*ydelta;
}

double NCBC2025::y_scndfresnel_lux( double x, double sintheta )
{
  double y0;
  if ( x < 0.02 ) {
    y0 = ( 1.0-x*(1.0-x*(0.88-x*(0.63915343915344-x*(0.39352481733434
           -0.2100936347603*x)))) );
  } else {
    if ( x > 1e3 )
      return y_scndfresnel_lux(1e3,sintheta)*std::sqrt(1e3/x);
    const double xs = std::sqrt(x);
    const double xp = (xs-1.0)/(xs+1.0);
    y0 = ( 0.4932694837076-xp*(0.9650580607045-xp*(0.2391707300191
           +xp*(1.231775501511-xp*(0.7162749186493+xp*(2.397392257122
           -xp*(1.286159487187+xp*(6.190338565047-xp*(2.15799495722
           +xp*(21.93947896347-xp*(11.66646771202+xp*(103.9329908487
           -xp*(110.3774765187+xp*(526.4911824684-xp*(699.0096620504
           +xp*(2257.913301253-xp*(2982.947865414+xp*(7410.986462373
           -xp*(9068.355396232+xp*(18053.43867394-xp*(20188.12192804
           +xp*(32202.67098807-xp*(33154.22980087+xp*(41377.0976758
           -xp*(39805.34287439+xp*(37198.46755365-xp*(33994.95728469
           +xp*(22166.48507268-xp*(19545.17021751+xp*(7852.794063735
           -xp*(6769.65264411+xp*(1249.922350984
           -1064.555223624*xp))))))))))))))))))))))))))))))) );
  }
  const double s = std::sqrt(sintheta);
  const double u = x*s;
  double ydelta;
  if ( u < 0.02 ) {
    ydelta = ( u*u*(0.44-u*(1.2783068783069-u*(2.5579113126732
               -4.5170131473465*u))) );
  } else {
    const double us = std::sqrt(u);
    const double up = (us-1.0)/(us+1.0);
    ydelta = ( 0.05839675969497+up*(0.1215403942474-up*(0.2327031821943
               +up*(0.6169197748584-up*(0.5576924427279+up*(1.966471723613
               -up*(0.9739495091226+up*(5.856571787922-up*(0.3468732618105
               +up*(19.93824999639+up*(5.12029535262-up*(85.66077650643
               +up*(6.233382554383-up*(410.7966508641-up*(135.2115982125
               +up*(1733.358492423-up*(1061.000115469+up*(5627.101189681
               -up*(4404.594142741+up*(13334.67031813-up*(12107.00455643
               +up*(22313.1140298-up*(23062.69191267+up*(25123.73631196
               -up*(30242.49736349+up*(17135.34812358-up*(26106.24295471
               +up*(4932.436511715-up*(13365.12956974-up*(1204.66538604
               +up*(3074.520821506
               -921.6011306758*up)))))))))))))))))))))))))))))) );
  }
  return y0 + sintheta*s*ydelta;
}
