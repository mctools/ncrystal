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

#include "NCrystal/NCrystal.hh"
#include <iostream>

namespace NC = NCrystal;

int main()
{
  const char * cfgstr = "stdlib::Al_sg225.ncmat;dcutoff=0.5;temp=25C";
  std::cout<<"Verifying cross sections of \""<<cfgstr<<"\""<<std::endl;

  auto pc = NC::createScatter( cfgstr ).clone().clone();

  std::vector<double> refxsvals = {
    1.39667, 1.39437, 1.38778, 1.37679, 1.36133, 1.342, 1.32158, 1.30743,
    1.30778, 1.32854, 1.37329, 1.37085, 1.38283, 1.37629, 1.31355, 1.37167,
    1.35602, 1.34481, 1.44825, 1.2065, 1.29534, 1.32806, 1.42857, 1.53696,
    1.50743, 1.11265, 1.1864, 1.26397, 1.34521, 1.00269, 1.06099, 1.12169,
    1.18474, 1.2501, 1.31772, 1.38759, 1.45968, 1.53396, 1.61043, 1.68906,
    1.76985, 1.18806, 1.2403, 1.29387, 1.34877, 1.40498, 1.46251, 0.143847,
    0.144769, 0.145723, 0.146725, 0.147773, 0.148863, 0.14999, 0.151153,
    0.152347, 0.153569, 0.154815, 0.156124, 0.157444
  };

  for (unsigned i=0;i<60;++i) {
    auto wl = NC::NeutronWavelength{i * 0.1};
    auto xsect = pc.crossSectionIsotropic( wl );
    //std::cout<<", "<<xsect.dbl();
    auto ref = NC::CrossSect{ refxsvals.at(i) };
    if ( std::fabs( ref.dbl() - xsect.dbl() ) > 0.01 )
      NCRYSTAL_THROW2(CalcError,"Unexpected cross section at "
                      <<wl<<": "<<xsect<<" (expected "<<ref<<")");
  }

  //At 3.5Aa we only have Bragg cones at 119.6 and 96.9 degrees, with scattering
  //components contributing the following fractions:
  const auto wl_sample = NC::NeutronWavelength{ 3.5 };
  constexpr double frac_bragg = 0.899618;;
  constexpr double frac_incelas = 0.00554598;
  constexpr double frac_inel = 0.0948364;
  constexpr double angle1_deg = 96.9;
  constexpr double angle2_deg = 119.6;

  //Angle1 has around 61% of the bragg diffraction, and angle2 the rest:
  constexpr double relfrac_angle1 = 0.6119430177909819;
  constexpr double relfrac_angle2 = 1.0 - relfrac_angle1;

  std::cout<<"Verifying sample outcomes of \""<<cfgstr<<"\" at "
           <<wl_sample<<std::endl;

  struct ScatStat {
    std::size_t ntot = 0;
    std::size_t n_inelas = 0;
    std::size_t n_elas_angle1 = 0;
    std::size_t n_elas_angle2 = 0;
  };
  ScatStat stat;
  const NC::NeutronEnergy ekin_sample( wl_sample );
  for (unsigned i=0;i<1000000;++i) {
    ++stat.ntot;
    auto outcome = pc.sampleScatterIsotropic( ekin_sample );
    if ( ekin_sample == outcome.ekin ) {
      //elastic
      double angle_deg = std::acos(outcome.mu.get())*NC::kToDeg;
      if ( std::fabs( angle_deg - angle1_deg ) < 0.05 )
        ++stat.n_elas_angle1;
      else if ( std::fabs( angle_deg - angle2_deg ) < 0.05 )
        ++stat.n_elas_angle2;
    } else {
      ++stat.n_inelas;
    }
  }

  auto validate = []( const char * binname,
                      std::size_t ntot,
                      std::size_t n,
                      double expected_fraction )
  {
    const double x = double(n) / ntot;
    const double xdev = std::sqrt( n ) / ntot;
    const double ndev = ( xdev ? ( x - expected_fraction ) / xdev : 1e99 );
    std::cout<<"Found "<<binname<<" fraction [%]: "
             <<NC::fmt(x*100,"%.5g")<<" +- "<<NC::fmt(xdev*100,"%.2g")
             <<" (expected "<<NC::fmt(expected_fraction*100,"%.5g")
             <<", deviation is "<<NC::fmt(ndev,"%.2g")<<" stddev)"<<std::endl;
    if ( xdev <= 0.0 || xdev > 0.005 )
      NCRYSTAL_THROW2(CalcError,"Too low statistics in bin "<<binname);
    if ( std::fabs( x - expected_fraction ) > 4.0 * xdev )
      NCRYSTAL_THROW2(CalcError,"Unexpected result in bin "<<binname);
  };

  std::size_t n_elas = stat.ntot - stat.n_inelas;
  std::size_t n_bragg = stat.n_elas_angle1 + stat.n_elas_angle2;
  std::size_t n_incelas = n_elas - n_bragg;

  validate( "inelas", stat.ntot, stat.n_inelas, frac_inel );
  validate( "incelas", stat.ntot, n_incelas, frac_incelas );
  //  validate( "bragg", stat.ntot, n_bragg, frac_bragg );
  validate( "bragg (angle 1)", stat.ntot, stat.n_elas_angle1,
            frac_bragg*relfrac_angle1 );
  validate( "bragg (angle 2)", stat.ntot, stat.n_elas_angle2,
            frac_bragg*relfrac_angle2 );


  //Now validate single crystal slightly:
  NC::MatCfg cfg("stdlib::Ge_sg227.ncmat;dcutoff=0.5");
  cfg.set_mos( NC::MosaicityFWHM{ 40.0*NC::kArcSec });
  cfg.set_dir1( NC::HKLPoint{5,1,1}, NC::LabAxis{0,0,1} );
  cfg.set_dir2( NC::HKLPoint{0,-1,1}, NC::LabAxis{0,1,0} );
  auto sc =  NC::createScatter( cfg );
  auto wlsc = NC::NeutronWavelength{1.540};
  {
    auto xsectsc = sc.crossSection( wlsc, { 0.0, 1.0, 1.0 } );
    std::cout << "singlecrystal Ge x-sect at "<<wlsc<<" Aa is "<<xsectsc
              <<" barn (orientation 1)"<<std::endl;
    auto expect = NC::CrossSect{ 591.026 };
    if ( std::fabs(xsectsc.dbl()-expect.dbl())>1.0 )
      NCRYSTAL_THROW2(CalcError,
                      "Unexpected cross section! (expected "<<expect<<")");
  }

  {
    auto xsectsc = sc.crossSection( wlsc, { 1.0, 1.0, 0.0 } );
    std::cout << "singlecrystal Ge x-sect at "<<wlsc<<" Aa is "<<xsectsc
              <<" barn (orientation 2)"<<std::endl;
    auto expect = NC::CrossSect{ 1.6676 };
    if ( std::fabs(xsectsc.dbl()-expect.dbl())>0.2 )
      NCRYSTAL_THROW2(CalcError,
                      "Unexpected cross section! (expected "<<expect<<")");

  }
  return 0;
}
