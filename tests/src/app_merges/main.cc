////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
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
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/pcbragg/NCPowderBragg.hh"
#include "NCrystal/internal/elincscatter/NCElIncScatter.hh"

#include <iostream>
#include <fstream>
namespace NC = NCrystal;

int main() {

  ///////////////////////////////// AbsOOV /////////////////////////////////////////////
  {
    auto abs1 = NC::createAbsorption("Al_sg225.ncmat");
    auto abs2 = NC::createAbsorption("Ni_sg225.ncmat");

    auto abs_merged = NC::Absorption( abs1.underlying().createMerged(abs2.underlying(),1.0,1.0) );

    auto xs1 = abs1.crossSectionIsotropic(NC::const_ekin_2200m_s);
    auto xs2 = abs2.crossSectionIsotropic(NC::const_ekin_2200m_s);
    auto xs_merged = abs_merged.crossSectionIsotropic(NC::const_ekin_2200m_s);
    std::cout<<"sigma_abs@2200m/s (1) = "<<xs1<<std::endl;
    std::cout<<"sigma_abs@2200m/s (2) = "<<xs2<<std::endl;
    std::cout<<"sigma_abs@2200m/s (merged) = "<<xs_merged<<std::endl;
    nc_assert_always( NC::ncabs( (xs1.dbl() + xs2.dbl()) - xs_merged.dbl() ) < 1e-15 );
  }

  ///////////////////////////////// PowderBragg /////////////////////////////////////////////
  {
    //Fake hkl lists (VectDFM is vector of (d-spacing,fsq*mult) pairs):
    const double v0_times_natoms_1  = 1.0;
    const double v0_times_natoms_2  = 1.2;
    NC::PowderBragg::VectDFM planes1 = { {0.2,500.0},  { 0.4, 20.0 }, { 2.0, 0.5 }, { 2.4+1e-12, 0.1 }, { 4.5, 0.3 } };
    NC::PowderBragg::VectDFM planes2 = { {0.05,30000.1}, { 2.4, 1.0 }, { 5.0, 0.1 } };

    auto pcb1 = NC::makeSO<NC::PowderBragg>( v0_times_natoms_1, std::move(planes1) );
    auto pcb2 = NC::makeSO<NC::PowderBragg>( v0_times_natoms_2, std::move(planes2) );
    auto pcb_merged = pcb1->createMerged(pcb2,1.0,1.0);

    std::string fn = "testpcbraggmerge.txt";
    std::cout<<"Writing "<<fn<<" (visualise with sb_ncperfval_plottxtfile script)"<<std::endl;
    std::ofstream fh(fn);
    fh << "# NB: Visualise this file with sb_ncperfval_plottxtfile\n";
    auto xs = [](const NC::ProcImpl::ProcPtr& sc,double thewl)
    {
      NC::CachePtr dummy;
      return sc->crossSectionIsotropic( dummy, NC::NeutronWavelength{thewl} ).dbl();
    };
    fh << "# colnames = PowderBragg1 ; PowderBragg2 ; PowderBragg1+PowderBragg2 ; PowderBragg_merged \n";
    fh << "# plotstyle = - \n";
    fh << "# alpha = 0.4 \n";
    fh << "# xlabel = angstrom \n";
    fh << "# ylabel = barn \n";
    for ( auto wl : NC::linspace(0.01,10.5,10000) ) {
      const double xs1 = xs(pcb1,wl);
      const double xs2 = xs(pcb2,wl);
      const double xs_merged = xs(pcb_merged,wl);
      nc_assert_always( xs1 >= 0.0 && xs2 >= 0.0 && xs_merged >= 0.0 );
      nc_assert( NC::ncabs( ( xs1+xs2 ) - xs_merged ) < 1e-10 );
      fh << wl << " "<< xs1 << " "<< xs2 << " "<<(xs1+xs2)<<" "<<xs_merged << std::endl;
    }
  }

  ///////////////////////////////// ElIncScatter /////////////////////////////////////////////
  {
    NC::VectD elem_msd1   = { 0.01, 0.005, 0.02 };
    NC::VectD elem_xs1    = { 1.9, 1.1, 4.1 };
    NC::VectD elem_scale1 = { 0.5, 1.5, 0.1 };

    NC::VectD elem_msd2   = { 0.005, 0.03123 };
    NC::VectD elem_xs2    = { 0.5,  1.7 };
    NC::VectD elem_scale2 = { 0.9,  1.23 };

    //Constructor similar to the ElIncXS constructor:
    auto elinc1 = NC::makeSO<NC::ElIncScatter>( elem_msd1, elem_xs1, elem_scale1 );

    //Constructor similar to the ElIncXS constructor:
    auto elinc2 = NC::makeSO<NC::ElIncScatter>( elem_msd2, elem_xs2, elem_scale2 );

    auto elinc_merged = elinc1->createMerged(elinc2,1.0,1.0);

    std::string fn = "testelincscattermerge.txt";
    std::cout<<"Writing "<<fn<<" (visualise with sb_ncperfval_plottxtfile script)"<<std::endl;
    std::ofstream fh(fn);
    fh << "# NB: Visualise this file with sb_ncperfval_plottxtfile\n";
    auto xs = [](const NC::ProcImpl::ProcPtr& sc,double thewl)
    {
      NC::CachePtr dummy;
      return sc->crossSectionIsotropic( dummy, NC::NeutronWavelength{thewl} ).dbl();
    };
    fh << "# colnames = ElInc1 ; ElInc2 ; ElInc1+ElInc2 ; ElInc_merged \n";
    fh << "# plotstyle = - \n";
    fh << "# alpha = 0.4 \n";
    fh << "# xlabel = angstrom \n";
    fh << "# ylabel = barn \n";
    for ( auto wl : NC::linspace(0.01,10.5,10000) ) {
      const double xs1 = xs(elinc1,wl);
      const double xs2 = xs(elinc2,wl);
      const double xs_merged = xs(elinc_merged,wl);
      nc_assert_always( xs1 >= 0.0 && xs2 >= 0.0 && xs_merged >= 0.0 );
      nc_assert( NC::ncabs( ( xs1+xs2 ) - xs_merged ) < 1e-10 );
      fh << wl << " "<< xs1 << " "<< xs2 << " "<<(xs1+xs2) << " "<<xs_merged<< std::endl;
    }
  }
  //TODO: unit test e.g. pcbragg empty cases, and scale factors.

  return 0;
}
