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
#include <iostream>
namespace NC = NCrystal;

void performSomeWork( NC::Scatter&& sc)
{
  //Some dummy usage of the scatter model.
  std::cout<<"-------------------------------------------------------------"<<std::endl;
  for ( double wavelength = 0.1; wavelength<5.0; wavelength += 0.4 ) {
    NC::NeutronEnergy ekin = NC::NeutronWavelength{wavelength};
    auto xsect = sc.crossSectionIsotropic(ekin);
    auto outcome = sc.sampleScatterIsotropic( ekin );
    std::cout<<"Sigma(lambda="<<wavelength<<"Aa)="<<xsect
             <<". Scattering gave theta="
             <<std::acos(outcome.mu.dbl())*NC::kToDeg<<"degree and dE="<<(outcome.ekin.dbl()-ekin.dbl())
             <<"eV"<<std::endl;
  }
}

int main() {

  {

    //First try with the stdlib one:
    std::string cfgstr = "Al_sg225.ncmat;dcutoff=0.5;inelas=0;incoh_elas=0;temp=25C";
    performSomeWork(NC::createScatter( cfgstr ));

    //Now register similarly named file in memory, but with slightly larger lattice parameters:
    NC::registerInMemoryStaticFileData( "Al_sg225.ncmat",
                                        "NCMAT v2\n"
                                        "@CELL\n"
                                        " lengths 5 5 5\n"//orig: 4.04958 4.04958 4.04958
                                        " angles 90 90 90\n"
                                        "@SPACEGROUP\n"
                                        "  225\n"
                                        "@ATOMPOSITIONS\n"
                                        "  Al 0   1/2 1/2\n"
                                        "  Al 0     0   0\n"
                                        "  Al 1/2 1/2   0\n"
                                        "  Al 1/2   0 1/2\n"
                                        "@DEBYETEMPERATURE\n"
                                        "  Al   410.4\n");

    performSomeWork(NC::createScatter( cfgstr ));

    //Again, but with more usual lattice parameters:
    NC::registerInMemoryStaticFileData( "Al_sg225.ncmat",
                                        "NCMAT v2\n"
                                        "@CELL\n"
                                        " lengths 4.04958 4.04958 4.04958\n"
                                        " angles 90 90 90\n"
                                        "@SPACEGROUP\n"
                                        "  225\n"
                                        "@ATOMPOSITIONS\n"
                                        "  Al 0   1/2 1/2\n"
                                        "  Al 0     0   0\n"
                                        "  Al 1/2 1/2   0\n"
                                        "  Al 1/2   0 1/2\n"
                                        "@DEBYETEMPERATURE\n"
                                        "  Al   410.4\n");

    performSomeWork(NC::createScatter( cfgstr ));


    //Now with in-file cfg:
    NC::registerInMemoryStaticFileData( "Al_sg225.ncmat",
                                        "NCMAT v2\n"
                                        "# NCRYSTALMATCFG[dcutoffup=1.0]\n"
                                        "@CELL\n"
                                        " lengths 4.04958 4.04958 4.04958\n"
                                        " angles 90 90 90\n"
                                        "@SPACEGROUP\n"
                                        "  225\n"
                                        "@ATOMPOSITIONS\n"
                                        "  Al 0   1/2 1/2\n"
                                        "  Al 0     0   0\n"
                                        "  Al 1/2 1/2   0\n"
                                        "  Al 1/2   0 1/2\n"
                                        "@DEBYETEMPERATURE\n"
                                        "  Al   410.4\n");

    performSomeWork(NC::createScatter( cfgstr ));
  }

  //For completely pure valgrind:
  NC::clearCaches();
  fclose( stdin );
  fclose( stdout );
  fclose( stderr );

  return 0;
}
