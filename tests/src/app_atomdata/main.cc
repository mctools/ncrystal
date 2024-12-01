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

#include "NCrystal/interfaces/NCAtomData.hh"
#include "NCrystal/internal/atomdb/NCAtomDB.hh"
#include "NCrystal/internal/utils/NCAtomUtils.hh"
#include <iostream>
namespace NC = NCrystal;

int main()
{
  const double fm = 0.1;
  auto H_direct = std::make_shared<NC::AtomData>( NC::SigmaBound{80.26}, -3.7390*fm, NC::SigmaAbsorption{0.3326}, NC::AtomMass{1.007825032241}, 1 );
  auto H1 = std::make_shared<NC::AtomData>( NC::SigmaBound{80.27}, -3.7406*fm, NC::SigmaAbsorption{0.3326}, NC::AtomMass{1.0078}, 1, 1 );
  auto H2 = std::make_shared<NC::AtomData>( NC::SigmaBound{2.05}, 6.671*fm, NC::SigmaAbsorption{0.000519}, NC::AtomMass{2.01410177811}, 1, 2 );
  auto H3 = std::make_shared<NC::AtomData>( NC::SigmaBound{0.14}, 4.792*fm, NC::SigmaAbsorption{0.0}, NC::AtomMass{3.01604928199}, 1, 3 );
  std::cout<< *H1<<std::endl;
  std::cout<< *H2<<std::endl;
  std::cout<< *H3<<std::endl;

  auto H_composed = std::make_shared<NC::AtomData>( NC::AtomData::ComponentList{ {0.99985,H1}, {0.00015,H2} } );
  std::cout<< *H_composed<<std::endl;
  std::cout<< *H_direct<<std::endl;

  auto Ti_direct = std::make_shared<NC::AtomData>( NC::SigmaBound{2.87}, -3.438*fm, NC::SigmaAbsorption{6.09}, NC::AtomMass{47.867}, 22 );
  auto Ti46 = std::make_shared<NC::AtomData>( NC::SigmaBound{0.0}, 4.93*fm, NC::SigmaAbsorption{0.59}, NC::AtomMass{45.9526316}, 22, 46 );
  auto Ti47 = std::make_shared<NC::AtomData>( NC::SigmaBound{1.5}, 3.63*fm, NC::SigmaAbsorption{1.7}, NC::AtomMass{46.9517631}, 22, 47 );
  auto Ti48 = std::make_shared<NC::AtomData>( NC::SigmaBound{0.}, -6.08*fm, NC::SigmaAbsorption{7.84}, NC::AtomMass{47.9479463}, 22, 48 );
  auto Ti49 = std::make_shared<NC::AtomData>( NC::SigmaBound{3.3}, 1.04*fm, NC::SigmaAbsorption{2.2}, NC::AtomMass{48.9478700}, 22, 49 );
  auto Ti50 = std::make_shared<NC::AtomData>( NC::SigmaBound{0.}, 6.18*fm, NC::SigmaAbsorption{0.179}, NC::AtomMass{49.9447912}, 22, 50 );
  auto Ti_composed = std::make_shared<NC::AtomData>( NC::AtomData::ComponentList{ {0.082,Ti46},
                                                                                  {0.074,Ti47},
                                                                                  {0.738,Ti48},
                                                                                  {0.054,Ti49},
                                                                                  {0.052,Ti50} } );
  std::cout<<*Ti46<<std::endl;
  std::cout<<*Ti47<<std::endl;
  std::cout<<*Ti48<<std::endl;
  std::cout<<*Ti49<<std::endl;
  std::cout<<*Ti50<<std::endl;
  std::cout<<*Ti_direct<<std::endl;
  std::cout<<*Ti_composed<<std::endl;

  std::cout << *NC::AtomDB::getNaturalElement( 1 ) << std::endl;
  std::cout << *NC::AtomDB::getNaturalElement( 2 ) << std::endl;
  std::cout << *NC::AtomDB::getNaturalElement( 3 ) << std::endl;
  std::cout << *NC::AtomDB::getNaturalElement( "H" ) << std::endl;
  std::cout << *NC::AtomDB::getNaturalElement( "He" ) << std::endl;
  std::cout << *NC::AtomDB::getNaturalElement( "Li" ) << std::endl;
  std::cout << *NC::AtomDB::getIsotope( 1, 2 ) << std::endl;
  std::cout << *NC::AtomDB::getIsotope( "H2" ) << std::endl;
  std::cout << *NC::AtomDB::getIsotope( "D" ) << std::endl;
  std::cout <<"--------------------------------------------------------------------"<<std::endl;
  for (auto za : NC::AtomDB::getAllEntries()) {
    if (za.second==0)
      std::cout << *NC::AtomDB::getNaturalElement( za.first ) << std::endl;
    else
      std::cout << *NC::AtomDB::getIsotope( za.first, za.second ) << std::endl;
  }

  //chem forms:
  auto testChemForm = [](const char* cf_str)
  {
    std::cout<<"Try decoding chem form \""<<cf_str<<"\""<<std::endl;
    auto cf = NC::tryDecodeSimpleChemicalFormula( cf_str );
    if (!cf.has_value()) {
      std::cout<<"  -> Decoding failed!!"<<std::endl;
    } else {
      for ( auto& n_smb : cf.value() ) {
        nc_assert_always(n_smb.second.isElement());
        std::cout<<"  -> Got "<<n_smb.first<<" * "<<NC::elementZToName(n_smb.second.Z())<<std::endl;
      }
    }
    return cf;
  };
  testChemForm("Al2O3");
  testChemForm("  \t Al2O3");
  testChemForm("Al2 O3");
  testChemForm("Al  2O 3   ");
  testChemForm("Al  2O2 3   ");
  testChemForm("Al");
  testChemForm("H2O");
  testChemForm("Bg4");
  testChemForm("Be4");
  testChemForm("2H");
  testChemForm("Na3Al1OBe2");
  testChemForm("N a3Al1OBe2");
  testChemForm("Na3Al1O'Be2");
  testChemForm("H0O");
  testChemForm("H02O");
  testChemForm("H102O");
  testChemForm("H2OH3");
  {
    auto cf = NC::tryDecodeSimpleChemicalFormula( "H2OH3" );
    nc_assert_always( cf.has_value() && cf.value().size() == 2 );
  }


  auto cf_h2o = testChemForm("H2O");
  auto cf_oh2 = testChemForm("OH2");
  nc_assert_always( cf_h2o == cf_oh2);

}
