#ifndef G4NCrystal_Install_hh
#define G4NCrystal_Install_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2023 NCrystal developers                                   //
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

namespace G4NCrystal {

  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // Call just after initialising the G4 run manager, in order to modify the  //
  // physics processes of neutrons and let any "NCrystal" properties          //
  // associated to the G4Materials take over the elastic hadronic physics     //
  // below 2eV (this is thus a run-time physics-list modification, and is an  //
  // alternative to the usual approach of hard-coded physics lists):          //
  //                                                                          //
  // NB: For this to work, your physics list must have installed exactly one  //
  //    active process derived from G4HadronElasticProcess for neutrons.      //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////

  void install();

  //Second version of install() which does nothing if no materials in the active
  //geometry has an "NCrystal" property (useful for framework implementers):

  void installOnDemand();

}

#endif
