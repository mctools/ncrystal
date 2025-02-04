
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

#include "NCTestPlugin.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCMath.hh"

void NCP::customPluginTest()
{
  //This function is called by NCrystal after the plugin is loaded, but only if
  //the NCRYSTAL_PLUGIN_RUNTESTS environment variable is set to "1". In case of
  //errors or anything about the plugin not working, simply throw an exception
  //(which is what the nc_assert_always function does below, but feel free to
  //simply throw an exception directly).

  //Note, emit all messages here and elsewhere in plugin code with NCPLUGIN_MSG
  //(or NCPLUGIN_WARN for warnings), never raw usage of std::cout or printf!

  NCPLUGIN_MSG("Testing plugin");

  //Create some test NCMAT data. For simplicity we will here base it on an
  //existing file, but add our @CUSTOMPLUGIN section (see NCPhysicsModel.cc for
  //a description of the format).:

  constexpr double custom_incohelas_sigma = 5.0;//barn
  constexpr double custom_incohelas_wl_threshold = 1.2;//angstrom
  const std::string base_data = "stdlib::Al_sg225.ncmat";

  std::string testdata;
  {
    testdata = NC::FactImpl::createTextData(base_data)->rawDataCopy();
    std::ostringstream customsection;
    customsection << "@CUSTOM_" << pluginNameUpperCase() << "\n"
                  << "  # sigma and wavelength threshold\n  "
                  << NC::fmt(custom_incohelas_sigma) << " "
                  << NC::fmt(custom_incohelas_wl_threshold) << "\n";
    testdata += customsection.str();
  }

  if ( true ) {
    NCPLUGIN_MSG("Test NCMAT data begin:");
    NCRYSTAL_RAWOUT(testdata);
    NCPLUGIN_MSG("Test NCMAT data end.");
  }

  //Let us load and exercise this testdata. We do NOT want to register the
  //testdata as a virtual file, since we want this test to not leave anything
  //registered after it is done running. So we simply put the raw data directly
  //into a MatCfg object (which is the C++ object representing a cfg-string):

  auto cfg = NC::MatCfg::createFromRawData( std::string(testdata) );
  auto cfg_incelasonly = NC::MatCfg::createFromRawData( std::string(testdata),
                                                        ";comp=incoh_elas" );
  auto scat_incelasonly = NC::createScatter( cfg_incelasonly );
  auto scat = NC::createScatter( cfg );
  auto scat_base = NC::createScatter( base_data );
  auto scat_base_noincelas = NC::createScatter( base_data + ";incoh_elas=0" );

  for ( auto& wlval : NC::linspace(0.1,3.1,16) ) {
    auto wl = NC::NeutronWavelength( wlval );
    auto xs = scat.crossSectionIsotropic( wl );
    auto xs_incelasonly = scat_incelasonly.crossSectionIsotropic( wl );

    auto xs_base = scat_base.crossSectionIsotropic( wl );
    auto xs_base_noincelas = scat_base_noincelas.crossSectionIsotropic( wl );
    NCPLUGIN_MSG( "xs @ "<<wl<<" : "
                  <<xs_base<<" (base)"<<" "
                  <<xs<<" (new)"<<" "
                  <<xs_incelasonly<<" (new inc elas)"
                  );

    //Test that we get the expected values:
    double expected_xs_incelasonly = ( wlval < custom_incohelas_wl_threshold
                                       ? custom_incohelas_sigma : 0.0 );
    double expected_xs = expected_xs_incelasonly + xs_base_noincelas.dbl();
    nc_assert_always( NC::floateq( expected_xs_incelasonly,
                                   xs_incelasonly.dbl() ) );
    nc_assert_always( NC::floateq( expected_xs, xs.dbl() ) );
  }

  //Also test the data file provided in the data/ subdirectory of our plugin:
  auto scat_somefile
    = NC::createScatter( "plugins::DummyPlugin/somefile.ncmat;comp=incoh_elas");

  for ( auto& wlval : NC::linspace(0.1,3.1,16) ) {
    auto wl = NC::NeutronWavelength( wlval );
    auto xs = scat_somefile.crossSectionIsotropic( wl );
    NCPLUGIN_MSG( "somefile.ncmat incoh-elas xs @ "<<wl<<" : "<<xs);
    double expected_xs = ( wl.dbl() <= 4.0 ? 6.0 : 0.0 );
    nc_assert_always( NC::floateq( expected_xs, xs.dbl() ) );
  }

  //Note: We do not test the scattering here for this simple dummy plugin.

  NCPLUGIN_MSG("All tests of plugin were successful!");
}
