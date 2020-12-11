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

//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
// This file illustrates with a small example how one can add a new physics     //
// model to NCrystal. It is mainly intended for experts, and it should be noted //
// that the APIs used might be less stable than those exposed to people who are //
// just using NCrystal rather than extending it.                                //
//                                                                              //
// For simplicity everything is done in a single file here. In a more realistic //
// scenario, one should of course split the various classes into separate       //
// .hh/.cc files as usual in C++ and develop the model as an actual plugin for  //
// NCrystal. See: https://github.com/mctools/ncrystal/wiki/PluginsDevelopment   //
//                                                                              //
// As a silly example, the new physics model "SimpleIncElas", replaces the      //
// incoherent elastic physics from NCrystal with a model which is completely    //
// isotropic and has a constant cross section (incidently, this is actually how //
// some people approximate incoherent elastic scattering!). It has a single     //
// free parameter, namely the cross section in barn.                            //
//                                                                              //
// All other physics (Bragg diffraction and inelastic scattering) should be     //
// untouched.                                                                   //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////


//Include entire public NCrystal C++ API:
#include "NCrystal/NCrystal.hh"

//This example uses the randIsotropicScatterAngle function which is only
//available in an internal NCrystal header file (this is mainly intended as an
//example of how such an internal file can be used, not because it was really
//needed). Note that these internal files provide a lot of useful utilities for
//math, numerical integration, vectors, matrices, random sampling of
//distributions, etc. But it is important to note that the NCrystal developers
//do not really guarantee any form of long-term API stability for these internal
//files (we simply do not have the manpower to do so):
#include "NCrystal/internal/NCRandUtils.hh"

#include <iostream>

namespace NC = NCrystal;

//////////////////////////////////////////////////////////////////////////////////
//First we implement the "great" physics model:

class SimpleIncElasScatter : public NC::ScatterIsotropic {
  public:

  SimpleIncElasScatter( double sigma )
    : ScatterIsotropic("SimpleIncElasScatter"),
      m_sigma(sigma)
  {
  }
  double crossSectionNonOriented( double /*ekin*/ ) const final
  {
    return m_sigma;
  }
  void generateScatteringNonOriented( double /*ekin*/, double& angle, double& delta_ekin ) const final
  {
    delta_ekin = 0;
    angle = randIsotropicScatterAngle(getRNG());
  }
protected:
  virtual ~SimpleIncElasScatter() = default;
private:
  double m_sigma;
};

//////////////////////////////////////////////////////////////////////////////////
//Then we implement a factory which can provide this model (along with other physics):

class SimpleIncElasScatterFactory : public NC::FactoryBase {

  const char * getName() const final
  {
    return "SimpleIncElasScatterFactory";
  }

  int canCreateScatter( const NC::MatCfg& cfg ) const final
  {
    //Must return value >0 if we should do something, and a value higher than
    //100 means that we take precedence over the standard NCrystal factory:
    if (!cfg.get_incoh_elas())
      return 0;//incoherent-elastic disabled, do nothing
    auto info = globalCreateInfo(cfg);
    unsigned n_simpleincelas = info->countCustomSections("SIMPLEINCELAS");
    if (n_simpleincelas==0)
      return 0;//input file did not have @CUSTOM_SIMPLEINCELAS section, do nothing
    if (n_simpleincelas>1) {
      //File had multiple @CUSTOM_SIMPLEINCELAS sections, which we (and therefore no-one) supports.
      NCRYSTAL_THROW(BadInput,"Multiple @CUSTOM_SIMPLEINCELAS sections not supported")
    }
    //Ok, all good. Tell the framework that we want to deal with this, with a
    //higher priority than the standard factory gives (which is 100):
    return 999;
  }

  NC::RCHolder<const NC::Scatter> createScatter( const NC::MatCfg& cfg ) const final
  {
    nc_assert(canCreateScatter(cfg));//should only be called when possible

    //Let us dig out the cross section parameter:
    double sigma = parseFileDataForSigma(cfg);

    //Ok, time to instantiate our new model:
    NC::RCHolder<const NC::Scatter> sc_simpleincelas(new SimpleIncElasScatter(sigma));
    //Now we just need to combine this with all the other physics
    //(i.e. Bragg+inelastic).  So ask the framework to set this up, except for
    //incoherent-elastic physics of course since we are now dealing with that
    //ourselves:
    auto cfg2 = cfg.clone();
    cfg2.set_incoh_elas(false);
    auto sc_std = globalCreateScatter(cfg2);

    //Combine:
    return combineScatterObjects(sc_std,sc_simpleincelas);
  }

private:
  double parseFileDataForSigma( const NC::MatCfg& cfg ) const
  {
    //Parse the @CUSTOM_SIMPLEINCELAS section to extract the cross section value
    //specified by the user. Throw BadInput exception in case of issues:
    auto info = globalCreateInfo(cfg);
    auto data = info->getCustomSection( "SIMPLEINCELAS" );
    if (data.size()!=1||data.at(0).size()!=1) {
      NCRYSTAL_THROW(BadInput,"Bad format of @CUSTOM_SIMPLEINCELAS section"
                     " (must consist of just one line with a single parameter");
    }
    double sigma = std::atof(data.at(0).at(0).c_str());
    if (!(sigma>0.0&&sigma<1e9))
      NCRYSTAL_THROW(BadInput,"Invalid sigma value specified in @CUSTOM_SIMPLEINCELAS section"
                     " (must be >0 barn and <1e9 barn");
    return sigma;
  }
};

void setupCustomDataLibrary()
{
  //Add NCMAT file in-memory (to keep this example self-contained). The file has
  //a custom section: "@CUSTOM_SIMPLEINCELAS", which includes data needed by the
  //new physics model. We make it a bit extreme and fix the incoherent-elastic
  //xsect to 100barn:
  NC::registerInMemoryStaticFileData( "Al_simpleincelas100barn.ncmat",
                                      "NCMAT v3                         \n"
                                      "@CUSTOM_SIMPLEINCELAS            \n"
                                      "  100.0                          \n"
                                      "@CELL                            \n"
                                      " lengths 4.04958 4.04958 4.04958 \n"
                                      " angles 90 90 90                 \n"
                                      "@SPACEGROUP                      \n"
                                      "  225                            \n"
                                      "@ATOMPOSITIONS                   \n"
                                      "  Al 0   1/2 1/2                 \n"
                                      "  Al 0     0   0                 \n"
                                      "  Al 1/2 1/2   0                 \n"
                                      "  Al 1/2   0 1/2                 \n"
                                      "@DEBYETEMPERATURE                \n"
                                      "  Al   410.4                     \n");
}

int main() {

  NC::libClashDetect();//Detect broken installation

  //Register our custom factory with NCrystal, thus ensuring it gets invoked
  //when users call createScatter(..).
  NC::registerFactory(std::make_unique<SimpleIncElasScatterFactory>());

  //Add the Al_sg225_simpleincelas100barn.ncmat (virtual) file which has a
  //@CUSTOM_SIMPLEINCELAS section:
  setupCustomDataLibrary();

  //Ok, time to try it out. First let's see what happens when someone invokes
  //createScatter(..). First, with a standard aluminium file (which has no
  //@CUSTOM_SIMPLEINCELAS section):

  NC::RCHolder<const NC::Scatter> scatter_std_Al(NC::createScatter("Al_sg225.ncmat;dcutoff=0.5"));
  std::cout<<"Scatter created from Al_sg225.ncmat has crossSection(5.0Aa) = "
           <<scatter_std_Al->crossSectionNonOriented(NC::wl2ekin(5.0))<<"barn"<<std::endl;

  //That should have printed out something around 0.1-0.15 barn.

  //Now let's see what happens with our file which has a @CUSTOM_SIMPLEINCELAS section:
  NC::RCHolder<const NC::Scatter> scatter_custom_Al(NC::createScatter("Al_simpleincelas100barn.ncmat;dcutoff=0.5"));
  std::cout<<"Scatter created from Al_simpleincelas100barn.ncmat.ncmat has crossSection(5.0Aa) = "
           <<scatter_custom_Al->crossSectionNonOriented(NC::wl2ekin(5.0))<<"barn"<<std::endl;

  //That should have printed out something around 100.1-100.15 barn, clearly the
  //effect of our silly 100barn model :-)

  return 0;
}
