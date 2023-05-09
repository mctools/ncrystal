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


//Include all NCrystal public headers:

#include "NCrystal/NCrystal.hh"

//This example uses the randIsotropicScatterMu function which is only available
//in an internal NCrystal header file (this is mainly intended as an example of
//how such an internal file can be used, not because it was really needed). Note
//that these internal files provide a lot of useful utilities for math,
//numerical integration, vectors, matrices, random sampling of distributions,
//etc. But it is important to note that the NCrystal developers do not really
//guarantee any form of long-term API stability for these internal files (we
//simply do not have the manpower to do so):

#include "NCrystal/internal/NCRandUtils.hh"

//We also use NCrystal::str2dbl from:
#include "NCrystal/internal/NCString.hh"

#include <iostream>

namespace NC = NCrystal;

//////////////////////////////////////////////////////////////////////////////////
//First we implement the "great" physics model:

class SimpleIncElasScatter final : public NC::ProcImpl::ScatterIsotropicMat {
public:

  SimpleIncElasScatter( NC::CrossSect sigma ) : m_sigma(sigma) {}

  const char * name() const noexcept override { return "SimpleIncElasScatter"; }

  NC::CrossSect crossSectionIsotropic( NC::CachePtr&, NC::NeutronEnergy ) const override
  {
    return m_sigma;
  }

  NC::ScatterOutcomeIsotropic sampleScatterIsotropic( NC::CachePtr&,
                                                      NC::RNG& rng,
                                                      NC::NeutronEnergy ekin ) const override
  {
    auto ekin_final = ekin; // elastic
    auto mu = randIsotropicScatterMu( rng ); //isotropic
    return { ekin_final, mu };
  }

private:
  NC::CrossSect m_sigma;
};

//////////////////////////////////////////////////////////////////////////////////
//Then we implement a factory which can provide this model (along with other physics):

class SimpleIncElasScatterFactory final : public NC::FactImpl::ScatterFactory {
public:

  const char * name() const noexcept override { return "SimpleIncElasScatterFactory"; }

  //Analyse cfg and determine if, and if so with what priority, the factory
  //can service the request:
  Priority query( const NC::FactImpl::ScatterRequest& cfg ) const override
  {
    //Must return Priority{value} if we can do something, and a value higher
    //than 100 means that we take precedence over the standard NCrystal
    //factory. If the request is not relevant, we return Priority{Unable}:

    if (!cfg.get_incoh_elas()) {
      //This factory is not relevant when incoherent-elastic is disabled.
      return Priority::Unable;
    }

    unsigned n_simpleincelas = cfg.info().countCustomSections("SIMPLEINCELAS");
    if (n_simpleincelas==0) {
      //This factory is only relevant if input file has @CUSTOM_SIMPLEINCELAS
      //section.
      return Priority::Unable;
    }

    if (n_simpleincelas>1) {
      //File had multiple @CUSTOM_SIMPLEINCELAS sections, which we do not support.
      NCRYSTAL_THROW(BadInput,"Multiple @CUSTOM_SIMPLEINCELAS sections not supported");
    }

    //Ok, all good. Tell the framework that we want to deal with this, with a
    //higher priority than the standard factory gives (which is 100):
    return Priority{999};
  }

  ProcPtr produce( const NC::FactImpl::ScatterRequest& cfg ) const override
  {
    //Service the request and produce the requested object (the framework will
    //only call this method after a previous call to query(..) indicated that
    //the production is possible - i.e. did not return Priority::Unable):

    //Let us dig out the cross section parameter (see the parseDataForSigma
    //method below for how it is done):
    NC::CrossSect sigma = parseDataForSigma(cfg);

    //Ok, time to instantiate our new model. We always use the makeSO function
    //instead of "new SimpleIncElasScatter" to do this (makeSO is similar to
    //std::make_shared):
    auto sc_simpleincelas = NC::makeSO<SimpleIncElasScatter>(sigma);

    //Now we just need to combine this with all the other physics
    //(i.e. Bragg+inelastic).  So ask the framework to set this up, except for
    //incoherent-elastic physics of course since we are now dealing with that
    //ourselves in sc_simpleincelas:
    auto sc_std = globalCreateScatter(cfg.modified("incoh_elas=0"));

    //Combine and return:
    return combineProcs( sc_std, sc_simpleincelas );
  }

private:
  NC::CrossSect parseDataForSigma( const NC::FactImpl::ScatterRequest& cfg ) const
  {
    //Parse the @CUSTOM_SIMPLEINCELAS section to extract the cross section value
    //specified by the user. Throw BadInput exception in case of issues:
    auto data = cfg.info().getCustomSection( "SIMPLEINCELAS" );
    if ( data.size() != 1 || data.at(0).size() != 1 ) {
      NCRYSTAL_THROW(BadInput,"Bad format of @CUSTOM_SIMPLEINCELAS section"
                     " (must consist of just one line with a single parameter");
    }
    double sigma = NC::str2dbl(data.at(0).at(0));
    if (!(sigma>0.0&&sigma<1e9))
      NCRYSTAL_THROW(BadInput,"Invalid sigma value specified in @CUSTOM_SIMPLEINCELAS section"
                     " (must be >0 barn and <1e9 barn");
    return NC::CrossSect{ sigma };
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
  NC::FactImpl::registerFactory(std::make_unique<SimpleIncElasScatterFactory>());

  //Add the Al_sg225_simpleincelas100barn.ncmat (virtual) file which has a
  //@CUSTOM_SIMPLEINCELAS section:
  setupCustomDataLibrary();

  //Ok, time to try it out. First let's see what happens when someone invokes
  //createScatter(..). First, with a standard aluminium file (which has no
  //@CUSTOM_SIMPLEINCELAS section):

  auto scatter_std_Al = NC::createScatter("Al_sg225.ncmat;dcutoff=0.5");
  std::cout<<"Scatter created from Al_sg225.ncmat has crossSection(5.0Aa) = "
           <<scatter_std_Al.crossSectionIsotropic(NC::NeutronWavelength{5.0})<<std::endl;

  //That should have printed out something around 0.1-0.15 barn.

  //Now let's see what happens with our file which has a @CUSTOM_SIMPLEINCELAS section:
  auto scatter_custom_Al = NC::createScatter("Al_simpleincelas100barn.ncmat;dcutoff=0.5");
  std::cout<<"Scatter created from Al_simpleincelas100barn.ncmat has crossSection(5.0Aa) = "
           <<scatter_custom_Al.crossSectionIsotropic(NC::NeutronWavelength{5.0})<<std::endl;

  //That should have printed out something around 100.1-100.15 barn, clearly the
  //effect of our silly 100barn model :-)

  return 0;
}
