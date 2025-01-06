
//Define exported symbols used by NCrystal to activate the plugin. This is
//mostly generic, except for the registerPlugin function which can potentially
//be modified in some use-cases.

#define NCPLUGIN_BOILERPLATE_CC
#include "NCrystal/NCPluginBoilerplate.hh"

#include "NCPluginFactory.hh"

void NCP::registerPlugin()
{
  //This function is required for the plugin to work. It should register
  //factories (or potentially other stuff, e.g. adding in-mem data files, etc.)
  //for the plugin.
  NC::FactImpl::registerFactory(std::make_unique<NCP::PluginFactory>());
};
