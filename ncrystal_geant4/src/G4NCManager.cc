////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2018 NCrystal developers                                   //
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

#include "G4NCrystal/G4NCManager.hh"
#include "NCrystal/NCScatter.hh"
#include "NCrystal/NCRandom.hh"
#include "NCrystal/NCFactory.hh"
#include "NCrystal/NCFactoryRegistry.hh"
#include "NCrystal/NCDefs.hh"

#include "G4Material.hh"
#include "Randomize.hh"
#include "globals.hh"

#include <sstream>
#include <iomanip>

namespace G4NCrystal {
  class G4WrappedRandGen  : public NCrystal::RandomBase {
    //Wraps G4's random stream for NCrystal usage.
  public:
    G4WrappedRandGen() {}
    virtual double generate() { return G4UniformRand(); }
  protected:
    virtual ~G4WrappedRandGen(){}
  };
}

G4NCrystal::Manager::Manager()
  : m_key("NCScat")
{
  NCrystal::setDefaultRandomGenerator(new G4WrappedRandGen);
}

G4NCrystal::Manager::~Manager()
{
  std::vector<const NCrystal::Scatter*>::iterator it, itE(m_scatters.end());
  for (it=m_scatters.begin();it!=itE;++it)
    (*it)->unref();
  NCrystal::setDefaultRandomGenerator(0);
}

G4NCrystal::Manager * G4NCrystal::Manager::s_mgr = 0;

G4NCrystal::Manager * G4NCrystal::Manager::getInstance()
{
  if (!s_mgr)
    s_mgr = new Manager;
  return s_mgr;
}
G4NCrystal::Manager * G4NCrystal::Manager::getInstanceNoInit()
{
  return s_mgr;
}

void G4NCrystal::Manager::addScatterProperty(G4Material* mat,const NCrystal::Scatter*scat)
{
  G4MaterialPropertiesTable* matprop = mat->GetMaterialPropertiesTable();
  if (matprop) {
    if (matprop->ConstPropertyExists(m_key.c_str())) {
      std::stringstream ss;
      ss<<"Setting property on material "<<mat->GetName()<<" more than once overrides previous content!";
      G4Exception ("G4NCrystal::Manager::addScatterProperty", "NCAddTwice",
                   JustWarning, ss.str().c_str());
    }
  } else {
    matprop = new G4MaterialPropertiesTable();
    mat->SetMaterialPropertiesTable(matprop);
  }

  unsigned idx(99999999);
  std::map<const NCrystal::Scatter*,unsigned>::const_iterator it = m_scat2idx.find(scat);
  if (it==m_scat2idx.end()) {
    idx = m_scatters.size();
    m_scat2idx[scat]=idx;
    m_scatters.push_back(scat);
    scat->ref();
  } else {
    //already known:
    idx = it->second;
  }
  nc_assert_always( unsigned(double(idx)) == idx );//make sure we can get the idx back out
  matprop->AddConstProperty(m_key.c_str(), idx);
}

void G4NCrystal::Manager::cleanup()
{
  if (s_mgr) {
    delete s_mgr;
    s_mgr=0;
  }
  NCrystal::clearInfoCaches();
  NCrystal::clearFactoryRegistry();
}

void G4NCrystal::Manager::handleError(const char* origin, unsigned id, NCrystal::Error::Exception& e)
{
  std::ostringstream s_code;
  s_code << "G4NC"<<e.getTypeName()<<std::setfill('0') << std::setw(3)<<id;
  G4Exception(origin, s_code.str().c_str(),FatalException, e.what());
}
