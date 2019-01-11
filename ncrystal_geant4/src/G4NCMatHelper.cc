////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2019 NCrystal developers                                   //
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

#include "G4NCrystal/G4NCMatHelper.hh"
#include "G4NCrystal/G4NCManager.hh"
#include "NCrystal/NCDefs.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCScatter.hh"
#include "NCrystal/NCFactory.hh"
#include "G4NistManager.hh"
#include <sstream>
#include <stdint.h>//cstdint hdr only in C++11

namespace G4NCrystal {

  unsigned greatestCommonDivisor(unsigned a, unsigned b) {
    while (b) {
      unsigned tmp = a % b;
      a = b;
      b = tmp;
    }
    return a;
  }

  typedef std::vector<std::pair<unsigned,unsigned> > ChemicalFormula;//(atomic_number,count) pairs

  void getChemicalFormula(const NCrystal::Info* info, ChemicalFormula& cf, unsigned& ntot) {

    //TODO for NC2: In principle the chemical formula derivation is not specific
    //to Geant4, and could move into NCrystal itself, if useful for some
    //reason. Info could for instance be added in NCInfo::objectDone().

    //extract raw info:
    if (!info->hasAtomInfo())
      NCRYSTAL_THROW(MissingInfo,"Selected crystal info source lacks info about atomic composition.");

    NCrystal::AtomList::const_iterator it = info->atomInfoBegin();
    NCrystal::AtomList::const_iterator itE = info->atomInfoEnd();
    nc_assert(it!=itE);
    cf.clear();
    cf.reserve(itE-it);
    ntot = 0;
    for ( ; it!=itE; ++it ) {
      if (it->atomic_number>150)
        NCRYSTAL_THROW2(BadInput,"invalid atomic number ("<<it->atomic_number<<")");
      if (it->number_per_unit_cell) {
        ntot += it->number_per_unit_cell;
        cf.push_back(ChemicalFormula::value_type(it->atomic_number,it->number_per_unit_cell));
      }
    }
    if (cf.empty())
      NCRYSTAL_THROW(BadInput,"Atomic composition info indicates an empty unit cell.");

    //reduce chemical formula by getting rid of greatest common divisor in
    //counts (e.g. Al2O3 instead of Al6O9):
    unsigned thegcd = cf.at(0).second;
    for (size_t i = 1; i<cf.size(); ++i)
      thegcd = greatestCommonDivisor(thegcd,cf.at(i).second);
    for (size_t i = 0; i<cf.size(); ++i)
      cf.at(i).second /= thegcd;
    ntot /= thegcd;

    //sort by atomic number:
    std::sort(cf.begin(),cf.end());

    //sanity check that each element only appears once:
    for (size_t i = 0; i+1<cf.size(); ++i)
      if (cf.at(i).first == cf.at(i+1).first)
        NCRYSTAL_THROW(BadInput,"Atomic composition info has duplicate entries for same atomic number.");
  }


  G4String getChemicalFormulaInHillSystemString( const ChemicalFormula& chemform ) {

    //Convert chemical formula to a string with element symbols and counts,
    //using the Hill System for providing the element order:
    //
    //  https://en.wikipedia.org/wiki/Chemical_formula#Hill_system
    //
    //Note that we do not implement any exceptions to the Hill system (e.g. O2
    //last in oxides, positive ion first in ionic compounds, etc.).
    //
    //This is done by constructing the following string keys to sort by:
    //
    // 1) Default to the element symbol ("Al","H","B","Be",...)
    // 2) If there is carbon anywhere in the entire formula, the sort key for
    //    carbon becomes "1" (i.e. it comes first) and for hydrogen it becomes
    //    "2" (i.e. it comes second)

    const std::vector<G4String>& symbol_db = G4NistManager::Instance()->GetNistElementNames();

    //sort key:
    std::vector<std::pair<std::string,std::string> > hillsystem_sort;
    ChemicalFormula::const_iterator itCF, itCF_End(chemform.end());

    bool any_carbon(false);
    for ( itCF = chemform.begin(); itCF != itCF_End; ++itCF ) {
      if (itCF->first==6) {
        any_carbon = true;
        break;
      }
    }
    for ( itCF = chemform.begin(); itCF != itCF_End; ++itCF ) {
      std::ostringstream symbandcount;
      std::string sortkey;
      if (any_carbon && (itCF->first==1||itCF->first==6) ) {
        sortkey = ( itCF->first==6 ? "1" : "2" );
      }
      if ( itCF->first < symbol_db.size() ) {
        symbandcount << symbol_db.at(itCF->first);
      } else {
        //Fall-back, unlikely to ever happen. Name as "Elem<xxx>" and put at end
        //of formula:
        symbandcount << "Elem<"<<itCF->first<<">";
        sortkey = "{"+symbandcount.str();
      }
      if (sortkey.empty())
        sortkey = symbandcount.str();
      if (itCF->second!=1)
        symbandcount << itCF->second;
      hillsystem_sort.push_back(std::pair<std::string,std::string>(sortkey,symbandcount.str()));
    }

    std::sort(hillsystem_sort.begin(),hillsystem_sort.end());
    std::vector<std::pair<std::string,std::string> >::const_iterator itHS;
    std::ostringstream ss;
    for ( itHS = hillsystem_sort.begin(); itHS != hillsystem_sort.end(); ++itHS )
      ss << itHS->second;
    return ss.str();
  }

  G4Material * getBaseG4MaterialWithCache( const NCrystal::Info* info ) {

    //Construct base material for a given relative atomic composition.

    //Figure out the (reduced) chemical formula, based on the atoms in the unit cell:
    ChemicalFormula chemform;
    unsigned chemform_ntot;
    getChemicalFormula(info, chemform, chemform_ntot);
    nc_assert( !chemform.empty() && chemform_ntot );

    //And the string representation:
    G4String chemform_str = getChemicalFormulaInHillSystemString(chemform);

    //Using chemform_str as key, check if we already have this material available:
    static std::map<std::string,size_t> cache;

    std::map<std::string,size_t>::iterator itCache = cache.find(chemform_str);
    if ( itCache != cache.end() ) {
      G4Material * cachedmat = G4Material::GetMaterialTable()->at(itCache->second);
      if (cachedmat)//might be NULL if material was deleted
        return cachedmat;
    }

    //Finally, create the base material based purely on the reduced chemical
    //formula (1atm, 293.15K, solid, 1g/cm3 - note 293.15 is the NCrystal default!):
    G4String matnameprefix("NCrystalBaseMat::");
    G4Material * mat = new G4Material( matnameprefix+chemform_str,
                                       1.0*CLHEP::gram/CLHEP::cm3,
                                       chemform.size(), kStateSolid,
                                       293.15 * CLHEP::kelvin,
                                       1.0 * CLHEP::atmosphere);

    ChemicalFormula::const_iterator itCF, itCF_End(chemform.end());
    for ( itCF = chemform.begin(); itCF != itCF_End; ++itCF ) {
      nc_assert( (uint64_t)(itCF->second) <= (uint64_t)std::numeric_limits<G4int>::max() );
      mat->AddElement( G4NistManager::Instance()->FindOrBuildElement(itCF->first,true),
                       G4int(itCF->second) );
    }
    mat->SetChemicalFormula(chemform_str);

    //Add to cache and return:
    cache[chemform_str] = mat->GetIndex();
    return mat;
  }

  G4Material * getG4MaterialWithCache( const NCrystal::MatCfg& cfg ) {

    //check if we already have this material available:


    std::string cache_key = cfg.toStrCfg(true,0);//NB: uses file name as
                                                 //specified, which is the
                                                 //correct thing to do (to avoid
                                                 //absolute paths in material
                                                 //names).
    static std::map<std::string,size_t> cache;
    std::map<std::string,size_t>::iterator itCache = cache.find(cache_key);
    if ( itCache != cache.end() ) {
      G4Material * cachedmat = G4Material::GetMaterialTable()->at(itCache->second);
      if (cachedmat)//might be NULL if material was deleted
        return cachedmat;
    }

    //First see if we are at all able to initialise an NCrystal scatter object
    //for the given configuration:
    NCrystal::RCHolder<const NCrystal::Scatter> scatter(NCrystal::createScatter(cfg));

    //Get NCrystal info object + G4 base material for the given chemical composition:
    NCrystal::RCHolder<const NCrystal::Info> info(NCrystal::createInfo(cfg));

    if (!info.obj()->hasDensity())
      NCRYSTAL_THROW(MissingInfo,"Selected crystal info source lacks info about material density.");

    G4Material * matBase = getBaseG4MaterialWithCache( info.obj() );

    //Create derived material with specific density, temperature and NCrystal scatter physics:

    //NB: Default temperature is same as NC::MatCfg's default (293.15) rather
    //than G4's STP (273.15). It will anyway always be overridden in the derived
    //material, but we do it like this to avoid two different temperatures even
    //when the user didn't specify anything.

    double temp = info.obj()->hasTemperature()
      ? info.obj()->getTemperature()*CLHEP::kelvin
      : 293.15*CLHEP::kelvin;

    G4String matnameprefix("NCrystal::");
    G4Material * mat = new G4Material( matnameprefix+cache_key,
                                       cfg.get_packfact()*info.obj()->getDensity()*(CLHEP::gram/CLHEP::cm3),
                                       matBase, kStateSolid, temp, 1.0 * CLHEP::atmosphere );

    G4NCrystal::Manager::getInstance()->addScatterProperty(mat,scatter.obj());

    //Add to cache and return:
    cache[cache_key] = mat->GetIndex();
    return mat;
  }
}

G4Material * G4NCrystal::createMaterial( const char * cfgstr )
{
  try {
    NCrystal::MatCfg cfg(cfgstr);
    return getG4MaterialWithCache(cfg);
  } catch ( NCrystal::Error::Exception& e ) {
    Manager::handleError("G4NCrystal::createMaterial",101,e);
  }
  return 0;
}

G4Material * G4NCrystal::createMaterial( const G4String& cfgstr )
{
  return createMaterial(cfgstr.c_str());
}

G4Material * G4NCrystal::createMaterial( const NCrystal::MatCfg&  cfg )
{
  try {
    return getG4MaterialWithCache(cfg);
  } catch ( NCrystal::Error::Exception& e ) {
    Manager::handleError("G4NCrystal::createMaterial",101,e);
  }
  return 0;
}

