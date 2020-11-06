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

#include "NCrystal/internal/NCAtomDBExtender.hh"
#include "NCrystal/internal/NCAtomDB.hh"
#include "NCrystal/internal/NCFactoryUtils.hh"
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCAtomUtils.hh"
#include "NCrystal/internal/NCMath.hh"
namespace NC = NCrystal;

void NC::AtomDBExtender::addData( const std::string& line, unsigned format_version )
{
  std::string l(line);
  trim(l);
  if (l.empty())
    NCRYSTAL_THROW(BadInput,"Invalid AtomDB specification (empty line)");
  if (!isSimpleASCII(line))
    NCRYSTAL_THROW2(BadInput,"Invalid AtomDB specification (must only contain simple ascii characters) :\""<<line<<"\"");
  VectS words;
  split(words,l);
  addData(words,format_version);
}

void NC::AtomDBExtender::addData( const NC::VectS& words, unsigned format_version )
{
  if (format_version==0)
    format_version = 9999;
  if (format_version<3)
    format_version=3;//earliest NCMAT version with any support for custom atomdb

  //Extensive validation (after this we know we are in correct format, and we
  //won't get any "nodefaults" lines):
  validateAtomDBLine( words, format_version );

  nc_assert(words.size()<std::numeric_limits<unsigned>::max());
  const unsigned nwords = static_cast<unsigned>(words.size());
  nc_assert(nwords>=3);

  auto lookupWithErrorMsg = [this](const std::string& n){
    auto data = lookupAtomData(n);
    if (!data)
      NCRYSTAL_THROW2(BadInput,"Invalid AtomDB specification (component \""
                      <<n<<"\" is not a known element, isotope, or mixture)");
    return data;
  };

  std::string label = words.at(0);

  //const double totfrac_tolerance = 1e-6;

  if (words.at(1)=="is") {
    if (nwords==3 || nwords == 4) {
      //simply alias, "X is Al" or "X is 1.0 Al".
      nc_assert( nwords==3 || str2dbl(words.at(2))==1.0 );
      populateDB(label,lookupWithErrorMsg(words.back()));
      return;
    }
    //Proper mixture:
    unsigned ncomponents = (nwords - 2) / 2;
    nc_assert(ncomponents>=2);
    AtomData::ComponentList components;
    components.resize(ncomponents);
    StableSum totalfraction;
    for (unsigned i = 0; i<ncomponents;++i) {
      std::size_t ifirstword = 2+2*i;
      auto& c = components.at(i);
      //important that next line is nc_assert_always, not just nc_assert!!!:
      nc_assert_always( safe_str2dbl(words.at(ifirstword),c.fraction) && !(c.fraction<=0) && !(c.fraction>1.0) );
      totalfraction.add(c.fraction);
      c.data = lookupWithErrorMsg(words.at(ifirstword+1));
    }
    //1e-9 here, 1e-10 in NCAtomUtils.cc (so safe to assert here):
    nc_assert(valueInInterval(1.0-1e-9,1.0+1e-9,totalfraction.sum()));
    double correction_factor = 1.0/totalfraction.sum();
    for (auto& c: components)
      c.fraction *= correction_factor;//"snap" to 1.0
    //Sort component list:
    std::stable_sort(components.begin(),components.end(),[](const AtomData::Component&a,const AtomData::Component&b)->bool
    {
      nc_assert(!!a.data  && !!b.data );
      if (a.fraction!=b.fraction)
        return a.fraction>b.fraction;
      return *a.data < *b.data;
    });
    populateDB(label,
               std::make_shared<const AtomData>(components));
    return;
  } else {
    nc_assert(endswith(words.at(1),"u"));
    auto getDblWithUnit = [](const std::string& s, const std::string& unit) -> double
    {
      nc_assert(endswith(s,unit));
      return str2dbl(s.substr(0,s.size()-unit.size()));
    };

    //Data entry:
    double mass = getDblWithUnit(words.at(1),"u"_s);
    double csl = getDblWithUnit(words.at(2),"fm"_s)*0.1;//fm=1e-15m -> sqrt(barn)=1e-14m
    double incxs = getDblWithUnit(words.at(3),"b"_s);
    double absxs = getDblWithUnit(words.at(4),"b"_s);
    nc_assert( !(mass<=0.0) );
    nc_assert( !(incxs<=0.0) );
    nc_assert( !(absxs<=0.0) );
    nc_assert( !ncisnan(mass) );
    AtomSymbol sbl(label);
    nc_assert(!sbl.isInvalid()&&!sbl.isCustomMarker());
    //Add entry for that element/isotope (note that label_trailingnumber is 0 if
    //we are not updating a specific isotope, which is as it should be for the
    //AtomData constructor):
    populateDB(label,
               std::make_shared<const AtomData>(SigmaBound{incxs},csl,absxs,mass,
                                                sbl.Z(),sbl.A()));
    return;
  }
  nc_assert(false);//should not get here
}

NC::AtomDataSP NC::AtomDBExtender::lookupAtomData(const std::string& lbl)
{
  auto it = m_db.find(lbl);
  if (it!=m_db.end())
    return it->second;
  if (m_allowInbuiltDB) {
    AtomDataSP ad = AtomDB::getIsotopeOrNatElem(lbl);
    if (!!ad)
      return ad;
  }
  AtomSymbol atomsymbol(lbl);
  NCRYSTAL_THROW2(BadInput,"Atom with label \""<<lbl<<"\" is unknown"
                  <<((atomsymbol.isIsotope()&&m_allowInbuiltDB)
                     ?". If it is a valid isotope which is simply missing in NCrystal's"
                     " internal database you must define it yourself":"")
                  <<(m_allowInbuiltDB?".":" (note that access to the inbuilt"
                     " database was disabled)."));
  return nullptr;
}

namespace NCrystal {
  namespace {
    //We group all entries in the global cache by their hash value. Whenever a
    //given AtomDBExtender instance populates the cache, it will check if it is
    //actually just replicating an AtomData instance which is already in the
    //cache. If so, it will populate it's own instance-specific cache with the
    //object from the global cache => same AtomData instance when e.g. loading
    //two different NCMAT files with the same new atom defined inside.
    static std::map<std::size_t,std::vector<AtomDataSP>> s_hash2atomdatas;
    static std::mutex s_hash2atomdatas_mutex;
  }
}

void NC::AtomDBExtender::clearGlobalCache()
{
  std::lock_guard<std::mutex> guard(s_hash2atomdatas_mutex);
  s_hash2atomdatas.clear();
}

void NC::AtomDBExtender::populateDB(const std::string& lbl, NC::AtomDataSP ad)
{
  nc_assert(!!ad);
  std::lock_guard<std::mutex> guard(s_hash2atomdatas_mutex);
  static bool first = true;
  if (first) {
    first = false;
    registerCacheCleanupFunction( clearGlobalCache );
  }

  std::size_t hashval = ad->hash();
  auto & v = s_hash2atomdatas[hashval];
  for (auto& existing_ad : v) {
    //The hashes match, so the additional check of sameValuesAs should ensure a
    //practically perfect comparison (we could of course also implement
    //operator==).
    if (ad->sameValuesAs(*existing_ad,1e-15,1e-15)) {
      //Found existing instance with exact same values, prefer that one:
      m_db[lbl] = existing_ad;
      return;
    }
  }
  v.push_back(ad);
  m_db[lbl] = std::move(ad);
}
