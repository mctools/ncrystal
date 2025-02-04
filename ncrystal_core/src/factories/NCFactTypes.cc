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

#include "NCrystal/factories/NCFactTypes.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCFileUtils.hh"
namespace NC = NCrystal;

NC::TextDataPath::TextDataPath( const std::string& fp )
{
  auto i = fp.find("::"_s);
  if ( i < fp.size() ) {
    m_fact = fp.substr(0,i);
    m_path = fp.substr(i+2);
    trim(m_path);
    trim(m_fact);
  } else {
    m_path = fp;
    trim(m_path);
  }
  if ( m_path.empty() )
    NCRYSTAL_THROW(BadInput,"TextDataPath constructed with empty path!");

  if ( contains(m_path,'~') ) {
    if ( startswith(m_path,"~/") ) {
      const char * ev = std::getenv("HOME");
      std::string homedir(ev?ev:"");
      if ( homedir.empty() )
        NCRYSTAL_THROW(BadInput,"Data paths are only allowed to start with \"~/\" when the environment variable HOME is set");
      if ( contains(homedir,'~') )
        NCRYSTAL_THROW(BadInput,"The environment variable $HOME contains a \"~\" character!");
      m_path = path_join(homedir,m_path.substr(2));
    }
    if ( contains(m_path,'~') )
      NCRYSTAL_THROW(BadInput,"Data paths are not allowed to contain \"~\" characters except for an initial \"~/\".");
  }

  if ( !m_fact.empty() ) {
    for ( auto c : m_fact ) {
      if ( !isAlphaNumeric(c) && c!= '_' && c!='-' )
        NCRYSTAL_THROW2(BadInput,"TextDataPath invalid character in factory name: "<<c);
    }
    //Absolute paths and paths starting with "./" can only be served by
    //abspath/relpath factories. If anyone wants to e.g. take over all absolute
    //paths (e.g. a virtual filesystem), they will have to deregister the
    //standard abspath factory first and then register their own factory with
    //the same name.
    if ( startswith(m_path,"./") && m_fact!="relpath" )
      NCRYSTAL_THROW2(BadInput,"Paths starting with \"./\" are always relative and served by the "
                      "\"relpath\" factory. Thus they can not be requested with the factory: \""<<m_fact<<"\"");
    if ( path_is_absolute(m_path) && m_fact!="abspath" )
      NCRYSTAL_THROW2(BadInput,"Absolute paths like \""<<m_path<<"\" are always served"
                      " by the \"abspath\" factory. Thus they can not be requested with the factory: \""<<m_fact<<"\"");
  }
}
