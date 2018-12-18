#ifndef NCrystal_ParseNCMAT_hh
#define NCrystal_ParseNCMAT_hh

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

#include "NCVector.hh"
#include <map>
#include <string>
#include <vector>
#include <fstream>

namespace NCrystal {

  class Vector;
  class Info;

  class NCMATParser {
  public:

    //Parse .ncmat files.

    //TODO for NC2: Needs update to be more robust, to support new data fields
    //and to optionally parse in-mem strings as well as files (to optionally
    //hardcode .ncmat contents in data files).

    NCMATParser(const char* fn);
    ~NCMATParser();

    std::map<std::string, unsigned> getAtomMap();//elementname -> number/cell
    double getDebyeTemp() const;//global debye temp
    const std::map<std::string, double>& getDebyeMap() const;//elementname -> per element debye temp
    const std::map<std::string, std::vector<Vector>  >& getAtomicPosMap() const;//elementname -> atomic coordinates
    unsigned getAtomPerCell () const;
    unsigned getSpacegroupNum() const;
    void getLatticeParameters(double &a, double &b, double &c, double &alpha, double &beta, double &gamma) const;

  private:
    std::string m_version;
    double m_debye_temp;
    std::map<std::string, double> m_debye_map;
    double m_a,m_b,m_c, m_alpha,m_beta,m_gamma;
    unsigned m_atomnum;
    unsigned m_sg;
    std::map<std::string, std::vector<Vector> > m_atomic_pos;
  };
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCrystal {
  inline double NCMATParser::getDebyeTemp() const {return m_debye_temp;}
  inline const std::map<std::string, double>& NCMATParser::getDebyeMap() const {return m_debye_map;}
  inline const std::map<std::string, std::vector<Vector>  >& NCMATParser::getAtomicPosMap() const { return m_atomic_pos; }
  inline unsigned NCMATParser::getAtomPerCell () const  {return m_atomnum;}
  inline unsigned NCMATParser::getSpacegroupNum() const { return m_sg;}
  inline void NCMATParser::getLatticeParameters(double &a, double &b, double &c, double &alpha, double &beta, double &gamma) const
  {a=m_a; b=m_b; c=m_c; alpha=m_alpha; beta=m_beta; gamma=m_gamma;  }
}

#endif
