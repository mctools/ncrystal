#ifndef NCrystal_NCMatLoader_hh
#define NCrystal_NCMatLoader_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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

#include <map>
#include <string>
#include <vector>
#include <fstream>
#include "NCRotMatrix.hh"

namespace NCrystal {

  class Vector;
  class Info;

  class NCMatLoader {
  public:
    struct Lattice {
      RotMatrix cell;
      RotMatrix reciprocal;//2pi times cell.inv()
      std::map<std::string, std::vector<Vector> > atomic_pos;
    };

    NCMatLoader(const char* fn);
    ~NCMatLoader();
    std::map<std::string, unsigned> getAtomMap();
    void fillHKL(Info &info,  const std::map<std::string, double>& msdmap, double min_ds, double max_ds, bool expandhkl) const;
    void getWhkl(std::vector<double>& result, const double ksq, const std::vector<double> & msd) const;

    //get
    double getDebyeTemp() const;
    const std::map<std::string, double>& getDebyeMap() const;
    unsigned getAtomPerCell () const;
    unsigned getSpacegroupNum() const;
    const Lattice& getLattice() const;
    void getLatticeParameters(double &a, double &b, double &c, double &alpha, double &beta, double &gamma) const;
    std::string m_version;

  private:
    double m_debye_temp;
    std::map<std::string, double> m_debye_map;
    double m_a,m_b,m_c, m_alpha,m_beta,m_gamma;
    unsigned m_atomnum;
    unsigned m_sg;
    Lattice m_lattice;


  };
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCrystal {
  inline double NCMatLoader::getDebyeTemp() const {return m_debye_temp;}
  inline const std::map<std::string, double>& NCMatLoader::getDebyeMap() const {return m_debye_map;}
  inline unsigned NCMatLoader::getAtomPerCell () const  {return m_atomnum;}
  inline unsigned NCMatLoader::getSpacegroupNum() const { return m_sg;}
  inline const NCMatLoader::Lattice& NCMatLoader::getLattice() const { return m_lattice; }
  inline void NCMatLoader::getLatticeParameters(double &a, double &b, double &c, double &alpha, double &beta, double &gamma) const
  {a=m_a; b=m_b; c=m_c; alpha=m_alpha; beta=m_beta; gamma=m_gamma;  }
}

#endif
