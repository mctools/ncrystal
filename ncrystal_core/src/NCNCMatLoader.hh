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
    struct Lattice{
      //TODO for NC2 for XX (general comment): do not prefix public data members with "m_" (only if private/protected).
      RotMatrix m_cell;
      RotMatrix m_reciprocal;//2pi times m_cell.inv()
      std::map<std::string, std::vector<Vector> > m_atomic_pos;
    };


    NCMatLoader(const char* fn);
    ~NCMatLoader();
    std::map<std::string, unsigned> getAtomMap();
    unsigned getNumAtom() const;
    unsigned getSpacegroupNum() const;
    const Lattice& getLattice() const;
    void fillHKL(Info &info,  const std::map<std::string, double>& msdmap, double min_ds, double max_ds, bool expandhkl) const;
    void getWhkl(std::vector<double>& result, const double ksq, const std::vector<double> & msd) const;


    //TODO for NC2: Most of these should be private and with accessor
    //methods. In any case, the "m_" prefix should not be used for public data
    //member (this should be fixed in a looot of files, not just here):

    std::string m_version;
    Lattice m_lattice;
    double m_a,m_b,m_c, m_alpha,m_beta,m_gamma;
    double m_debye_temp;
    std::map<std::string, double> m_debye_map;
    unsigned m_atomnum;
    unsigned m_sg;


  private:
    std::ifstream& ignoreCharNTimes(std::ifstream& file, unsigned num, const char& c='\n');

  };
}

#endif
