#ifndef NCrystal_NeutronSCL_hh
#define NCrystal_NeutronSCL_hh

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

#include "NCrystal/NCDefs.hh"

namespace NCrystal {
  // a data table class for neutron scattering length (SCL) and its power (i.e. cross section)
  // scattering data are the recommended values in [1] according to [2,3]
  //[1] " Table of coherent scattering lengths and cross sections",
  //available at http://www.ati.ac.at/~neutropt/scattering/Scattering_lengths_table_20010419.pdf
  //[2] LANDOLT-BORNSTEIN, New Series I/16A (Ed. H. Schopper) Chap.6, Springer, Berlin 2000
  //[3] Atomic Data Nuclear Data Tables 49 (1991) 65

  class NeutronSCL {
  public:
    struct IsotopeData  {
      unsigned atomic_num; //Z
      double mass;    //in atomic mass units (Dalton)
      double coh_sl;  //in (barn)^0.5
      double inc_xs;  //in barn
      double cap_xs;  //in barn
      IsotopeData(unsigned z, double m, double csl, double ixs, double cxs)
      :atomic_num(z), mass(m), coh_sl(csl), inc_xs(ixs), cap_xs(cxs) {};
    };
    typedef std::map<const std::string, const IsotopeData> DataMap;
    static const NeutronSCL *instance();

    DataMap::const_iterator begin() const { return m_natural_elements.begin(); };
    DataMap::const_iterator end() const { return m_natural_elements.end(); };

    double getNeutronWeightedMass(const std::string& element) const;// mass divided by neutron mass
    double getAtomicMass(const std::string& element) const;//amu
    unsigned getAtomicNumber(const std::string& element) const;

    //in barn
    double getFreeXS(const std::string& element) const;
    double getBoundXS(const std::string& element) const;
    double getCaptureXS(const std::string& element) const;
    double getCoherentXS(const std::string& element) const;
    double getIncoherentXS(const std::string& element) const;

    //in sqrt(barn)
    double getIncoherentSL(const std::string& element) const;
    double getCoherentSL(const std::string& element) const;

    const std::string& getAtomName(unsigned z) const;

  private:
    NeutronSCL();
    ~NeutronSCL(){}

    NeutronSCL(NeutronSCL const&);
    void operator=(NeutronSCL const&);
    static NeutronSCL* m_instance;
    void addData(std::string name, unsigned z, double m, double csl, double ixs, double cxs);

    //TODO for NC2: use Z as key and argument in all member functions
    DataMap m_natural_elements;
    std::map<unsigned,std::string> m_z_to_elemname;

  };
}
#endif
