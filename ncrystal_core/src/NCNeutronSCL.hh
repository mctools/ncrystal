#ifndef NCrystal_NeutronSCL_hh
#define NCrystal_NeutronSCL_hh

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
      unsigned atomic_num;
      unsigned mass_num; //0 for natural element
      double mass;
      double coh_sl;  //in (barn)^0.5
      double inc_xs;  //in barn
      double cap_xs;  //in barn
      IsotopeData(unsigned z, unsigned a, double m, double csl, double ixs, double cxs)
      :atomic_num(z), mass_num(a), mass(m), coh_sl(csl), inc_xs(ixs), cap_xs(cxs) {};
    };

    static NeutronSCL *instance();

    std::map<const std::string, const IsotopeData>::const_iterator begin() { return m_natural_elements.begin(); };
    std::map<const std::string, const IsotopeData>::const_iterator end() { return m_natural_elements.end(); };

    double getNeutronWeightedMass(const std::string& element);
    double getAtomicMass(const std::string& element);
    unsigned getAtomicNumber(const std::string& element);

    //in barn
    double getFreeXS(const std::string& element);
    double getBoundXS(const std::string& element);
    double getCaptureXS(const std::string& element);
    double getCaptureXS_eV(const std::string& element, double kieV);
    double getCoherentXS(const std::string& element);
    double getIncoherentXS(const std::string& element);


    //in sqrt(barn)
    double getIncoherentSL(const std::string& element);
    double getCoherentSL(const std::string& element);

    private:
    NeutronSCL();
    ~NeutronSCL(){}

    NeutronSCL(NeutronSCL const&);
    void operator=(NeutronSCL const&);
    static NeutronSCL* m_instance;
    void addData(std::string name, unsigned z, unsigned a, double m, double csl, double ixs, double cxs);

    //TODO for NC2: use pair<z,a> as key, so all the member functions will use pair instead of std::string
    std::map<const std::string, const IsotopeData> m_natural_elements;

  };
}
#endif
