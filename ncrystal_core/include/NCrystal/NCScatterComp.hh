#ifndef NCrystal_ScatterComp_hh
#define NCrystal_ScatterComp_hh

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

#include "NCrystal/NCScatter.hh"
#include <vector>
#include <utility>

/////////////////////////////////////////////////////////////////////
// Composition class which combines a list of scatter calculators  //
// into one, more complete picture.                                //
//                                                                 //
/////////////////////////////////////////////////////////////////////

namespace NCrystal {

  class ScatterComp : public Scatter {
  public:

    ScatterComp(const char * calculator_type_name = "ScatterComp");

    //Must add at least one component:
    void addComponent(Scatter*, double scale = 1.0 );

    size_t nComponents() const { return m_calcs.size(); }
    const Scatter * component(size_t i) const { return m_calcs.at(i).scatter; }
    Scatter * component(size_t i) { return m_calcs.at(i).scatter; }

    virtual double crossSection(double ekin, const double (&neutron_direction)[3] ) const;

    virtual void domain(double& ekin_low, double& ekin_high) const {
      ekin_low = m_threshold_lower; ekin_high = m_threshold_upper;
    }

    virtual void generateScattering( double ekin, const double (&neutron_direction)[3],
                                     double (&resulting_neutron_direction)[3], double& delta_ekin ) const;

    virtual bool isOriented() const;

  protected:
    virtual ~ScatterComp();
    struct Component {
      double threshold_lower;
      double threshold_upper;
      double scale;
      Scatter* scatter;
      bool operator<(const Component&) const;
    };
    std::vector<Component> m_calcs;//thresholds and calcs
    double m_threshold_lower;
    double m_threshold_upper;
    mutable int m_isOriented;
    void checkIsOriented() const;
  };

}

#endif
