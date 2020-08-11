#ifndef NCrystal_SABData_hh
#define NCrystal_SABData_hh

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

  class SABData : public UniqueID {
  public:

    //Immutable data structure defining an S(alpha,beta) scattering kernel.

    //Access data:
    const VectD& alphaGrid() const { return m_a; }
    const VectD& betaGrid() const { return m_b; }
    const VectD& sab() const { return m_sab; }
    double temperature() const { return m_t; }
    SigmaBound boundXS() const { return m_bxs; }
    double elementMassAMU() const { return m_m; }
    double suggestedEmax() const { return m_sem; }

    //Constructors etc. (all expensive operations forbidden):
    SABData( VectD&& alphaGrid, VectD&& betaGrid, VectD&& sab,
             double temperature, SigmaBound boundXS, double elementMassAMU,
             double suggestedEmax = 0 );
    SABData ( SABData && ) = default;
    SABData & operator= ( SABData && ) = default;
    SABData ( const SABData & ) = delete;
    SABData & operator= ( const SABData & ) = delete;
    SABData () = delete;
    ~SABData () = default;

  private:
    VectD m_a, m_b, m_sab;
    double m_t, m_m, m_sem;
    SigmaBound m_bxs;
  };

  class VDOSData : public UniqueID {
  public:

    //Immutable data structure defining a VDOS (phonon density spectrum).

    //Access data:
    const PairDD& vdos_egrid() const { return m_e; }
    const VectD& vdos_density() const { return m_d; }
    double temperature() const { return m_t; }
    SigmaBound boundXS() const { return m_bxs; }
    double elementMassAMU() const { return m_m; }

    //Constructors etc. (all expensive operations forbidden):
    VDOSData( PairDD egrid, VectD&& density,
              double temperature, SigmaBound boundXS, double elementMassAMU );
    VDOSData ( VDOSData && ) = default;
    VDOSData & operator= ( VDOSData && ) = default;
    VDOSData ( const VDOSData & ) = delete;
    VDOSData & operator= ( const VDOSData & ) = delete;
    VDOSData () = delete;
    ~VDOSData () = default;

  private:
    PairDD m_e;
    VectD m_d;
    double m_t, m_m;
    SigmaBound m_bxs;
  };

}


#endif
