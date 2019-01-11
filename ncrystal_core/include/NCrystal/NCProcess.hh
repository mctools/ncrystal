#ifndef NCrystal_Process_hh
#define NCrystal_Process_hh

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

#include "NCrystal/NCCalcBase.hh"

/////////////////////////////////////////////////////////////////////
// Base class for calculations of processes in materials.          //
//                                                                 //
// Note that for maximum compatibility vector directions are       //
// specified via double arrays of length 3.                        //
//                                                                 //
// Note that the unit for kinetic energy in the calls below is eV  //
// (electronvolt).                                                 //
/////////////////////////////////////////////////////////////////////

namespace NCrystal {

  class NCRYSTAL_API Process : public CalcBase {
  public:

    Process(const char * calculator_type_name);

    //The process cross-section (in barns):
    virtual double crossSection(double ekin, const double (&neutron_direction)[3] ) const = 0;

    //Override if the process is guaranteed to always produce vanishing
    //cross-sections outside a given domain (the special case
    //ekin_low=ekin_high=infinity is used to indicate a somewhat unusual process
    //with vanishing cross-section everywhere):
    virtual void domain(double& ekin_low, double& ekin_high) const { ekin_low = 0.0; ekin_high = infinity; }

    //Check if process is oriented. For non-oriented processes, the results do
    //not depend on the incident direction of the neutron, and outcomes such as
    //scatterings are phi-symmetric:
    virtual bool isOriented() const { return true; };

    //For non-oriented processes, callers can either use the method above to
    //access cross-sections (with free choice of reference frame for the
    //direction vectors), or use the following method for convenience:
    virtual double crossSectionNonOriented( double ekin ) const;

    virtual void validate();//call to perform a quick (incomplete) validation
                            //that cross sections are vanishing outside
                            //domain(..).

    bool isNull() const;
  protected:
    virtual ~Process();
  };

}

//Convenience function which allows most custom vector-classes (those without
//virtual methods and whose first 24 bytes holds x, y and z coordinates in three
//doubles) to be passed easily to the methods above by enclosing the variable
//name in NC_VECTOR_CAST(..) in the call, without the need for first copying
//values (nb: repeated in NCVector.hh for convenience):
#ifndef NC_VECTOR_CAST
#  define NC_VECTOR_CAST(v) (reinterpret_cast<double(&)[3]>(v))
#  define NC_CVECTOR_CAST(v) (reinterpret_cast<const double(&)[3]>(v))
#endif

#endif
