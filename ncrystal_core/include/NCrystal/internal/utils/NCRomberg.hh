#ifndef NCrystal_Romberg_hh
#define NCrystal_Romberg_hh

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

#include "NCrystal/core/NCDefs.hh"

namespace NCRYSTAL_NAMESPACE {

  class Romberg {
  public:

    //Efficient implementation of Romberg integration. Mostly for efficiency,
    //but also for robustness, it ALWAYS calculates at least R(4,4), which means
    //that at least 17 function evaluations will be required.
    //
    Romberg(){}
    virtual ~Romberg(){}

    //In the simplest use-case, override and implement your function in the evalFunc method:
    virtual double evalFunc(double) const = 0;

    //However, as integration will internally always evaluate the function at 16
    //or more equally spaced points at once, it might occasionally be more
    //efficient to evaluate the function at all of these points in a coherent
    //manner. In case this is desired, the user must override evalFuncMany and
    //evalFuncSum instead (due to the nature of C++, evalFunc will still have to
    //be implemented - however this can be a dummy implementation throwing an
    //exception in case of anyone calling it accidentally):

    virtual void evalFuncMany(double* fvals, unsigned n, double offset, double delta) const;
    virtual double evalFuncManySum(unsigned n, double offset, double delta) const;

    //Users should override this in order to change heuristics of when to stop
    //integration because the result is precise enough (or the level is getting
    //too deep and expensive).:
    virtual bool accept( unsigned level, double prev_estimate, double estimate,
                         double a, double b ) const;

    //Default behaviour in case of convertion issues is to throw an exception
    //after first dumping the function to a file for debugging. Client code can
    //override this method to change this behaviour (if the overriding methods
    //returns without exceptions, integration will return the best estimate
    //which may or may not be useful):
    virtual void convergenceError(double a, double b) const;

    //Perform integration (NB: integrand might be evaluated at b+epsilon where
    //epsilon represent a tiny numerical error):
    double integrate(double a, double b) const;


    //Dump the function to a file for inspection:
    void writeFctToFile(const std::string& filename, double a, double b, unsigned npts) const;


  };

}

#endif
