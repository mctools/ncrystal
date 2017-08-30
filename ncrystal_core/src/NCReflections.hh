#ifndef NCrystal_Reflections_hh
#define NCrystal_Reflections_hh

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

#include "NCVector.hh"
#include <vector>

//ReflectionFamily and Reflections hold parameters needed to describe Bragg
//diffraction physics. A reflection "family" is here defined as a group of
//planes with identical Fsquared and d-spacings.
//
//They are mainly used in the single crystal Bragg diffraction model.


namespace NCrystal {

  struct ReflectionFamily {

    std::vector<Vector> eqhkl_normals;
    double Fsqr;
    double inv2d;
    double wlthr;//cached wavelength cut, below which nothing in the family can contribute.
    double d;
    double h;
    double k;
    double l;

    ReflectionFamily(double ch, double ck, double cl, double fsqr, double dspacing, double wt = 0.0)
      : Fsqr(fsqr), inv2d(0.5/dspacing), wlthr(wt), d(dspacing), h(ch), k(ck), l(cl) {}

    ~ReflectionFamily() {}

#if __cplusplus >= 201103L
    //enable move assignment/construction:
    ReflectionFamily & operator= ( ReflectionFamily && ) = default;
    ReflectionFamily( ReflectionFamily && ) = default;
    //disable expensive assignment/copy construction:
    ReflectionFamily & operator= ( const ReflectionFamily & ) = delete;
    ReflectionFamily( const ReflectionFamily & ) = delete;
#else
    //Old C++ - have to implement more expensive assignment/copy construction:
    ReflectionFamily & operator= ( const ReflectionFamily & o )
    {
      eqhkl_normals = o.eqhkl_normals;//potentially expensive!
      Fsqr=o.Fsqr; d=o.d; inv2d=o.inv2d; h=o.h; k=o.k; l=o.l; wlthr=o.wlthr;
      return *this;
    }
    ReflectionFamily( const ReflectionFamily & o ) { *this = o; }
#endif
    bool operator < ( const ReflectionFamily & o ) const
    {
      //sort by d-spacing (fall-back to other variables for reproducibility):
      if ( o.d!=d ) return o.d < d;
      if ( o.Fsqr != Fsqr) return o.Fsqr < Fsqr;
      if ( o.h != Fsqr) return o.h < h;
      if ( o.k != Fsqr) return o.k < k;
      return o.l < l;
    }
  };

  struct Reflections{
    Reflections() : numAtom(0), V0(0) {}
    ~Reflections(){}
    std::vector<ReflectionFamily> family_list;
    double numAtom;
    double V0;
    void shrink_to_fit() {
      if (family_list.size()==family_list.capacity())
        return;
#if __cplusplus >= 201103L
      family_list.shrink_to_fit();//not 100% guarantee, but good enough hopefully.
#else
      //100% guarantee in both c++98 and c++11 (but we deleted copy constructors
      //in c++11 above to avoid mistakes elsewhere):
      std::vector<ReflectionFamily>(family_list.begin(),family_list.end()).swap(family_list);
#endif
    }
  };

}

#endif
