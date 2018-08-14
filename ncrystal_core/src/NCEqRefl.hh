#ifndef NCrystal_EqRefl_hh
#define NCrystal_EqRefl_hh

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

#include <set>

//Class EqRefl provides symmetry-equivalent reflections for a given space group
//number, by providing a list of all (h,k,l) indices symmetry-equivalent to a
//given (h,k,l) index.
//
//Note that while two equivalent (h,k,l) indices will always have identical
//values of d-spacing and structure factors (F^2), the opposite is not in
//general true. So (h,k,l) groupings based on (d,F^2) values might contain more
//than one symmetry-equivalency group in each (d,F^2) group.
//
//TODO for NC2: SCBragg currently assumes Info objects to either have
//demi-normals or contain equivalency groupings. This is currently true, but
//seems a bit fragile (e.g. the NCMAT factory produces (d,F^2) groups, not
//equivalency groups - but does fortunately provide demi-normals). A more robust
//situation would be desirable (either all factories should produce hkl lists
//grouped by equivalency, or SCBragg code should become more intelligent).

namespace NCrystal {

  class EqRefl
  {
  public:

    EqRefl(int spacegroup);
    ~EqRefl();

    struct HKL {
      HKL(int hh, int kk, int ll) : h(hh), k(kk), l(ll) {}
      int h;
      int k;
      int l;
      bool operator<(const HKL&o) const {
        return ( h!=o.h ? h<o.h : ( k!=o.k ? k<o.k : l<o.l ) );
      }
    };

    const std::set<HKL>& getEquivalentReflections(int h, int k, int l);

  private:
    void calc_Triclinic_1_2(int,int,int);
    void calc_Monoclinic_3_15(int,int,int);
    void calc_Orthorhombic_16_74(int,int,int);
    void calc_Tetragonal_75_88(int,int,int);
    void calc_Tetragonal_89_142(int,int,int);
    void calc_Trigonal_143_148(int,int,int);
    void calc_Trigonal_149_167(int,int,int);
    void calc_Hexagonal_168_176(int,int,int);
    void calc_Hexagonal_177_194(int,int,int);
    void calc_Cubic_195_206(int,int,int);
    void calc_Cubic_207_230(int,int,int);

    void add(int h,int k,int l) {
      HKL a(h,k,l);
      HKL am(-h,-k,-l);
      //only insert one deminormal, not both am and a:
      if (!m_planes.count(a)&&!m_planes.count(am))
        m_planes.insert(am<a?a:am);
    }

    std::set<HKL> m_planes;
    void (EqRefl::*m_calc) (int,int,int);

  };


}

#endif
