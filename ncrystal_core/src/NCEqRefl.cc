////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2021 NCrystal developers                                   //
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

#include "NCrystal/internal/NCEqRefl.hh"
#include "NCrystal/NCDefs.hh"

NCrystal::EqRefl::EqRefl(int sg)
{
  if (sg<1||sg>230)
    NCRYSTAL_THROW(BadInput,"Space group number is not in the range 1 to 230");

  if (sg<149) {
    if (sg<75) {
      if (sg<3)
        m_calc = &EqRefl::calc_Triclinic_1_2;
      else if (sg<16)
        m_calc = &EqRefl::calc_Monoclinic_3_15;
      else
        m_calc = &EqRefl::calc_Orthorhombic_16_74;
    } else {
      if (sg<89)
        m_calc = &EqRefl::calc_Tetragonal_75_88;
      else if (sg<143)
        m_calc = &EqRefl::calc_Tetragonal_89_142;
      else
        m_calc = &EqRefl::calc_Trigonal_143_148;
    }
  } else {
    if (sg<195) {
      if (sg<168)
        m_calc = &EqRefl::calc_Trigonal_149_167;
      else if (sg<177)
        m_calc = &EqRefl::calc_Hexagonal_168_176;
      else
        m_calc = &EqRefl::calc_Hexagonal_177_194;
    } else {
      if (sg<207)
        m_calc = &EqRefl::calc_Cubic_195_206;
      else
        m_calc = &EqRefl::calc_Cubic_207_230;
    }
  }

}

NCrystal::EqRefl::~EqRefl()
{
}

const std::set<NCrystal::EqRefl::HKL>& NCrystal::EqRefl::getEquivalentReflections(int h, int k, int l)
{
  m_planes.clear();
  add(h,k,l);//this one is common for all spacegroups
  (this->*m_calc)(h,k,l);
  return m_planes;
}

void NCrystal::EqRefl::calc_Triclinic_1_2(int, int, int)
{
}

void NCrystal::EqRefl::calc_Monoclinic_3_15(int h, int k, int l)
{
  add(h,-k,l);
}

void NCrystal::EqRefl::calc_Orthorhombic_16_74(int h, int k, int l)
{
  add(h,-k,-l); add(h,-k,l); add(h,k,-l);
}

void NCrystal::EqRefl::calc_Tetragonal_75_88(int h, int k, int l)
{
  add(k,-h,-l); add(h,k,-l); add(k,-h,l);
}

void NCrystal::EqRefl::calc_Tetragonal_89_142(int h, int k, int l)
{
  add(k,h,l); add(k,-h,-l); add(h,k,-l);
  add(k,-h,l); add(h,-k,-l); add(k,h,-l); add(h,-k,l);
}

void NCrystal::EqRefl::calc_Trigonal_143_148(int h, int k, int l)
{
  add(h+k,-h,-l); add(k,-h-k,l);
}

void NCrystal::EqRefl::calc_Trigonal_149_167(int h, int k, int l)
{
  add(h+k,-h,-l); add(k,-h-k,l); add(k,h,-l); add(h+k,-k,l); add(h,-h-k,-l);
}

void NCrystal::EqRefl::calc_Hexagonal_168_176(int h, int k, int l)
{
  add(h,k,-l); add(k,-h-k,-l); add(k,-h-k,l); add(h+k,-h,l); add(h+k,-h,-l);
}

void NCrystal::EqRefl::calc_Hexagonal_177_194(int h, int k, int l)
{
  add(k,-h-k,-l); add(h+k,-h,-l); add(h,k,-l); add(k,-h-k,l);
  add(h+k,-h,l); add(k,h,l); add(h+k,-k,l); add(h,-h-k,l);
  add(k,h,-l); add(h+k,-k,-l); add(h,-h-k,-l);
}

void NCrystal::EqRefl::calc_Cubic_195_206(int h, int k, int l)
{
  add(h,-k,-l); add(h,-k,l); add(h,k,-l);
  add(k,l,h); add(k,-l,-h); add(k,-l,h); add(k,l,-h);
  add(l,h,k); add(l,-h,-k); add(l,-h,k); add(l,h,-k);
}

void NCrystal::EqRefl::calc_Cubic_207_230(int h, int k, int l)
{
  add(h,-k,-l); add(h,-k,l); add(h,k,-l);
  add(k,l,h); add(k,-l,-h); add(k,-l,h); add(k,l,-h);
  add(l,h,k); add(l,-h,-k); add(l,-h,k); add(l,h,-k);
  add(k,h,l); add(k,-h,-l); add(k,-h,l); add(k,h,-l);
  add(l,k,h); add(l,-k,-h); add(l,-k,h); add(l,k,-h);
  add(h,l,k); add(h,-l,-k); add(h,-l,k); add(h,l,-k);
}
