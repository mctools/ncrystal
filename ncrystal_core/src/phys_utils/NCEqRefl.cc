////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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

#include "NCrystal/internal/phys_utils/NCEqRefl.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
namespace NC = NCrystal;

struct NC::EqRefl::Helper {

  using ERL = EquivReflList;

  static ERL calc_Triclinic_1_2(int h, int k, int l)
  {
    ERL e; e.add(h,k,l);
    return e;
  }

  static ERL calc_Monoclinic_3_15(int h, int k, int l)
  {
    ERL e; e.add(h,k,l);
    e.add(h,-k,l);
    return e;
  }

  static ERL calc_Orthorhombic_16_74(int h, int k, int l)
  {
    ERL e; e.add(h,k,l);
    e.add(h,-k,-l); e.add(h,-k,l); e.add(h,k,-l);
    return e;
  }

  static ERL calc_Tetragonal_75_88(int h, int k, int l)
  {
    ERL e; e.add(h,k,l);
    e.add(k,-h,-l); e.add(h,k,-l); e.add(k,-h,l);
    return e;
  }

  static ERL calc_Tetragonal_89_142(int h, int k, int l)
  {
    ERL e; e.add(h,k,l);
    e.add(k,h,l); e.add(k,-h,-l); e.add(h,k,-l);
    e.add(k,-h,l); e.add(h,-k,-l); e.add(k,h,-l); e.add(h,-k,l);
    return e;
  }

  static ERL calc_Trigonal_143_148(int h, int k, int l)
  {
    ERL e; e.add(h,k,l);
    e.add(h+k,-h,-l); e.add(k,-h-k,l);
    return e;
  }

  static ERL calc_Trigonal_149_167(int h, int k, int l)
  {
    //Laue class -3m1 (and -3m for rhombohedral groups in hexagonal axes):
    ERL e; e.add(h,k,l);
    e.add(h+k,-h,-l); e.add(k,-h-k,l); e.add(k,h,-l); e.add(h+k,-k,l); e.add(h,-h-k,-l);
    return e;
  }

  static ERL calc_Trigonal_31m(int h, int k, int l)
  {
    //Laue class -31m (SG 149, 151, 153, 157, 159, 162, 163):
    ERL e; e.add(h,k,l);
    e.add(h+k,-h,-l); e.add(k,-h-k,l); e.add(k,h,l); e.add(h+k,-k,-l); e.add(h,-h-k,l);
    return e;
  }

  static ERL calc_Rhombohedral_3(int h, int k, int l)
  {
    //Laue class -3 in rhombohedral axes (SG 146, 148):
    ERL e; e.add(h,k,l);
    e.add(k,l,h); e.add(l,h,k);
    return e;
  }

  static ERL calc_Rhombohedral_3m(int h, int k, int l)
  {
    //Laue class -3m in rhombohedral axes (SG 155, 160, 161, 166, 167):
    ERL e; e.add(h,k,l);
    e.add(k,l,h); e.add(l,h,k); e.add(k,h,l); e.add(h,l,k); e.add(l,k,h);
    return e;
  }

  static bool isTrigonal_31m( int sg )
  {
    return ( sg==149 || sg==151 || sg==153 || sg==157 || sg==159
             || sg==162 || sg==163 );
  }

  static ERL calc_Hexagonal_168_176(int h, int k, int l)
  {
    ERL e; e.add(h,k,l);
    e.add(h,k,-l); e.add(k,-h-k,-l); e.add(k,-h-k,l); e.add(h+k,-h,l); e.add(h+k,-h,-l);
    return e;
  }

  static ERL calc_Hexagonal_177_194(int h, int k, int l)
  {
    ERL e; e.add(h,k,l);
    e.add(k,-h-k,-l); e.add(h+k,-h,-l); e.add(h,k,-l); e.add(k,-h-k,l);
    e.add(h+k,-h,l); e.add(k,h,l); e.add(h+k,-k,l); e.add(h,-h-k,l);
    e.add(k,h,-l); e.add(h+k,-k,-l); e.add(h,-h-k,-l);
    return e;
  }

  static ERL calc_Cubic_195_206(int h, int k, int l)
  {
    ERL e; e.add(h,k,l);
    e.add(h,-k,-l); e.add(h,-k,l); e.add(h,k,-l);
    e.add(k,l,h); e.add(k,-l,-h); e.add(k,-l,h); e.add(k,l,-h);
    e.add(l,h,k); e.add(l,-h,-k); e.add(l,-h,k); e.add(l,h,-k);
    return e;
  }

  static ERL calc_Cubic_207_230(int h, int k, int l)
  {
    ERL e; e.add(h,k,l);
    e.add(h,-k,-l); e.add(h,-k,l); e.add(h,k,-l);
    e.add(k,l,h); e.add(k,-l,-h); e.add(k,-l,h); e.add(k,l,-h);
    e.add(l,h,k); e.add(l,-h,-k); e.add(l,-h,k); e.add(l,h,-k);
    e.add(k,h,l); e.add(k,-h,-l); e.add(k,-h,l); e.add(k,h,-l);
    e.add(l,k,h); e.add(l,-k,-h); e.add(l,-k,h); e.add(l,k,-h);
    e.add(h,l,k); e.add(h,-l,-k); e.add(h,-l,k); e.add(h,l,-k);
    return e;
  }
};

NC::EqRefl::EqRefl( int sg, bool rhombohedral_axes )
{
  if (sg<1||sg>230)
    NCRYSTAL_THROW(BadInput,"Space group number is not in the range 1 to 230");
  if ( rhombohedral_axes ) {
    if ( !isRhombohedralSpaceGroup( sg ) )
      NCRYSTAL_THROW2(BadInput,"Rhombohedral axes requested for space group "
                      <<sg<<" which is not rhombohedral.");
    m_calc = ( ( sg == 146 || sg == 148 )
               ? &Helper::calc_Rhombohedral_3
               : &Helper::calc_Rhombohedral_3m );
    return;
  }
  if (sg<149) {
    if (sg<75) {
      if (sg<3)
        m_calc = &Helper::calc_Triclinic_1_2;
      else if (sg<16)
        m_calc = &Helper::calc_Monoclinic_3_15;
      else
        m_calc = &Helper::calc_Orthorhombic_16_74;
    } else {
      if (sg<89)
        m_calc = &Helper::calc_Tetragonal_75_88;
      else if (sg<143)
        m_calc = &Helper::calc_Tetragonal_89_142;
      else
        m_calc = &Helper::calc_Trigonal_143_148;
    }
  } else {
    if (sg<195) {
      if (sg<168)
        m_calc = ( Helper::isTrigonal_31m(sg)
                   ? &Helper::calc_Trigonal_31m
                   : &Helper::calc_Trigonal_149_167 );
      else if (sg<177)
        m_calc = &Helper::calc_Hexagonal_168_176;
      else
        m_calc = &Helper::calc_Hexagonal_177_194;
    } else {
      if (sg<207)
        m_calc = &Helper::calc_Cubic_195_206;
      else
        m_calc = &Helper::calc_Cubic_207_230;
    }
  }
}
