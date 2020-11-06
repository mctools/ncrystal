#ifndef NCrystal_AtomUtils_hh
#define NCrystal_AtomUtils_hh

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

  /////////////////////////////
  // Z<->elementName mapping //
  /////////////////////////////

  const std::string& elementZToName(unsigned z);//returns empty string if invalid
  unsigned elementNameToZ(const std::string&);//returns 0 if invalid

  ////////////////////////////
  // Decode atomic symbols. //
  ////////////////////////////

  class AtomSymbol final {
  public:
    //Valid atomic symbols are:
    //
    //1) Elements: "H", "He", "Li", ...
    //2) Isotopes: "H1", "H2", "He3", "Li6", ... (impossible ones with A<Z are invalid)
    //   NB: special aliases: "D" for "H2" and "T" for "H3".
    //3) Custom markers: "X", "X1", "X2", ..., "X99" (and no others).

    AtomSymbol(const std::string& s);
    ~AtomSymbol() = default;

    bool isInvalid() const;
    bool isCustomMarker() const;
    bool isElement() const;
    bool isIsotope() const;

    //If isElement() or isIsotope() it is safe to call (NB: A() returns 0 for
    //elements):
    unsigned Z() const;
    unsigned A() const;

  private:
    unsigned m_z = 0, m_a = 0;
    //(Z,0)=element, (Z,A>Z)=isotope, (0,1-100)=custom("X{i-1}") (0,0)=invalid
    void longInit(const std::string&);
  };

  /////////////////////////////////////////////////////////////////////////////////
  // Validate the syntex of atomdb lines.  This should be kept in sync with the  //
  // NCMAT data spec in ncmat_doc.md and the cfg-str description in NCMatCfg.hh. //
  // NCMatCfg.hh. The special syntax available/required in cfg-strings (with     //
  // semicolones and @) is not handled here. In case of problems, a BadInput     //
  // exception is thrown, which means that any code implementing subsequent      //
  // parsing of the lines is allowed to simply assume correct format. Note that  //
  // it only validates a single line at a time, and therefore does not verify    //
  // that a line with the word "nodefaults" comes at the beginning. It also does //
  // not verify if a line is missing previous definitions (i.e. if it refers to  //
  // X5 it can't check if X5 was defined previously):                            //
  /////////////////////////////////////////////////////////////////////////////////

  void validateAtomDBLine(const VectS& words, unsigned ncmat_version = 3);

}


////////////////////////////
// Inline implementations //
////////////////////////////

inline NCrystal::AtomSymbol::AtomSymbol(const std::string& s)
  : m_z(elementNameToZ(s)),
    m_a(0)
{
  //Usual case of regular element names is inlined, other cases needs
  //longer init:
  if (m_z==0)
    longInit(s);
}

inline bool NCrystal::AtomSymbol::isInvalid() const { return m_z == 0 && m_a ==0; }
inline bool NCrystal::AtomSymbol::isCustomMarker() const { return m_z == 0 && m_a !=0; }
inline bool NCrystal::AtomSymbol::isElement() const { return m_z != 0 && m_a ==0; }
inline bool NCrystal::AtomSymbol::isIsotope() const { return m_z != 0 && m_a !=0; }
inline unsigned NCrystal::AtomSymbol::Z() const { nc_assert(m_z != 0); return m_z; }
inline unsigned NCrystal::AtomSymbol::A() const { nc_assert(m_z != 0); return m_a; }

#endif
