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

#include "NCrystal/internal/utils/NCAtomUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMath.hh"
namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    static const std::string s_natelemlist[] = { "H"_s, "He"_s, "Li"_s, "Be"_s, "B"_s, "C"_s,
      "N"_s, "O"_s, "F"_s, "Ne"_s, "Na"_s, "Mg"_s, "Al"_s, "Si"_s, "P"_s, "S"_s, "Cl"_s,
      "Ar"_s, "K"_s, "Ca"_s, "Sc"_s, "Ti"_s, "V"_s, "Cr"_s, "Mn"_s, "Fe"_s, "Co"_s,
      "Ni"_s, "Cu"_s, "Zn"_s, "Ga"_s, "Ge"_s, "As"_s, "Se"_s, "Br"_s, "Kr"_s, "Rb"_s,
      "Sr"_s, "Y"_s, "Zr"_s, "Nb"_s, "Mo"_s, "Tc"_s, "Ru"_s, "Rh"_s, "Pd"_s, "Ag"_s,
      "Cd"_s, "In"_s, "Sn"_s, "Sb"_s, "Te"_s, "I"_s, "Xe"_s, "Cs"_s, "Ba"_s, "La"_s,
      "Ce"_s, "Pr"_s, "Nd"_s, "Pm"_s, "Sm"_s, "Eu"_s, "Gd"_s, "Tb"_s, "Dy"_s, "Ho"_s,
      "Er"_s, "Tm"_s, "Yb"_s, "Lu"_s, "Hf"_s, "Ta"_s, "W"_s, "Re"_s, "Os"_s, "Ir"_s,
      "Pt"_s, "Au"_s, "Hg"_s, "Tl"_s, "Pb"_s, "Bi"_s, "Po"_s, "At"_s, "Rn"_s, "Fr"_s,
      "Ra"_s, "Ac"_s, "Th"_s, "Pa"_s, "U"_s, "Np"_s, "Pu"_s, "Am"_s, "Cm"_s, "Bk"_s,
      "Cf"_s, "Es"_s, "Fm"_s, "Md"_s, "No"_s, "Lr"_s, "Rf"_s, "Db"_s, "Sg"_s, "Bh"_s,
      "Hs"_s, "Mt"_s, "Ds"_s, "Rg"_s, "Cn"_s, "Nh"_s, "Fl"_s, "Mc"_s, "Lv"_s, "Ts"_s,
      "Og"_s };
    static const std::map<std::string,unsigned> s_natelem_name2z_map = []()
    {
      std::map<std::string,unsigned> mm;
      constexpr unsigned zmax = (sizeof(s_natelemlist)/sizeof(s_natelemlist[0]));
      for (unsigned zm1 = 0; zm1<zmax; ++zm1)
        mm[s_natelemlist[zm1]] = zm1+1;
      return mm;
    }();

  }
}

const std::string& NC::elementZToName(unsigned z)
{
  constexpr unsigned zmax = (sizeof(s_natelemlist)/sizeof(s_natelemlist[0]));
  if ( z < 1 || z > zmax ) {
    static std::string str_empty;
    return str_empty;
  }
  return s_natelemlist[z-1];
}

unsigned NC::elementNameToZ(const std::string& name) {
  auto it = s_natelem_name2z_map.find(name);
  return ( it == s_natelem_name2z_map.end() ? 0 : it->second );
}

void NC::AtomSymbol::longInit(const std::string& symbol)
{
  //NB: Keep in synch with describeError above.
  nc_assert(m_z==0&&m_a==0);
  std::string tmp = symbol;
  trim(tmp);
  std::string lbl;
  std::string digits;
  std::tie(lbl,digits) = decomposeStrWithTrailingDigits( tmp );
  if ( digits.empty() ) {
    //Special cases first:
    if (lbl=="D") {
      m_z = 1;
      m_a = 2;
    } else if (lbl=="T") {
      m_z = 1;
      m_a = 3;
    } else if (lbl=="X") {
      m_a = 1;
    } else {
      m_z = elementNameToZ(lbl);
      //if m_z is 0, we are invalid, otherwise an element
    }
    return;
  }

  if ( digits.at(0) == '0' ) {
    //Error: Trailing digit can not start with 0!
    return;
  }
  unsigned digitval = static_cast<unsigned>(digits.size()>3 ? 9999 : str2int(digits));
  if (digitval>300) {
    //Error: Trailing digit too high
    return;
  }
  nc_assert(digitval>0);
  if (lbl=="X") {
    if (digitval>99)
      return;//Error: too high trailing digit for custom marker
    m_a = digitval+1;
    return;
  }
  m_z = elementNameToZ(lbl);
  if (m_z) {
    if ( digitval < m_z ) {
      m_z = 0;//Error: invalid isotope (has A < Z)
    } else {
      m_a = digitval;//isotope!
    }
    return;
  } else {
    //Error: Element not recognised (perhaps custom error for lbl=T/D?)
    return;
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    static const std::string s_lowerabc = "abcdefghijklmnopqrstuvwxyz";
    static const std::string s_upperabc = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  }
}

void NC::validateAtomDBLine(const VectS& words, unsigned ncmat_version )
{
  nc_assert(ncmat_version>=3);//for now all NCMAT versions >= v3 have the same
                              //ATOMDB format
  (void)ncmat_version;

  auto isDblWithUnit = [](const std::string& s, const std::string& unit, double& val) -> bool
  {
    if (!endswith(s,unit))
      return false;
    auto ss = s.substr(0,s.size()-unit.size());
    return safe_str2dbl(ss,val);
  };

  if (words.size()>10000)//not documented in ncmat_doc.md, but is really just a sanity check.
    NCRYSTAL_THROW2(BadInput,"Invalid specification. Too many words in line!");

  const unsigned nwords = static_cast<unsigned>(words.size());

  static const std::string allowed_special_chars="+-."_s;
  static const std::string allowed_chars = s_upperabc+s_lowerabc+"0123456789+-%*"_s+allowed_special_chars;
  for ( auto w : words ) {
    if (w.empty())
      NCRYSTAL_THROW(BadInput,"Invalid specification (empty part)");
    if (!isSimpleASCII(w,AllowTabs::No,AllowNewLine::No))
      NCRYSTAL_THROW2(BadInput,"Invalid specification (must only contain simple ascii characters) :\""<<w<<"\")");
    if (!contains_only(w,allowed_chars))
      NCRYSTAL_THROW2(BadInput,"Invalid specification (must only contain a-zA-Z0-9 and "<<allowed_special_chars<<") :\""<<w<<"\"");
  }

  if (nwords==1&&words.front()=="nodefaults")
    return;//ok (if first line)

  auto joinedline = [&words](){ return joinstr(words); };

  auto throwGeneralError = [&joinedline]() { NCRYSTAL_THROW2(BadInput,"Invalid specification. Line"
                                                             " has unsupported format: \""<<joinedline()<<"\""); };

  bool mode_is = (nwords==3||(nwords>=4&&nwords%2==0)) && words.at(1) == "is";
  double val_mass,val_csl, val_ixs,val_axs;
  bool mode_dataentry = ( nwords==5
                          && isDblWithUnit(words.at(1),"u"_s,val_mass)
                          && isDblWithUnit(words.at(2),"fm"_s,val_csl)
                          && isDblWithUnit(words.at(3),"b"_s,val_ixs)
                          && isDblWithUnit(words.at(4),"b"_s,val_axs) );
  if ( mode_dataentry ) {
    if ( ! (val_mass>0.0 && val_mass < 10000.0 ) )
      NCRYSTAL_THROW2(BadInput,"Invalid specification (mass value out of range in \""<<joinedline()<<"\")");
    if ( ! (val_csl>-1e7 && val_csl < 1e7 ) )
      NCRYSTAL_THROW2(BadInput,"Invalid specification (coherent scattering length out of range in \""<<joinedline()<<"\")");
    if ( ! (val_ixs>=0.0 && val_ixs < 1e7 ) )
      NCRYSTAL_THROW2(BadInput,"Invalid specification (incoherent cross section out of range in \""<<joinedline()<<"\")");
    if ( ! (val_axs>=0.0 && val_axs < 1e7 ) )
      NCRYSTAL_THROW2(BadInput,"Invalid specification (absorption cross section out of range in \""<<joinedline()<<"\")");
  }

  nc_assert(! (mode_is && mode_dataentry) );
  if ( !mode_is && !mode_dataentry )
    throwGeneralError();

  //Check label at first position:
  std::string label = words.at(0);
  AtomSymbol atomsymbol(label);
  auto invalidLabelError = [](const std::string&lbl) { NCRYSTAL_THROW2(BadInput,"Invalid specification. The label \""<<lbl
                                                  <<"\" is neither a standard element name (e.g. Al, H), an isotope"
                                                  " (e.g. Li6), or a custom marker (X, X1, X2, ..., X99)."); };
  if (atomsymbol.isInvalid())
    invalidLabelError(label);

  //Data entry lines:
  if ( mode_dataentry ) {
    if (atomsymbol.isCustomMarker() ) {
      NCRYSTAL_THROW2(BadInput,"Invalid specification (label \""<<label
                      <<"\" is a custom marker, but those are not allowed in data entry lines)");
    }
    return;
  }
  nc_assert(mode_is);
  if ( atomsymbol.isIsotope() ) {
    NCRYSTAL_THROW2(BadInput,"Invalid AtomDB specification (isotope labels like \""
                    <<label<<"\" are not allowed as alias or name of mixtures)");
  }
  if ( nwords==3 ) {
    //Lines like "X is Al"
    if ( AtomSymbol(words.at(2)).isInvalid() )
      invalidLabelError( words.at(2) );
  } else {
    //Lines like "X is 0.9 Al 0.1 Cr"
    nc_assert(nwords%2==0&&nwords>=4);
    unsigned ncomp = (nwords-2)/2;
    StableSum totfrac;
    for ( unsigned icomp = 0; icomp < ncomp; ++icomp ) {
      //Fraction:
      double frac;
      if (!safe_str2dbl(words.at(2+icomp*2),frac)||frac<=0.0 || !(frac<=1.0) )
        NCRYSTAL_THROW2(BadInput,"Invalid specification (invalid fraction: \""<<words.at(2+icomp*2)<<"\" )");
      totfrac.add(frac);
      if (AtomSymbol(words.at(3+icomp*2)).isInvalid())
        invalidLabelError( words.at(3+icomp*2) );
    }
    if ( ncabs(1.0-totfrac.sum())>1e-10 )//1e-10 here, 1e-9 in NCAtomDBExtender (which should be slightly more relaxed).
      NCRYSTAL_THROW2(BadInput,"Invalid specification (fractions do not add up to 1: \""<<joinedline()<<"\" )");
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    bool readNextChemFormEntry( const char*& it, const char* itE, DecodedChemForm& cf ) {
      //returns false in case of errors
      auto skipWhiteSpace = [](const char*& it2, const char* it2E)
      {
        while ( it2!=it2E && isWhiteSpace(*it2) )
          ++it2;
      };
      skipWhiteSpace(it,itE);
      if ( it == itE )
        return false;
      if (!contains(s_upperabc,*it))
        return false;//must start with upper case character.
      auto itSymbE = std::next(it);
      while ( itSymbE != itE && contains(s_lowerabc,*itSymbE) )
        ++itSymbE;//include all immediately following lower case chars in symbol

      AtomSymbol symbol(std::string(it,std::distance(it,itSymbE)));
      if (!symbol.isElement())
        return false;

      //skip chars used for symbol and any leading whitespace:
      it = itSymbE;
      skipWhiteSpace(it,itE);

      //See if we are followed by a number:
      if ( it==itE || contains(s_upperabc,*it ) ) {
        //end of string or encountered new element name, implicit count of 1:
        cf.emplace_back( 1, symbol );
        return true;
      }
      //Must be followed by a number here. Allowed forms are for now always
      //integral and positive, i.e. pure digits 0-9!
      auto isDigit = [](char c) { return c>='0'&&c<='9'; };
      std::uint_least32_t digitVal = 0;
      auto addDigitChar = [&digitVal](char c) {
        constexpr std::uint_least32_t maxAllowed = 1000000000;
        constexpr std::uint_least32_t maxAllowedDiv10 = maxAllowed/10;
        nc_assert(c>='0'&&c<='9');
        if ( digitVal > maxAllowedDiv10 )
          return false;
        digitVal *= 10;
        digitVal += (int(c)-int('0'));
        return digitVal <= maxAllowed;
      };
      if ( !isDigit(*it) )
        return false;
      if ( !addDigitChar(*it) )
        return false;
      //include all immediately following digits:
      auto itDigitE = std::next(it);
      while ( itDigitE != itE && isDigit(*itDigitE) ) {
        if ( !addDigitChar(*itDigitE) )
          return false;
        ++itDigitE;
      }
      if (!digitVal)
        return false;//prevent a count of 0
      cf.emplace_back( digitVal, symbol );
      it = itDigitE;
      skipWhiteSpace(it,itE);
      return true;
    }

  }
}

NC::Optional<NC::DecodedChemForm> NC::tryDecodeSimpleChemicalFormula( std::string s )
{
  DecodedChemForm cf;
  {
    const char * it = s.data();
    auto itE = it + s.size();
    while ( it != itE ) {
      if (!readNextChemFormEntry(it,itE,cf))
        return NullOpt;
    }
  }
  if ( cf.empty() )
    return NullOpt;
  if ( cf.size() == 1 )
    return cf;
  return normaliseSimpleChemicalFormula( cf );
}

NC::DecodedChemForm NC::decodeSimpleChemicalFormula( std::string s ) {
  auto cf = tryDecodeSimpleChemicalFormula(s);
  if (!cf.has_value())
    NCRYSTAL_THROW2(BadInput,"Invalid chemical formula: "<<s);
  return std::move(cf.value());
}

void NC::streamSimpleChemicalFormula( std::ostream& os, const DecodedChemForm& cf )
{
  nc_assert( !cf.empty() );
  for ( auto& e : cf ) {
    nc_assert(e.second.isElement());
    os << elementZToName(e.second.Z());
    if ( e.first != 1 )
      os << e.first;
  }
}

NC::DecodedChemForm NC::normaliseSimpleChemicalFormula( const DecodedChemForm& input )
{
  //Merge duplicates:
  NC::DecodedChemForm res;
  for ( const auto& e : input ) {
    nc_assert( e.second.isElement() );
    bool found = false;
    for ( auto& er : res ) {
      if ( er.second.Z() == e.second.Z() ) {
        er.first += e.first;
        found = true;
        break;
      }
    }
    if ( !found )
      res.push_back( e );
  }
  //Sort according to the hill system. First we must know if there are any
  //carbon atoms:
  const bool hasCarbon = [&res]()
  {
    for ( auto& e : res )
      if ( e.second.Z() == 6 )
        return true;
    return false;
  }();

  //Sort:
  auto sortValue = [&hasCarbon]( const AtomSymbol& atom ) -> const std::string&
  {
    nc_assert( atom.isElement() );
    const auto zval = atom.Z();
    const std::string& name = elementZToName(zval);
    nc_assert(!name.empty());

    if ( !hasCarbon || ( zval!=1 && zval!=6 ) )
      return name;
    static const std::string hillC("Aa");
    static const std::string hillH("Ab");
    nc_assert( hillC < hillH );
    nc_assert( hillC < "Ac"_s );
    nc_assert( "Ac"_s < "Ag"_s );
    nc_assert( "Ac"_s < "B"_s );
    return zval == 6 ? hillC : hillH;
  };
  std::stable_sort(res.begin(),res.end(),
                   [&sortValue]( const DecodedChemForm::value_type& a,
                                 const DecodedChemForm::value_type& b )
                   {
                     const std::string& sortA = sortValue(a.second);
                     const std::string& sortB = sortValue(b.second);
                     if ( sortA != sortB )
                       return sortA < sortB;
                     return a.first < b.first;
                   });
  return res;
}

