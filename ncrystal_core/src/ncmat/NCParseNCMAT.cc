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

#include "NCrystal/internal/ncmat/NCParseNCMAT.hh"
#include "NCrystal/core/NCException.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <sstream>
#if nc_cplusplus >= 201703L
#  include <functional>//for std::invoke
#else
#  define NCPARSENCMAT_CALL_MEMBER_FN(objectaddr,ptrToMember)  ((*objectaddr).*(ptrToMember))
#endif
namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {

  namespace {

    constexpr bool highBitIsSet( char c ) noexcept
    {
      //Remember char might be signed or unsigned! So this might mean a value
      //>=128, or it might mean a negative value.
      return c & '\x80';
    }
    static_assert( highBitIsSet( (char)(128) ), "" );
    static_assert( highBitIsSet( (char)(255) ), "" );
    static_assert( highBitIsSet( (char)(140) ), "" );
    static_assert( highBitIsSet( (char)(221) ), "" );
    static_assert( !highBitIsSet( (char)(0) ), "" );
    static_assert( !highBitIsSet( (char)(127) ), "" );
    static_assert( !highBitIsSet( (char)(1) ), "" );
    static_assert( !highBitIsSet( (char)(32) ), "" );
    static_assert( !highBitIsSet( (char)(90) ), "" );
    static_assert( !highBitIsSet( '\x7F' ), "" );
    static_assert( highBitIsSet( '\xFF' ), "" );
  }

  class NCMATParser {
  public:

    //Parse .ncmat files.

    //Parse input. Will throw BadInput exceptions in case of problems. It will
    //do some rudimentary syntax checking (including presence/absence of data
    //sections and will call NCMATData::validate), but not a full validation of
    //data (a more complete validation is typically carried out afterwards by
    //the NCMAT Loader code). It will always clear the input pointer
    //(i.e. release/close the resource).
    NCMATParser(const TextData&);
    ~NCMATParser() = default;

    NCMATData&& getData() { return std::move(m_data); }

  private:

    typedef VectS Parts;
    void parseFile( TextData::Iterator itLine, TextData::Iterator itLineE );
    void parseLine( const std::string&, Parts&, unsigned linenumber ) const;
    void validateElementName(const std::string& s, unsigned lineno) const;
    double str2dbl_withfractions(const std::string&) const;

    //Section handling:
    typedef void (NCMATParser::*handleSectionDataFn)(const Parts&,unsigned);
    void handleSectionData_HEAD(const Parts&,unsigned);
    void handleSectionData_CELL(const Parts&,unsigned);
    void handleSectionData_ATOMPOSITIONS(const Parts&,unsigned);
    void handleSectionData_SPACEGROUP(const Parts&,unsigned);
    void handleSectionData_DEBYETEMPERATURE(const Parts&,unsigned);
    void handleSectionData_DYNINFO(const Parts&,unsigned);
    void handleSectionData_DENSITY(const Parts&,unsigned);
    void handleSectionData_ATOMDB(const Parts&,unsigned);
    void handleSectionData_STATEOFMATTER(const Parts&,unsigned);
    void handleSectionData_CUSTOM(const Parts&,unsigned);
    void handleSectionData_OTHERPHASES(const Parts&,unsigned);
    void handleSectionData_TEMPERATURE(const Parts&,unsigned);

    //Collected data:
    NCMATData m_data;

    //Long vectors of data in @DYNINFO sections are kept in:
    NCMATData::DynInfo * m_active_dyninfo;
    VectD * m_dyninfo_active_vector_field;
    bool m_dyninfo_active_vector_field_allownegative;

    //Handle "cubic" keyword in @CELL section:
    Optional<double> m_cell_cubic;

    //For error messages:
    std::string descr() const
    {
      std::string s;
      s.reserve(m_data.sourceDescription.str().size()+2);
      s += '"';
      s+= m_data.sourceDescription.str();
      s+='"';
      return s;
    }
  };

  NCMATData parseNCMATData( const TextData& text, bool doFinalValidation )
  {
    NCMATParser parser( text );
    if (!doFinalValidation)
      return parser.getData();
    NCMATData data = parser.getData();
    data.validate();
    return data;
  }

}

double NC::NCMATParser::str2dbl_withfractions(const std::string& ss) const
{
if (!contains(ss,'/'))
  return str2dbl(ss);
 if (m_data.version==1)
   NCRYSTAL_THROW2(BadInput,"specification with fractions not supported in"
                   " NCMAT v1 files (offending parameter is \""<<ss<<"\")");

 VectS parts;
 split(parts,ss,0,'/');
 if (parts.size()!=2)
   NCRYSTAL_THROW2(BadInput,"multiple fractions in numbers are not supported so could not parse \""<<ss<<"\"");
 for (auto&e: parts)
   if (e.empty())
     NCRYSTAL_THROW2(BadInput,"empty denominator or numerator so could not parse \""<<ss<<"\"");
 double a = str2dbl(parts.at(0));
 double b = str2dbl(parts.at(1));
 if (ncisnan(a)||ncisnan(b)||ncisinf(a)||ncisinf(b))
   NCRYSTAL_THROW2(BadInput,"invalid division attempted in \""<<ss<<"\"");
 if (!b)
   NCRYSTAL_THROW2(BadInput,"division by zero attempted in \""<<ss<<"\"");
 return a/b;
}

NC::NCMATParser::NCMATParser( const TextData& input )
  : m_active_dyninfo(nullptr),
    m_dyninfo_active_vector_field(nullptr),
    m_dyninfo_active_vector_field_allownegative(false)
{
  //Setup source description strings first as it is used in error messages:
  m_data.sourceDescription = input.dataSourceName();

  //Inspect first line to ensure format is NCMAT and extract version:
  auto itLine = input.begin();

  if ( itLine == input.end() )
    NCRYSTAL_THROW2(BadInput,"Empty data: "<<descr());
  const std::string& line = *itLine;

  //First line is special, we want the file to start with "NCMAT" with no
  //whitespace in front, so we explicitly test this before invoking the more
  //generic parseLine machinery below:
  if (!startswith(line,"NCMAT"))
    NCRYSTAL_THROW2(BadInput,descr()<<": is not in NCMAT format: The first 5 characters in the first line must be \"NCMAT\"");

  //Parse first line to get file format version:
  Parts parts;
  parseLine(line,parts,1);
  if ( parts.size() == 2 ) {
    if ( parts.at(1) == "v1" ) {
      m_data.version = 1;
      if (contains(line,'#'))
        NCRYSTAL_THROW2(BadInput,descr()<<": has comments in the first line, which is not allowed in the NCMAT v1 format");
    } else if ( parts.at(1) == "v2" ) {
      m_data.version = 2;
    } else if ( parts.at(1) == "v3" ) {
      m_data.version = 3;
    } else if ( parts.at(1) == "v4" ) {
      m_data.version = 4;
    } else if ( parts.at(1) == "v5" ) {
      m_data.version = 5;
    } else if ( parts.at(1) == "v6" ) {
      m_data.version = 6;
    } else if ( parts.at(1) == "v7" ) {
      m_data.version = 7;
    } else {
      NCRYSTAL_THROW2(BadInput,descr()<<": is in an NCMAT format version, \""<<parts.at(1)<<"\", which is not recognised by this installation of NCrystal");
      static_assert( supported_ncmat_format_version_max == 7, "");
    }
  }
  if (!m_data.version)
    NCRYSTAL_THROW2(BadInput,descr()<<": is missing clear NCMAT format version designation in the first line, which should look like e.g. \"NCMAT v1\".");

  //Initial song and dance to classify source and format is now done, so proceed to parse rest of file:
  parseFile( ++itLine, input.end() );

  //Unalias element names:
  m_data.unaliasElementNames();

  //Check that spacegroup is absent or 195-230 if cubic keyword was used (we
  //can't do this in NCMatData::validate since it doesn't know if cubic keyword
  //was used, or lengths/angles just happen to look cubic - because the user
  //should in principle be allowed to specify a lower compatible symmetry
  //(e.g. spacegroup=1) if not actually using the "cubic" keyword:
  if ( m_cell_cubic.has_value() && m_data.hasSpaceGroup() && m_data.spacegroup<195 )
    NCRYSTAL_THROW2(BadInput,descr()<<": The \"cubic\" keyword is not allowed in the @CELL section"
                    " if the @SPACEGROUP section indicates a non-cubic group (cubic space group numbers are"
                    " 195..230 which does not include the provided number: "<<m_data.spacegroup<<")");

}

void NC::NCMATParser::parseFile( TextData::Iterator itLine, TextData::Iterator itLineE )
{
  //Setup map which will be used to delegate parsing of individual sections (NB:
  //Must also update error reporting code below when adding/remove section names
  //in new format versions):
  std::map<std::string,handleSectionDataFn> section2handler;
  section2handler["HEAD"] = &NCMATParser::handleSectionData_HEAD;
  section2handler["CELL"] = &NCMATParser::handleSectionData_CELL;
  section2handler["ATOMPOSITIONS"] = &NCMATParser::handleSectionData_ATOMPOSITIONS;
  section2handler["SPACEGROUP"] = &NCMATParser::handleSectionData_SPACEGROUP;
  section2handler["DEBYETEMPERATURE"] = &NCMATParser::handleSectionData_DEBYETEMPERATURE;
  if (m_data.version>=2) {
    section2handler["DYNINFO"] = &NCMATParser::handleSectionData_DYNINFO;
    section2handler["DENSITY"] = &NCMATParser::handleSectionData_DENSITY;
  }
  if (m_data.version>=3) {
    section2handler["ATOMDB"] = &NCMATParser::handleSectionData_ATOMDB;
    section2handler["CUSTOM"] = &NCMATParser::handleSectionData_CUSTOM;
  }
  if (m_data.version>=5) {
    section2handler["STATEOFMATTER"] = &NCMATParser::handleSectionData_STATEOFMATTER;
  }
  if (m_data.version>=6) {
    section2handler["OTHERPHASES"] = &NCMATParser::handleSectionData_OTHERPHASES;
  }
  if (m_data.version>=7) {
    section2handler["TEMPERATURE"] = &NCMATParser::handleSectionData_TEMPERATURE;
  }

  //Technically handle the part before the first section ("@SECTIONNAME") by the
  //same code as all other parts of the file, by putting it in a "HEAD" section:
  std::string current_section = "HEAD";
  std::map<std::string,handleSectionDataFn>::const_iterator itSection = section2handler.find(current_section);
  std::set<std::string> sections_seen;

  //Actual parsing (starting from the second line of input in this function):
  unsigned lineno(1);
  Parts parts;
  parts.reserve(16);

  bool sawAnySection = false;
  for ( ; itLine != itLineE; ++itLine ) {
    const std::string& line = *itLine;
    parseLine(line,parts,++lineno);

    if (m_data.version==1 && contains(line,'#')) {
      if (sawAnySection||(!parts.empty()&&parts.at(0)[0]=='@')||line.at(0)!='#')
        NCRYSTAL_THROW2(BadInput,descr()<<": has comments in a place which "
                        "is not allowed in the NCMAT v1 format"" (must only appear "
                        "before the first data section and with the # marker at the"
                        " beginning of the line).");
    }

    //ignore lines which are empty or only whitespace and comments:
    if (parts.empty())
      continue;

    if (parts.at(0)[0]=='@') {
      //New section marker! First check that the syntax of this line is valid:
      sawAnySection = true;
      if (parts.size()>1)
        NCRYSTAL_THROW2(BadInput,descr()<<": should not have non-comment entries after a section marker"
                        " (found \""<<parts.at(1)<<"\" after \""<<parts.at(0)<<"\" in line "<<lineno<<")");
      if (line[0]!='@')
        NCRYSTAL_THROW2(BadInput,descr()<<": should not have whitespace before a section marker"
                        " (problem with indented \""<<parts.at(0)<<"\" in line "<<lineno<<")");

      std::string new_section(&parts.at(0)[1]);
      if (new_section.empty())
        NCRYSTAL_THROW2(BadInput,descr()<<": has missing section name after '@' symbol in line "<<lineno<<")");

      const bool is_custom_section = startswith(new_section,"CUSTOM_");

      //Close current section by sending it an empty parts list (and previous line number where that section ended).
      parts.clear();
#if nc_cplusplus >= 201703L
      std::invoke(itSection->second,*this,parts,lineno-1);
#else
      NCPARSENCMAT_CALL_MEMBER_FN(this,itSection->second)(parts,lineno-1);
#endif

      //Guard against repeating an existing section (unless DYNINFO or custom sections, where it is allowed)
      bool multiple_sections_allowed = ( is_custom_section || new_section=="DYNINFO" );
      if ( !multiple_sections_allowed ) {
        if (sections_seen.count(new_section) )
          NCRYSTAL_THROW2(BadInput,descr()<<": multiple @"<<new_section<<" sections are not allowed (line "<<lineno<<")");
        sections_seen.insert(new_section);
      }

      //Try to switch to new section
      std::swap(current_section,new_section);
      itSection = section2handler.find( is_custom_section ? "CUSTOM"_s : current_section );

      nc_assert( m_data.version>=1 && m_data.version <= 7 );
      if ( itSection == section2handler.end() ) {
        //Unsupported section name. For better error messages, first check if it
        //is due to file version:
        if (m_data.version==1 && (current_section == "DYNINFO"||current_section=="DENSITY") ) {
          NCRYSTAL_THROW2(BadInput,descr()<<": has @"<<current_section<<" section which is not supported in the indicated"
                          " NCMAT format version, \"NCMAT v1\". It is only available starting with \"NCMAT v2\".");
        }
        if ( m_data.version<3 && (is_custom_section||current_section=="ATOMDB") ) {
          NCRYSTAL_THROW2(BadInput,descr()<<": has @"<<current_section<<" section which is not supported in the indicated"
                          " NCMAT format version, \"NCMAT v"<<m_data.version<<"\". It is only available starting with \"NCMAT v3\".");
        }
        if ( m_data.version<5 && (current_section=="STATEOFMATTER") ) {
          NCRYSTAL_THROW2(BadInput,descr()<<": has @STATEOFMATTER section which is not supported in the indicated"
                          " NCMAT format version, \"NCMAT v"<<m_data.version<<"\". It is only available starting with \"NCMAT v5\".");
        }
        if ( m_data.version<6 && (current_section=="OTHERPHASES") ) {
          NCRYSTAL_THROW2(BadInput,descr()<<": has @OTHERPHASES section which is not supported in the indicated"
                          " NCMAT format version, \"NCMAT v"<<m_data.version<<"\". It is only available starting with \"NCMAT v6\".");
        }
        if ( m_data.version<7 && (current_section=="TEMPERATURE") ) {
          NCRYSTAL_THROW2(BadInput,descr()<<": has @TEMPERATURE section which is not supported in the indicated"
                          " NCMAT format version, \"NCMAT v"<<m_data.version<<"\". It is only available starting with \"NCMAT v7\".");
        }
        NCRYSTAL_THROW2(BadInput,descr()<<": has @"<<current_section<<" section which is not a supported section name.");
      }
      //Succesfully switched to the new section, proceed to next line (after
      //adding entry in customSections in case of a custom section):
      if ( is_custom_section ) {
        if ( current_section.size() <= 7 )
          NCRYSTAL_THROW2(BadInput,descr()<<": has @"<<current_section
                          <<" section (needs additional characters after \"CUSTOM_\").");
        m_data.customSections.emplace_back(current_section.substr(7),NCMATData::CustomSectionData());
      }
      continue;
    }

    //Line inside active section was succesfully parsed.
#if nc_cplusplus >= 201703L
    std::invoke(itSection->second,*this,parts,lineno);
#else
    NCPARSENCMAT_CALL_MEMBER_FN(this,itSection->second)(parts,lineno);
#endif
  }

  //End of input. Close current section by sending it an empty parts list.
  parts.clear();
#if nc_cplusplus >= 201703L
  std::invoke(itSection->second,*this,parts,lineno);
#else
  NCPARSENCMAT_CALL_MEMBER_FN(this,itSection->second)(parts,lineno);
#endif

}


void NC::NCMATParser::parseLine( const std::string& line,
                                 Parts& parts,
                                 unsigned lineno ) const
{
  //Ignore trailing comments and split line on all whitespace to return the
  //actual parts in a vector. This function is a bit like
  //line.split('#',1)[0].split() in Python, but also checks encoding which is
  //different for comments and outside comments.

  //We only allow pure ASCII in the non-comment parts, and UTF-8 in
  //comments. For the ASCII parts, we don't allow any control characters except
  //\n\r\t. For efficiency reasons, we don't currently check that the comments
  //are actually UTF-8, but we could in principle perform a few checks easily:
  //UTF-8 strings never have null bytes, and any multi-byte character will have
  //all bytes with values >=128 (i.e. the high bit is set in all bits of
  //multi-byte characters).

  //Relevant parts of ASCII byte values for our purposes:
  //
  //0-31 forbidden control chars, except \t (9), \n (10), \r (13)
  //32 : space
  //33 : ! (valid char)
  //34 : " (valid char)
  //35 : # (valid char)
  //36-126 : all valid chars
  //127: forbidden control char.
  //negative or 128-255: forbidden (but could indicate UTF-8 multibyte char).

  parts.clear();
  const char * c = &line[0];
  const char * cE = c + line.size();
  const char * partbegin = nullptr;
  for (;c!=cE;++c) {
    if ( *c > 32 && *c < 127 && *c != '#') {
      //A regular character which should go in the parts vector
      if (!partbegin)
        partbegin = c;
      continue;
    }
    if ( isOneOf(*c,' ','\t') ) {
      //A whitespace character (we don't support silly stuff like vertical tabs,
      //and we kind of only grudgingly and silent accept tabs as well)
      if (partbegin) {
        parts.emplace_back(partbegin,c-partbegin);
        partbegin=0;
      }
      continue;
    }
    if ( isOneOf(*c,'\n','\r','#') ) {
      //EOL or comment begin. Only allow \r if it is in \r\n combination
      //(e.g. "DOS line endings"). A standalone \r can hide the line leading up
      //to it in printouts:
      //
      // #include <iostream>
      // int main() {
      //   std::cout<< "@CELL \n\r#comment\n..."<<std::endl<<std::endl;
      //   std::cout<< "@CELL \r\n#comment\n..."<<std::endl<<std::endl;
      //   std::cout<< "@CELL \n#comment\n..."<<std::endl<<std::endl;
      //   std::cout<< "@CELL \r#comment\n..."<<std::endl<<std::endl;
      // }
      //
      //Gives the output when run in a terminal (notice the missing @CELL in the fourth case!):
      //
      //@CELL
      //#comment
      //...
      //
      //@CELL
      //#comment
      //...
      //
      //@CELL
      //#comment
      //...
      //
      //#comment
      //...
      //
      // -> we don't really care that '\r' was used in ancient Macintosh systems...)
      // -> we also live with the fact that '\r' could still hide within comments
      //    (just too inefficient to search all comments)
      //
      if (*c=='\r') {
        if ( (c+1)!=cE && *(c+1)!='\n' ) {
          NCRYSTAL_THROW2(BadInput,descr()<<": contains invalid character at position "
                          <<(c-&line[0])<<" in line "<<lineno<<". Carriage return codes (aka \\r) "
                          " are not allowed unless used as part of DOS line endings.");
        }
      }

      break;
    }
    //Only reach here in case of errors:
    NCRYSTAL_THROW2(BadInput,descr()<<": contains invalid character at position "
                    <<(c-&line[0])<<" in line "<<lineno<<". Only regular ASCII characters"
                    " (including spaces) are allowed outside comments (comments can be UTF-8)");
  }
  if (partbegin) {
    //still need to add last part
    parts.emplace_back(partbegin,c-partbegin);
    partbegin=0;
  }

  //Check no illegal control codes occur in comments:
  for (;c!=cE;++c) {
    if ( *c>=32 && *c!=127 )
      continue;//ascii printable character

    if ( highBitIsSet( *c ) )
      continue;//possibly part of a multibyte utf-8 char, and not a control
               //char. More advanced utf-8 analysis needs rather complicated
               //code (we could do it of course...).

    if ( isOneOf(*c,'\t','\n') )
      continue;//ok

    if (*c=='\r') {
      if ( (c+1)!=cE && *(c+1)!='\n' ) {
        NCRYSTAL_THROW2(BadInput,descr()<<": contains invalid character at position "
                        <<(c-&line[0])<<" in line "<<lineno<<". Carriage return codes (aka \\r) "
                        " are not allowed unless used as part of DOS line endings.");
      }
      continue;
    }
    NCRYSTAL_THROW2(BadInput,descr()<<": contains illegal control code character in line "<<lineno);
  }
}

void NC::NCMATParser::handleSectionData_HEAD(const Parts& parts, unsigned lineno)
{
  if (parts.empty())
    return;
  //The HEAD pseudo-section should not have any actual contents.
  NCRYSTAL_THROW2(BadInput,descr()<<": should not have non-comment entries before the"
                  " first section (found \""<<parts.at(0)<<"\" in line "<<lineno<<")");
}

void NC::NCMATParser::handleSectionData_CELL(const Parts& parts, unsigned lineno)
{
  if (parts.empty()) {
    //finish up, apply "cubic" and validate.
    if ( m_cell_cubic.has_value() ) {
      nc_assert( m_data.cell.lengths[0]==0.0 && m_data.cell.lengths[1]==0.0 && m_data.cell.lengths[2]==0.0 );
      nc_assert( m_data.cell.angles[0]==0.0 && m_data.cell.angles[1]==0.0 && m_data.cell.angles[2]==0.0 );
      m_data.cell.lengths = { m_cell_cubic.value(),m_cell_cubic.value(),m_cell_cubic.value() };
      m_data.cell.angles = { 90.0, 90.0, 90.0 };
    }
    try {
      m_data.validateCell();
    } catch (Error::BadInput&e) {
      NCRYSTAL_THROW2(BadInput,e.what()<<" (problem in the @CELL section ending in line "<<lineno<<")");
    }
    return;
  }
  const auto& keyword = parts.at(0);

  if ( keyword=="cubic" ) {
    if ( m_data.version < 4 )
      NCRYSTAL_THROW2(BadInput,descr()<<": \"cubic\" keyword in @CELL section requires NCMAT v4 or later. Problem in line "<<lineno);
    const bool has_lengths = !( m_data.cell.lengths[0]==0.0 && m_data.cell.lengths[1]==0.0 && m_data.cell.lengths[2]==0.0 );
    const bool has_angles = !( m_data.cell.angles[0]==0.0 && m_data.cell.angles[1]==0.0 && m_data.cell.angles[2]==0.0 );
    if ( has_lengths || has_angles ) {
      NCRYSTAL_THROW2(BadInput,descr()<<": The \"cubic\" keyword can not be provided at the same time as the \""
                      <<(has_lengths?"lengths":"angles")
                      <<"\" keyword in the @CELL section in line "<<lineno);
    }
    if ( m_cell_cubic.has_value() )
      NCRYSTAL_THROW2(BadInput,descr()<<": repeated keyword \"cubic\" in line "<<lineno);
    if ( parts.size() != 2 )
      NCRYSTAL_THROW2(BadInput,descr()<<": wrong number of data entries after \"cubic\" keyword in line "<<lineno<<" (expected a single number)");
    try {
      double cubic_val = str2dbl(parts.at(1));
      m_cell_cubic = cubic_val;
      if ( !(cubic_val>0.0) || cubic_val>1e4 )
        NCRYSTAL_THROW(BadInput,"invalid value or value out of range");
    } catch (Error::BadInput&e) {
      NCRYSTAL_THROW2(BadInput,descr()<<": problem while decoding \"cubic\" parameter in line "<<lineno<<" : "<<e.what());
    }
    return;
  }
  if ( !isOneOf(keyword,"lengths","angles") ) {
    NCRYSTAL_THROW2(BadInput,descr()<<": found \""<<keyword<<"\" where \"lengths\""
                    <<(m_data.version>=4?",  \"angles\", or  \"cubic\"":" or \"angles\"")
                    <<" keyword was expected in @CELL section in line "<<lineno);
  }
  if ( parts.size() != 4 ) {
    NCRYSTAL_THROW2(BadInput,descr()<<": wrong number of data entries after \""<<keyword<<"\" keyword in line "<<lineno<<" (expected three numbers)");
  }
  std::array<double,3>& targetvector = ( keyword=="lengths" ? m_data.cell.lengths : m_data.cell.angles );
  if (!(targetvector[0]==0.&&targetvector[1]==0.&&targetvector[2]==0.)) {
    NCRYSTAL_THROW2(BadInput,descr()<<": repeated keyword \""<<keyword<<"\" in line "<<lineno);
  }
  std::array<double,3> v;
  for (unsigned i = 0; i<3; ++i) {
    if ( parts.at(i+1) == "!!" ) {
      if ( keyword!="lengths" )
        NCRYSTAL_THROW2(BadInput,descr()<<": Usage of \"!!\" to repeat previous value can only be used for \"lengths\" keyword, not \""<<keyword<<"\" (in line "<<lineno<<")");
      if ( i == 0 )
        NCRYSTAL_THROW2(BadInput,descr()<<": Usage of \"!!\" to repeat previous length value can not be used for the first value (in line "<<lineno<<")");
      if ( m_data.version < 4 )
        NCRYSTAL_THROW2(BadInput,descr()<<": Usage of \"!!\" to repeat previous length value requires NCMAT v4 or later (in line "<<lineno<<")");
      v.at(i) = v.at(i-1);
      continue;
    }
    try {
      v[i] = str2dbl(parts.at(i+1));
    } catch (Error::BadInput&e) {
      NCRYSTAL_THROW2(BadInput,descr()<<": problem while decoding \""<<keyword<<"\" parameter #"<<i+1<<" in line "<<lineno<<" : "<<e.what());
    }
  }
  targetvector = v;
  if ( targetvector[0]==0. && targetvector[1]==0. && targetvector[2]==0. ) {
    NCRYSTAL_THROW2(BadInput,descr()<<": vector \""<<keyword<<"\" is a null-vector in line "<<lineno);
  }
}

void NC::NCMATParser::validateElementName(const std::string& s, unsigned lineno) const
{
  try{
    NCMATData::validateElementNameByVersion(s,m_data.version);
  } catch (Error::BadInput&e) {
    NCRYSTAL_THROW2(BadInput,descr()<<": "<<e.what()<<" [in line "<<lineno<<"]");
  }
}

void NC::NCMATParser::handleSectionData_ATOMPOSITIONS(const Parts& parts, unsigned lineno)
{
  if (parts.empty()) {
    if (m_data.atompos.empty())
      NCRYSTAL_THROW2(BadInput,descr()<<": no element positions specified in @ATOMPOSITIONS section (expected in line "<<lineno<<")");
    try {
      m_data.validateAtomPos();
    } catch (Error::BadInput&e) {
      NCRYSTAL_THROW2(BadInput,e.what()<<" (problem in the @ATOMPOSITIONS section ending in line "<<lineno<<")");
    }
    return;
  }
  validateElementName(parts.at(0),lineno);
  if (parts.size()!=4)
    NCRYSTAL_THROW2(BadInput,descr()<<": wrong number of data entries after element name \""<<parts.at(0)<<"\" in line "<<lineno<<" (expected three numbers)");
  std::array<double,3> v;
  for (unsigned i = 0; i<3; ++i) {
    try {
      v[i] = str2dbl_withfractions(parts.at(i+1));
    } catch (Error::BadInput&e) {
      NCRYSTAL_THROW2(BadInput,descr()<<": problem while decoding position parameter #"<<i+1<<" for element \""<<parts.at(0)<<"\" in line "<<lineno<<" : "<<e.what());
    }
  }
  m_data.atompos.emplace_back(parts.at(0),v);
}

void NC::NCMATParser::handleSectionData_SPACEGROUP(const Parts& parts, unsigned lineno)
{
  if (parts.empty()) {
    if (m_data.spacegroup == 0)
      NCRYSTAL_THROW2(BadInput,descr()<<": no spacegroup number specified in @SPACEGROUP section (expected in line "<<lineno<<")");
    try {
      m_data.validateSpaceGroup();
    } catch (Error::BadInput&e) {
      NCRYSTAL_THROW2(BadInput,e.what()<<" (problem in the @SPACEGROUP section ending in line "<<lineno<<")");
    }
    return;
  }
  if ( m_data.spacegroup != 0 || parts.size()>1 )
    NCRYSTAL_THROW2(BadInput,descr()<<": multiple entries specified in @SPACEGROUP section in line "<<lineno<<" (requires just a single number)");
  int sg(0);
  try {
    sg = str2int(parts.at(0));
  } catch (Error::BadInput&e) {
    NCRYSTAL_THROW2(BadInput,descr()<<": problem while decoding spacegroup parameter in line "<<lineno<<" : "<<e.what());
  }
  m_data.spacegroup = sg;
}

void NC::NCMATParser::handleSectionData_DEBYETEMPERATURE(const Parts& parts, unsigned lineno)
{
  if (parts.empty()) {
    if (!m_data.hasDebyeTemperature())
      NCRYSTAL_THROW2(BadInput,descr()<<": missing data in @DEBYETEMPERATURE section (expected in line "<<lineno<<")");
    try {
      m_data.validateDebyeTemperature();
    } catch (Error::BadInput&e) {
      NCRYSTAL_THROW2(BadInput,e.what()<<" (problem in the @DEBYETEMPERATURE section ending in line "<<lineno<<")");
    }
    return;
  }

  if ( m_data.debyetemp_global.has_value() )
    NCRYSTAL_THROW2(BadInput,descr()<<": invalid entries found after global Debye temperature was already specified (offending entries are in line "<<lineno<<")");

  if (parts.size()==1) {
    if ( !m_data.debyetemp_perelement.empty() )//global DT not supported if per-element entries already seen (or if NCMAT version is v4 or later, but that is handled below).
      NCRYSTAL_THROW2(BadInput,descr()<<": invalid entries found in line "<<lineno<<" (missing element name or temperature?)");
    try {
      double tmp = str2dbl(parts.at(0));
      m_data.debyetemp_global = DebyeTemperature{ tmp };
    } catch (Error::BadInput&e) {
      NCRYSTAL_THROW2(BadInput,descr()<<": problem while decoding global Debye temperature in line "<<lineno<<" : "<<e.what());
    }
    if ( m_data.debyetemp_global.has_value() && m_data.version>=4 ) {
      m_data.debyetemp_global.reset();
      NCRYSTAL_THROW2(BadInput,descr()<<": Global Debye temperatures are not allowed in NCMAT v4 or later (problem in line "<<lineno<<")");
    }
  } else if (parts.size()==2) {
    validateElementName(parts.at(0),lineno);
    DebyeTemperature dt;
    try {
      dt = DebyeTemperature{ str2dbl(parts.at(1)) };
    } catch (Error::BadInput&e) {
      NCRYSTAL_THROW2(BadInput,descr()<<": problem while decoding debye temperature for element \""<<parts.at(0)<<"\" in line "<<lineno<<" : "<<e.what());
    }
    m_data.debyetemp_perelement.emplace_back(parts.at(0),dt);
  } else {
    NCRYSTAL_THROW2(BadInput,descr()<<": wrong number of data entries in line "<<lineno);
  }
}

void NC::NCMATParser::handleSectionData_DYNINFO(const Parts& parts, unsigned lineno)
{
  std::string e1 = descr();
  if (parts.empty()) {
    if (!m_active_dyninfo)
      NCRYSTAL_THROW2(BadInput,e1<<": no input found in @DYNINFO section (expected in line "<<lineno<<")");
    try {
      m_active_dyninfo->validate( m_data.version );
    } catch (Error::BadInput&e) {
      NCRYSTAL_THROW2(BadInput,e.what()<<" (problem found in the @DYNINFO section ending in line "<<lineno<<")");
    }

    //Simple validation passed. Clear and return (after squeezing memory):
    for (auto& e : m_active_dyninfo->fields)
      e.second.shrink_to_fit();

    m_dyninfo_active_vector_field = nullptr;
    m_dyninfo_active_vector_field_allownegative = false;
    m_active_dyninfo = nullptr;
    return;
  }
  if (!m_active_dyninfo) {
    m_data.dyninfos.push_back(NCMATData::DynInfo());
    m_active_dyninfo = &m_data.dyninfos.back();
  }

  //all keywords use lowercase characters + '_' and must start with a lower-case
  //letter:

  VectD * parse_target = nullptr;
  NCMATData::DynInfo& di = *m_active_dyninfo;
  Parts::const_iterator itParseToVect(parts.begin()), itParseToVectE(parts.end());
  const std::string& p0 = parts.at(0);

  static_assert('A'<'a'&&'0'<'a'&&'_'<'a',"");
  if ( p0[0] >= 'a' && contains_only(p0,"abcdefghijklmnopqrstuvwxyz_") ) {

    ////////////////////////////
    //line begins with a keyword

    if (parts.size()<2)
      NCRYSTAL_THROW2(BadInput,e1<<": provides no arguments for keyword \""<<p0<<"\" in line "<<lineno);

    m_dyninfo_active_vector_field = nullptr;//new keyword, deactivate active field.
    m_dyninfo_active_vector_field_allownegative = false;//forbid negative numbers except where we explicitly allow them
    ++itParseToVect;//skip keyword if later parsing values into vector
    const std::string& p1 = parts.at(1);

    if ( isOneOf(p0,"fraction","element","type") ) {
      itParseToVect = itParseToVectE;//handle argument parsing here

      /////////////////////////////////////////////////////
      //Handle common fields "fraction", "element", "type":

      if (parts.size()!=2)
        NCRYSTAL_THROW2(BadInput,e1<<": does not provide exactly one argument to keyword \""<<p0<<"\" in line "<<lineno);
      if ( ( p0 == "fraction" && di.fraction != -1.0 )
           || ( p0 == "element" && !di.element_name.empty() )
           || ( p0 == "type" && di.dyninfo_type != NCMATData::DynInfo::Undefined ) )
        NCRYSTAL_THROW2(BadInput,e1<<": keyword \""<<p0<<"\" is specified a second time in line "<<lineno);

      //Specific handling of each:
      if ( p0 == "fraction" ) {
        double fr(-1.0);
        try {
          fr = str2dbl_withfractions(p1);
        } catch (Error::BadInput&e) {
          NCRYSTAL_THROW2(BadInput,e1<<": problem while decoding fraction parameter in line "<<lineno<<" : "<<e.what());
        }
        if ( !(fr<=1.0) || !(fr>0.0) )//this also tests for NaN
          NCRYSTAL_THROW2(BadInput,e1<<": problem while decoding fraction parameter in line "<<lineno<<" (must result in a number greater than 0.0 and at most 1.0)");
        di.fraction = fr;
      } else if ( p0 == "element" ) {
        validateElementName(p1,lineno);
        di.element_name = p1;
      } else if ( p0 == "type" ) {
        if ( p1 == "scatknl" )
          di.dyninfo_type = NCMATData::DynInfo::ScatKnl;
        else if ( p1 == "vdos" )
          di.dyninfo_type = NCMATData::DynInfo::VDOS;
        else if ( p1 == "vdosdebye" )
          di.dyninfo_type = NCMATData::DynInfo::VDOSDebye;
        else if ( p1 == "freegas" )
          di.dyninfo_type = NCMATData::DynInfo::FreeGas;
        else if ( p1 == "sterile" )
          di.dyninfo_type = NCMATData::DynInfo::Sterile;
        else
          NCRYSTAL_THROW2(BadInput,e1<<": invalid @DYNINFO type specified in line "
                          <<lineno<<" (must be one of \"scatknl\", \"vdos\", \"vdosdebye\", \"freegas\", \"sterile\")");
      }
      return;
    }

    //////////////////////////////////////////////////////////////
    //Not a common field, parse into generic DynInfo::fields map :

    if ( di.fields.find(p0) != di.fields.end() )
      NCRYSTAL_THROW2(BadInput,e1<<": keyword \""<<p0<<"\" is specified a second time in line "<<lineno);

    //Setup new vector for parsing into:
    di.fields[p0] = VectD();
    parse_target = &di.fields[p0];
    //Check if supports entry over multiple lines (mostly for keywords
    //potentially needing large number of arguments):
    if ( isOneOf(p0,"sab","sab_scaled","sqw","alphagrid","betagrid","qgrid",
                 "omegagrid","egrid","vdos_egrid", "vdos_density") ) {
      parse_target->reserve(256);//will be squeezed later
      if ( isOneOf(p0,"sqw", "qgrid", "omegagrid") )
        NCRYSTAL_THROW2(BadInput,descr()<<": support for kernels in S(q,w) format and the keyword \""<<p0<<"\" in line "
                        <<lineno<<" is not yet supported (but is planned for inclusion in later NCMAT format versions)");
      m_dyninfo_active_vector_field = parse_target;
      m_dyninfo_active_vector_field_allownegative = (p0=="betagrid"||p0=="omegagrid");
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // Only get here if we have to parse numbers into a vector in DynInfo.fields:

  if ( m_dyninfo_active_vector_field )
    parse_target = m_dyninfo_active_vector_field;
  if ( !parse_target )
    NCRYSTAL_THROW2(BadInput,descr()<<": Unexpected content in line "<<lineno<<": "<<parts.front());
  nc_assert_always( itParseToVect != itParseToVectE );
  std::size_t idx = (itParseToVect-parts.begin());
  std::string tmp_strcache0, tmp_strcache1;
  for (; itParseToVect!=itParseToVectE; ++itParseToVect,++idx) {
    double val;
    const std::string * srcnumstr = &(*itParseToVect);
    const std::string * srcrepeatstr = nullptr;
    //First check for compact notation of repeated entries:
    auto idx_repeat_marker = srcnumstr->find('r');
    if (idx_repeat_marker != std::string::npos) {
      tmp_strcache0.assign(*srcnumstr,0,idx_repeat_marker);
      tmp_strcache1.assign(*srcnumstr,idx_repeat_marker+1,srcnumstr->size()-(idx_repeat_marker+1));
      srcnumstr = &tmp_strcache0;
      srcrepeatstr = &tmp_strcache1;
    }

    unsigned repeat_count = 1;
    try {
      if (srcrepeatstr) {
        int irc = str2int(*srcrepeatstr);
        if (irc<2)
          NCRYSTAL_THROW2(BadInput,"repeated entry count parameter must be >= 2");
        repeat_count = irc;
      }
      val = str2dbl(*srcnumstr);
    } catch (Error::BadInput&e) {
      NCRYSTAL_THROW2(BadInput,e1<<": problem while decoding vector entry #"<<1+(itParseToVect-parts.begin())<<" in line "<<lineno<<" : "<<e.what());
    }
    if (ncisnan(val)||ncisinf(val))
      NCRYSTAL_THROW2(BadInput,e1<<": problem while decoding vector entry #"<<1+(itParseToVect-parts.begin())<<" in line "<<lineno<<" : NaN or infinite number");
    if ( !m_dyninfo_active_vector_field_allownegative && val<0.0 )
      NCRYSTAL_THROW2(BadInput,e1<<": problem while decoding vector entry #"<<1+(itParseToVect-parts.begin())<<" in line "<<lineno<<" : Negative number");
    while (repeat_count--)
      parse_target->push_back(val);
  }
}

void NC::NCMATParser::handleSectionData_DENSITY(const Parts& parts, unsigned lineno)
{
  if (parts.empty()) {
    if (!m_data.density)
      NCRYSTAL_THROW2(BadInput,descr()<<": no input found in @DENSITY section (expected in line "<<lineno<<")");
    try {
      m_data.validateDensity();
    } catch (Error::BadInput&e) {
      NCRYSTAL_THROW2(BadInput,e.what()<<" (problem in the @DENSITY section ending in line "<<lineno<<")");
    }
    return;
  }
  if ( m_data.density>0.0 )
    NCRYSTAL_THROW2(BadInput,descr()<<": too many lines in @DENSITY section in line "<<lineno);
  if (parts.size()!=2)
    NCRYSTAL_THROW2(BadInput,descr()<<": wrong number of entries on line "<<lineno<<" in @DENSITY section");
  double density_val;
  try {
    density_val = str2dbl(parts.at(0));
  } catch (Error::BadInput&e) {
    NCRYSTAL_THROW2(BadInput,descr()<<": problem while decoding density value in line "<<lineno<<" : "<<e.what());
  }
  if (parts.at(1)=="atoms_per_aa3") {
    m_data.density_unit = NCMATData::ATOMS_PER_AA3;
    m_data.density = density_val;
  } else if (parts.at(1)=="kg_per_m3") {
    m_data.density_unit = NCMATData::KG_PER_M3;
    m_data.density = density_val;
  } else if (parts.at(1)=="g_per_cm3") {
    m_data.density_unit = NCMATData::KG_PER_M3;
    m_data.density = density_val * 1000.0;
  } else {
    NCRYSTAL_THROW2(BadInput,descr()<<": invalid density unit in line "<<lineno);
  }
  if ( !(m_data.density>0.0) )
    NCRYSTAL_THROW2(BadInput,descr()<<": invalid density value in line "<<lineno);
}

void NC::NCMATParser::handleSectionData_STATEOFMATTER(const Parts& parts, unsigned lineno)
{
  if (parts.empty()) {
    if (!m_data.stateOfMatter.has_value())
      NCRYSTAL_THROW2(BadInput,descr()<<": no input found in @STATEOFMATTER section (expected in line "<<lineno<<")");
    return;
  }
  if ( m_data.stateOfMatter.has_value() )
    NCRYSTAL_THROW2(BadInput,descr()<<": too many lines in @STATEOFMATTER section in line "<<lineno);
  if (parts.size()!=1)
    NCRYSTAL_THROW2(BadInput,descr()<<": wrong number of entries on line "<<lineno<<" in @STATEOFMATTER section");
  if (parts.at(0)=="solid") {
    m_data.stateOfMatter = NCMATData::StateOfMatter::Solid;
  } else if (parts.at(0)=="liquid") {
    m_data.stateOfMatter = NCMATData::StateOfMatter::Liquid;
  } else if (parts.at(0)=="gas") {
    m_data.stateOfMatter = NCMATData::StateOfMatter::Gas;
  } else {
    NCRYSTAL_THROW2(BadInput,descr()<<": invalid state of matter type specified in @STATEOFMATTER section in line "
                    <<lineno<<" (must be \"solid\", \"liquid\", or \"gas\")");
  }
}

void NC::NCMATParser::handleSectionData_OTHERPHASES(const Parts& parts, unsigned lineno)
{

  if (parts.empty()) {
    if (m_data.otherPhases.empty())
      NCRYSTAL_THROW2(BadInput,descr()<<": no input found in @OTHERPHASES section (expected in line "<<lineno<<")");
    return;
  }
  if (parts.size()<2)
    NCRYSTAL_THROW2(BadInput,descr()<<": wrong number of entries on line "<<lineno<<" in @OTHERPHASES section");
  auto volfrac = StrView(parts.at(0)).toDbl();
  if ( !volfrac.has_value() || !(volfrac.value()>0.0) || !(volfrac.value()<1.0) )
    NCRYSTAL_THROW2(BadInput,descr()<<": invalid volume fraction \""<<parts.at(0)<<"\" specified in @OTHERPHASES section in line "
                    <<lineno<<" (must be a floating point number greater than 0.0 and less than 1.0)");
  std::string cfgstr = parts.at(1);
  for ( auto i : ncrange(2,(int)parts.size()) ) {
    cfgstr += ' ';//normalise whitespace to single space (as documented in the NCMAT doc)
    cfgstr += parts.at(i);
  }

  m_data.otherPhases.emplace_back(volfrac.value(),cfgstr);
}

void NC::NCMATParser::handleSectionData_TEMPERATURE(const Parts& parts, unsigned lineno)
{
  if (parts.empty()) {
    if ( !m_data.temperature.has_value())
      NCRYSTAL_THROW2(BadInput,descr()<<": no input found in @TEMPERATURE section (expected in line "<<lineno<<")");
    try {
      m_data.validateTemperature();
    } catch (Error::BadInput&e) {
      NCRYSTAL_THROW2(BadInput,e.what()<<" (problem in the @TEMPERATURE section ending in line "<<lineno<<")");
    }
    return;
  }
  if ( m_data.temperature.has_value() )
    NCRYSTAL_THROW2(BadInput,descr()<<": too many lines in @TEMPERATURE section in line "<<lineno);
  if ( !isOneOf((int)parts.size(),1,2) )
    NCRYSTAL_THROW2(BadInput,descr()<<": wrong number of entries on line "<<lineno<<" in @TEMPERATURE section");
  auto temperature_value = StrView(parts.back()).toDbl();
  if ( !temperature_value.has_value() )
    NCRYSTAL_THROW2(BadInput,descr()<<": problem decoding temperature value in line "<<lineno);
  if ( !(temperature_value.value()>0.0) || !(temperature_value.value()<=1e6) )//NB: use same thresholds in NCNCMATData.cc
    NCRYSTAL_THROW2(BadInput,descr()<<": out of range temperature value in line "<<lineno);
  if ( parts.size() == 2 && parts.front() != "default" )
    NCRYSTAL_THROW2(BadInput,descr()<<": Entry in line "<<lineno
                    <<" must be a temperature value or the keyword \"default\" followed by a temperature value");
  NCMATData::TemperatureType temptype = ( parts.size() == 1
                                          ? NCMATData::TemperatureType::Fixed
                                          : NCMATData::TemperatureType::Default );
  m_data.temperature.emplace( Temperature{ temperature_value.value() }, temptype );

}


void NC::NCMATParser::handleSectionData_ATOMDB(const Parts& parts, unsigned lineno)
{
  if (parts.empty())
    return;//end of section, nothing to do
  if ( parts.at(0)!="nodefaults" )
    validateElementName(parts.at(0),lineno);
  m_data.atomDBLines.emplace_back(parts);
}

void NC::NCMATParser::handleSectionData_CUSTOM(const Parts& parts, unsigned)
{
  if (parts.empty())
    return;//end of section, nothing to do
  nc_assert(!m_data.customSections.empty());
  m_data.customSections.back().second.push_back(parts);
}
