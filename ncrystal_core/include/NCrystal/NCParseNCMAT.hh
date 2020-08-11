#ifndef NCrystal_ParseNCMAT_hh
#define NCrystal_ParseNCMAT_hh

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

#include "NCrystal/NCNCMATData.hh"
#include "NCrystal/NCFile.hh"

namespace NCrystal {

  class NCMATParser {
  public:

    //Parse .ncmat files.

    //Parse input. Will throw BadInput exceptions in case of problems. It will
    //do some rudimentary syntax checking (including presence/absence of data
    //sections and will call NCMATData::validate), but not a full validation of
    //data (a more complete validation is typically carried out afterwards by
    //the NCMAT Loader code). It will always clear the input pointer
    //(i.e. release/close the resource).
    NCMATParser(UniquePtr<TextInputStream>& input);
    ~NCMATParser();

    //Access parsed data (this transfers the contained data to the passed
    //object, so can only be done once):
    void getData(NCMATData&);

  private:

    typedef std::vector<std::string> Parts;
    void parseFile( TextInputStream* );
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

    //Collected data:
    NCMATData m_data;

    //Long vectors of data in @DYNINFO sections are kept in:
    NCMATData::DynInfo * m_active_dyninfo;
    VectD * m_dyninfo_active_vector_field;
    bool m_dyninfo_active_vector_field_allownegative;

  };
}

#endif
