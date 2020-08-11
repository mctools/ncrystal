#ifndef NCrystal_NCMATData_hh
#define NCrystal_NCMATData_hh

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

  class NCMATData {

    //Data structure which holds the parsed information equivalent to what is
    //specified in an .ncmat file (typically it will be filled out by using the
    //NCMATParser class to parse such file). A simple basic validation of
    //sections can be performed by calling the validate method, which will throw
    //a BadInput exception in case of missing mandatory sections, partially
    //filled sections, or obviously wrong input (e.g. numbers out of range,
    //etc.). The validation here is simple and can of course not catch all
    //issues. Some might only be found when actually processing the data
    //further, transforming it to NCInfo object, etc. There are also
    //section-specific validation method, which are used both internally in the
    //validate() method, as well as inside the NCMATParser class, to validate
    //each section as it is read from file (resulting in more specific error
    //messages - nb. the DynInfo::validate method does not include
    //sourceFullDescr in error messages, all others do).

  public:
    NCMATData() { clear(); }
    ~NCMATData();

    void validate() const;
    void clear();
    void swap(NCMATData& other);

    //Metadata
    int version;
    std::string sourceDescription;
    std::string sourceType;
    std::string sourceFullDescr;

    //convenience (for a validated instance, this is the same as hasCell or hasAtomPos):
    bool hasUnitCell() const;

    //@CELL
    struct Cell {
      std::array<double,3> lengths;
      std::array<double,3> angles;//in degrees
    } cell;
    bool hasCell() const;
    void validateCell() const;

    //@ATOMPOSITIONS
    using Pos = std::array<double,3>;
    std::vector<std::pair<std::string,Pos> > atompos;
    bool hasAtomPos() const { return !atompos.empty(); }
    void validateAtomPos() const;

    //@SPACEGROUP
    int spacegroup;
    bool hasSpaceGroup() const { return spacegroup>0; }
    void validateSpaceGroup() const;

    //@DEBYETEMPERATURE
    double debyetemp_global;
    std::vector<std::pair<std::string,double> > debyetemp_perelement;
    bool hasDebyeTemperature() const { return debyetemp_global || !debyetemp_perelement.empty(); }
    void fillPerElementDebyeTempMap(std::map<std::string,double>&) const;
    void validateDebyeTemperature() const;

    //@DYNINFO
    struct DynInfo {
      DynInfo();
      ~DynInfo();
      enum DynInfoType { Sterile, FreeGas, VDOSDebye, VDOS, ScatKnl, Undefined };
      DynInfoType dyninfo_type;
      std::string element_name;
      double fraction;
      //For simplicity, keep all numerical entries (except the always present
      //"fraction") as key->vector<double>, even though some keys only
      //correspond to single numbers. The validate() method will take care of
      //checking correct number and presence of all values, and is also here
      //kept synchronised with the NCMAT loader code.
      typedef std::map<std::string,VectD > FieldMapT;
      FieldMapT fields;//Stuff like "temperature", "alphagrid", "betagrid", ...
      void validate() const;//throws BadInput in case of problems specific to this DynInfo object (missing/wrong fields, etc.)
    };
    std::vector<DynInfo> dyninfos;
    bool hasDynInfo() const { return !dyninfos.empty(); }

    //@DENSITY
    double density;
    enum DensityUnit { ATOMS_PER_AA3, KG_PER_M3 } density_unit;
    bool hasDensity() const { return density!=0; }
    void validateDensity() const;

    //Misc helpers;
    static bool couldBeElementName(const std::string&);
  private:
    void validateElementName( const std::string& s ) const;
  };

}

#endif
