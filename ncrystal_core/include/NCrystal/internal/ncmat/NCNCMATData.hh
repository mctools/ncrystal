#ifndef NCrystal_NCMATData_hh
#define NCrystal_NCMATData_hh

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

#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  class NCRYSTAL_API NCMATData : private MoveOnly {

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
    //Validate. NB: If the aim is to support T and D as aliases for H3 and H2
    //(it usually is), you should call unaliasElementNames() before calling
    //validate():
    void validate() const;

    //Metadata
    int version = 0;
    constexpr static int latest_version = static_cast<int>(supported_ncmat_format_version_max);
    DataSourceName sourceDescription;

    //convenience (for a validated instance, this is the same as hasCell or hasAtomPos):
    bool hasUnitCell() const;

    //State of matter.
    enum class StateOfMatter { Solid, Gas, Liquid };
    Optional<StateOfMatter> stateOfMatter;
    static const char * stateOfMatter2Str(StateOfMatter);

    //@CELL
    struct Cell {
      std::array<double,3> lengths = std::array<double,3>{0.,0.,0.};
      std::array<double,3> angles = std::array<double,3>{0.,0.,0.};//in degrees
    } cell;
    bool hasCell() const;
    void validateCell() const;

    //@ATOMPOSITIONS
    using Pos = std::array<double,3>;
    std::vector<std::pair<std::string,Pos> > atompos;
    bool hasAtomPos() const { return !atompos.empty(); }
    void validateAtomPos() const;

    //@TEMPERATURE (v7+)
    enum class TemperatureType { Default, Fixed };
    Optional<std::pair<Temperature,TemperatureType>> temperature;
    void validateTemperature() const;

    //@SPACEGROUP
    int spacegroup = 0;
    bool hasSpaceGroup() const { return spacegroup>0; }
    void validateSpaceGroup() const;

    //@DEBYETEMPERATURE (note that if absent, VDOSDebye DynInfo sections can have their own debye_temp field)
    Optional<DebyeTemperature> debyetemp_global;
    std::vector<std::pair<std::string,DebyeTemperature> > debyetemp_perelement;
    bool hasDebyeTemperature() const { return debyetemp_global.has_value() || !debyetemp_perelement.empty(); }
    void validateDebyeTemperature() const;

    //@DYNINFO
    struct NCRYSTAL_API DynInfo {
      enum DynInfoType { Sterile, FreeGas, VDOSDebye, VDOS, ScatKnl, Undefined };
      static const char * diType2Str(DynInfoType);
      const char * typeStr() const { return diType2Str(dyninfo_type); }
      DynInfoType dyninfo_type = Undefined;
      std::string element_name;
      double fraction = -1.0;
      //For simplicity, keep all numerical entries (except the always present
      //"fraction") as key->vector<double>, even though some keys only
      //correspond to single numbers. The validate() method will take care of
      //checking correct number and presence of all values, and is also here
      //kept synchronised with the NCMAT loader code.
      typedef std::map<std::string,VectD> FieldMapT;
      FieldMapT fields;//Stuff like "temperature", "alphagrid", "betagrid", ...
      void validate( int version ) const;//throws BadInput in case of problems specific to this DynInfo object (missing/wrong fields, etc.)
    };
    std::vector<DynInfo> dyninfos;
    bool hasDynInfo() const { return !dyninfos.empty(); }

    //@DENSITY
    enum DensityUnit { ATOMS_PER_AA3, KG_PER_M3 } density_unit = ATOMS_PER_AA3;
    double density = 0.0;
    bool hasDensity() const { return density!=0.0; }
    void validateDensity() const;

    //@ATOMDB
    typedef VectS AtomDBLine;
    typedef std::vector<AtomDBLine> AtomDBLines;
    AtomDBLines atomDBLines;
    bool hasAtomDB() const { return !atomDBLines.empty(); }
    void validateAtomDB() const;

    //@OTHERPHASES
    std::vector<std::pair<double,std::string>> otherPhases;//{volfrac,cfgstr}
    bool hasOtherPhases() const { return !otherPhases.empty(); }
    void validateOtherPhases() const;

    //@CUSTOM_xxx section contents. Kept in order of appearance in file.
    typedef VectS CustomLine;
    typedef std::vector<CustomLine> CustomSectionData;
    typedef std::string CustomSectionName;
    std::vector<std::pair<CustomSectionName,CustomSectionData>> customSections;

    //Helper for validating "element" names (throws BadInput in case of problems):
    static void validateElementNameByVersion(const std::string&, unsigned version );

    //Call this once after filling fields (and before calling validate()!) in
    //order to unalias all element names (D->H2, T->H3).
    void unaliasElementNames();

    //Plumbing:
    NCMATData() = default;
    ~NCMATData() = default;
    NCMATData( const NCMATData& ) = delete;
    NCMATData& operator=( const NCMATData& ) = delete;
    NCMATData( NCMATData&& ) = default;
    NCMATData& operator=( NCMATData&& ) = default;

    // Experimental and INCOMPLETE support for dumping the raw NCMAT data into JSON format:
    void toJSON( std::ostream& ) const;

  private:
    void validateElementName( const std::string& ) const;
  };

}

#endif
