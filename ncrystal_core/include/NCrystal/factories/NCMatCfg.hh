#ifndef NCrystal_MatCfg_hh
#define NCrystal_MatCfg_hh

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
#include "NCrystal/text/NCTextData.hh"

namespace NCRYSTAL_NAMESPACE {

  class SCOrientation;
  class MatInfoCfg;

  class NCRYSTAL_API MatCfg {
  public:

    //////////////////////////////////////////////////////////////////////////////
    //                                                                          //
    // Class which is used to encapsulate the high-level configuration of a     //
    // given material.                                                          //
    //                                                                          //
    // The configuration can either be set up and edited via the API below, or  //
    // simply by providing a so-called configuration string in the              //
    // constructor. Once a MatCfg object is set up, it is possible to use it    //
    // to load corresponding physics objects, by passing it to the global       //
    // createInfo, createScatter, or createAbsorption functions (cf. NCFact.hh).//
    //                                                                          //
    // In the simplest case, the configuration will combine one source of input //
    // (text) data (for instance an NCMAT file) with a set of configuration     //
    // parameters (for instance the material temperature or a single crystal    //
    // orientation). In more advanced cases, it is also possible for a MatCfg   //
    // object to define a multi-phase material. Such a MatCfg object does not   //
    // contain input data or (with a few exceptions) configuration parameters,  //
    // but instead simply provides a list of child MatCfg objects, one for each //
    // phase.                                                                   //
    //                                                                          //
    // NB: Even if the MatCfg object is not "multi-phase", and only refers to   //
    // a single source of input data, the physics objects loaded with it might  //
    // still represent a multi-phase material. This follows from the fact that  //
    // for instance a single NCMAT file can provide a multi-phase material, but //
    // the MatCfg object itself can not know whether or not this is the case.   //
    //                                                                          //
    //////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    // Cfg-string constructor:                                                 //
    /////////////////////////////////////////////////////////////////////////////
    //

    MatCfg( const std::string& datafile_and_parameters );
    MatCfg( const char* datafile_and_parameters );

    // Construct material configuration by providing a configuration
    // string. Such a cfg-string must at least contain the name identifying a
    // piece of input data (in any supported format). This can for instance be a
    // filename, and more generally it must be anything which can be used look
    // up a TextData object via the infrastructure in NCFactImpl.hh, after
    // constructing a TextDataPath object with it (see also NCDataSources.hh for
    // how to modify this search).
    //
    // Additionally, it is optionally possible to set configuration parameters
    // by appending one or more occurances of "<name>=<value>" to the file name,
    // all separated with semi-colons (;). Alternatively, the parameters can
    // also be set via specific set_xxx() methods below.
    //
    // Thus the single line:
    //
    //    MatCfg cfg("myfile.ncmat;temp=77.0;density=2.3gcm3");
    //
    // Has the same effect as the three lines:
    //
    //    MatCfg cfg("myfile.ncmat");
    //    cfg.set_temp( Temperature{77.0} );
    //    cfg.set_density( Density{2.3} );
    //
    // See notes below concerning multiphase materials, in-file embedding of
    // parameters, data caching, and instantiation directly from text data.
    //
    // The supported parameter names are (for the latest NCrystal release)
    // documented at:
    //
    //   https://github.com/mctools/ncrystal/wiki/CfgRefDoc
    //
    // In case your version of NCrystal differs from the one mentioned on the
    // wiki page above, you can instead generate the same documentation by
    // invoking one of the following:
    //
    //   nctool --doc                              #From the command line
    //   NCrystal::MatCfg::genDoc(std::cout);    //#From C++ code.
    //   print( NCrystal.generateCfgStrDoc() )     #From Python code
    //
    // C++ savvy users might also wish to take a look at the internal header
    // file, NCCfgVars.hh, where the various parameters are actually defined.
    //
    // Multiphase materials can be composed directly in cfg-strings by simple
    // additive composition from existing materials, as will be explained in the
    // following. Note that very complicated multiphase materials, or those with
    // nanostructures resulting in SANS physics, might be better described in
    // NCMAT data. In that case, resulting physics objects might be multiphased,
    // although the MatCfg object will seem to be single phase.
    //
    // The configuration string syntax used to define a multiphase material with
    // N phases is:
    //
    //   "phases<FRAC1*CFG1&..&FRACN*CFGN>[;COMMONCFG]"
    //
    // Here FRACi is the (by-volume) fraction of the i'th component, and
    // COMMONCFG, CFG1, ..., CFGN are cfg strings. The final cfg string of the
    // i'th component is "CFG1;COMMONCFG", so COMMONCFG will be the most
    // convenient place to change parameters which should affect all phases,
    // such as material temperature.  One special exception is the use of
    // "phasechoice" and "density" parameters (see below) in the COMMONCFG
    // string, which will always be applied only to the top-level MatCfg object,
    // and not be passed on to the child phases (doing so would defeat the
    // purpose of using phasechoice parameters).
    //
    // Examples:
    //
    // * Equal (by volume) mix of Al and Mg powder:
    //   "phases<0.5*Al_sg225.ncmat&0.5*Mg_sg194.ncmat>"
    // * Same at T=20K:
    //   "phases<0.5*Al_sg225.ncmat&0.5*Mg_sg194.ncmat>;temp=20K"
    // * Aluminium with voids:
    //   "phases < 0.99 * Al_sg225.ncmat & 0.01 * void.ncmat >"
    // * Alumina powder suspended in water:
    //   "phases<0.01*Al2O3_sg167_Corundum.ncmat&0.99*LiquidWaterH2O_T293.6K.ncmat>"
    //
    // Also note that in case where a material is defined using the "phases<..>"
    // multiphase syntax, but only one phase is actually specified - the
    // resulting MatCfg object will end up as a regular single phase object.
    //
    // A multi-phase MatCfg object is primarily useful since it can be passed to
    // createInfo, createAbsorption, or createScatter as any MatCfg object and
    // receive Info, Absorption, or Scatter physics objects in return. But note
    // that some methods on multi-phase MatCfg objects are unavailable
    // (i.e. will throw exceptions). This can for instance happen when asking
    // for the value of a parameter which is not the same in all phases. Setting
    // parameters is on the other hand always possible, and merely leads to the
    // value being set for all phases (with the exception of phasechoice/density
    // parameters as noted above). Additionally, a multi-phase MatCfg object
    // never has text data directly available, since such data is associated
    // with the individual phases.
    //
    // In addition to the cfg-string syntax above, it is also possible to setup
    // multi-phase configurations with the following constructors:

    using Phase = std::pair<double,MatCfg>;
    using PhaseList = std::vector<Phase>;
    MatCfg( PhaseList&& );
    MatCfg( const PhaseList& );

    ////////////////////////////////////////////////////////////////////////////
    // For additional flexiblity, it is also possible to instantiate MatCfg
    // objects by providing the bulk text data directly instead of a filename -
    // in the form of either a TextData object or a string with the content. In
    // either case, there should be enough information available that the
    // framework can figure out the type of data provided. In order of
    // preference, this will be done by looking at #1 the explicit dataType
    // string provided as a parameter to the constructor (if non-empty), #2 the
    // dataType from the TextData object, #3 if the TextData content starts with
    // the 5 characters "NCMAT", the type will be inferred to be "ncmat", and
    // finally as a last resort, #4 the file extension of the filename provided,
    // if any. If all these options fail to determine the data format, a
    // BadInput exception will be thrown.
    //
    // Note that providing data manually as a string will result in a new
    // TextData object being constructed internally with a new dataUID(), with
    // implications for downstream-caching.

    explicit MatCfg( TextDataSP, std::string parameters = {} );

    //To avoid overload ambiguities, the constructor accepting raw data directly
    //is encapsulated in a static creation method:
    static MatCfg createFromRawData( std::string&& data,
                                     std::string parameters = {},
                                     std::string dataType = {} );

    //Generate cfg-string documentation, either full output including
    //explanations, short output with one variable per line, or as a JSON data
    //structure:
    enum class GenDocMode { TXT_SHORT, TXT_FULL, JSON };
    static void genDoc( std::ostream&, GenDocMode = GenDocMode::TXT_FULL );

    /////////////////////////////////////////////////////////////////////////////
    // Methods for setting parameters:                                         //
    /////////////////////////////////////////////////////////////////////////////
    //
    // Directly set from C++ code:
    void set_temp( Temperature );
    void set_temp( double );
    void set_dcutoff( double );
    void set_dcutoffup( double );
    void set_mos( MosaicityFWHM );
    void set_mosprec( double );
    void set_sccutoff( double );
    void set_dirtol( double );
    void set_coh_elas( bool );
    void set_incoh_elas( bool );
    void set_sans( bool );
    void set_inelas( const std::string& );
    void set_infofactory( const std::string& );
    void set_scatfactory( const std::string& );
    void set_absnfactory( const std::string& );
    void set_lcmode( std::int_least32_t );
    void set_ucnmode( const Optional<UCNMode>& );
    void set_vdoslux( int );
    void set_atomdb( const std::string& );
    void set_lcaxis( const LCAxis& );
    void set_dir1( const HKLPoint&, const LabAxis& );
    void set_dir1( const CrystalAxis&, const LabAxis& );
    void set_dir2( const HKLPoint&, const LabAxis& );
    void set_dir2( const CrystalAxis&, const LabAxis& );
    void set_dir1( const OrientDir& );
    void set_dir2( const OrientDir& );
    void setOrientation( const SCOrientation& );//sets dir1+dir2+dirtol

    ///////////////////////////////////////////////////////////////////////////////////
    //As explained above, MatCfg objects optionally retain a "density state",
    //which can be either a scaling of, or an absolute override of, the density
    //that will otherwise be created from the configuration. Note that final
    //density value must be obtained from the loaded Info objects, not the
    //MatCfg objects which can't know the underlying material density which it
    //is possibly scaling.

    void set_density( const Density& );
    void set_density( const NumberDensity& );
    void scale_density( double );//scales existing state (NB: repeated calls have commulative effect!)
    void set_density( const DensityState& );// (NB: repeated calls have commulative effect if argument is scale factor!)
    DensityState get_density() const;
    bool hasDensityOverride() const;//True if get_density() does not yield {SCALEFACTOR,1.0}.

    //
    // Set parameters from a string, using the same format as that supported by
    // the constructor, e.g. "par1=val1;...;parn=valn":
    void applyStrCfg( const std::string&  );

    /////////////////////////////////////////////////////////////////////////////
    // Methods for accessing parameters or derived information:                //
    /////////////////////////////////////////////////////////////////////////////
    //
    // Directly access from C++ (might error for multiphase cfgs in case of
    // inconsistent values across phases):

    Temperature get_temp() const;
    double get_dcutoff() const;
    double get_dcutoffup() const;
    MosaicityFWHM get_mos() const;
    double get_mosprec() const;
    double get_sccutoff() const;
    double get_dirtol() const;
    const LCAxis& get_lcaxis() const;
    std::string get_infofactory() const;
    std::string get_scatfactory() const;
    std::string get_absnfactory() const;
    std::int_least32_t get_lcmode() const;
    StrView get_ucnmode_str() const;
    Optional<UCNMode> get_ucnmode() const;
    int get_vdoslux() const;
    std::string get_atomdb() const;
    std::vector<VectS> get_atomdb_parsed() const;
    bool get_coh_elas() const;
    bool get_incoh_elas() const;
    bool get_sans() const;
    std::string get_inelas() const;
    OrientDir get_dir1() const;
    OrientDir get_dir2() const;

    // Check if single/layered-single crystal, or access dir1+dir2+dirtol as
    // SCOrientation (do not call any of these for multiphase cfgs):
    bool isSingleCrystal() const;//true if mos or orientation parameters are set
    bool isLayeredCrystal() const;//true if lcaxis parameter is set
    SCOrientation createSCOrientation() const;

    //Serialise in various forms. Note that if the MatCfg object was constructed
    //from anonymous text data and not named data (filename), the data returned
    //can not necessarily be used for an exact recreation of an equivalent
    //MatCfg object.
    //Note that for multiphase MatCfg objects, it is not allowed to call
    //toEmbeddableCfg() or toStrCfg(include_datafile=false).
    std::string toStrCfg( bool include_datafile = true ) const;
    std::string toEmbeddableCfg() const;//Produces a string like "NCRYSTALMATCFG[temp=500.0;dcutoff=0.2]"
                                        //which can be embedded in data files.
    std::string toJSONCfg() const;//same info as toStrCfg but as JSON struct.
    void dump( std::ostream &out, bool add_endl=true ) const;

    /////////////////////////////////////////////////////////////////////////////
    // Associated text data (content and meta data of data file).              //
    /////////////////////////////////////////////////////////////////////////////
    //
    // For non-multiphase cfg objects, access associated text and meta data:
    const TextData& textData() const;
    TextDataSP textDataSP() const;
    const TextDataUID textDataUID() const;

    //For non-multiphase cfg objects, access data type of text data (usually the
    //extension of the data file). Calling this method should be the only manner
    //of determining the format of the associated textData.
    const std::string& getDataType() const;

    //For non-multiphase cfg objects, access data source name (for output only,
    //code should not parse this!):
    const DataSourceName& getDataSourceName() const;

    /////////////////////////////////////////////////////////////////////////////
    // Methods for multi-phase configurations:                                 //
    /////////////////////////////////////////////////////////////////////////////
    //
    bool isSinglePhase() const;
    bool isMultiPhase() const;
    const PhaseList& phases() const;//empty list if isSinglePhase.


    /////////////////////////////////////////////////////////////////////////////
    // Phasechoices                                                            //
    /////////////////////////////////////////////////////////////////////////////

    //Pick out child phases by index in Info::getPhases() list. Adding phase
    //choices here (or with the "phasechoice" pseudo-parameter) can be used to
    //ensure that passing the MatCfg object to high-level createXXX functions
    //will result in the physics object corresponding to the given phase will be
    //returned rather than the top-level objects (NOTE: this parameter picks out
    //children at the Info object level, not at the MatCfg-level child
    //phases). As child phase Info objects might themselves be multiphase, the
    //method can be called more than once to navigate further into the Info
    //object tree of child phases:
    using PhaseChoices = SmallVector_IC<unsigned,4,SVMode::LOWFOOTPRINT>;
    const PhaseChoices& getPhaseChoices() const;
    void appendPhaseChoices( const PhaseChoices& );
    void appendPhaseChoice( unsigned );

    /////////////////////////////////////////////////////////////////////////////
    // Copy/assign/move/clone/destruct                                         //
    /////////////////////////////////////////////////////////////////////////////
    //
    // Copy/assignment/cloning is allowed and is a priori cheap, since internal
    // data structures are shared until modified (aka copy-on-write):
    //
    MatCfg(const MatCfg&);
    MatCfg& operator=(const MatCfg&);
    MatCfg( MatCfg&& );
    MatCfg& operator=(MatCfg&&);
    ~MatCfg();
    MatCfg clone() const;

    /////////////////////////////////////////////////////////////////////////////
    // Interface for NCrystal factory infrastructure:                          //
    /////////////////////////////////////////////////////////////////////////////

    // Advanced interface, which is intended solely for use by the NCrystal
    // factory infrastructure, for validation, factory selection, and to avoid
    // unnecessary re-initialisation of expensive objects.

    //Create a "thinned" clone which does not keep a strong reference to the
    //associated TextData object, but which still retains all other information
    //(including the textDataUID) needed for comparisons, dumps, etc. Thus, such
    //a clone is suitable for e.g. cache keys which must be kept around for
    //comparison but which should not retain strong references to potentially
    //large text data.
    MatCfg cloneThinned() const;
    bool isThinned() const;

    //Clone with phase choices or density discarded:
    MatCfg cloneWithoutPhaseChoices() const;
    MatCfg cloneWithoutDensityState() const;

    //Trivial means not multiphase, no phase choices, no density state override:
    bool isTrivial() const;

    //Verify that the parameter values are not inconsistent or
    //incomplete. Throws BadInput exception if they are.
    void checkConsistency() const;

    //Unspecified ordering (for usage as map keys):
    bool operator<(const MatCfg&) const;

    //Low-level access to parameter data:
    const Cfg::CfgData& rawCfgData() const;//NB: singlephase only fct
    void apply( const Cfg::CfgData& );//NB: applies to all phases if multiphase

    //////////////////////////////////////////////////////////////////////
    // Obsolete functions kept temporarily for backwards compatibility. //
    //////////////////////////////////////////////////////////////////////

    //The get_packfact function was replaced with get_density which plays a
    //slightly different role:
    double get_packfact() const { return 1.0; }
    //The ignorefilecfg keyword is no longer supported (it was too complicated
    //to support):
    bool ignoredEmbeddedConfig() const { return false; }

  private:
    struct constructor_args;
    explicit MatCfg( constructor_args&& );
    struct from_raw_t {};
    explicit MatCfg( from_raw_t, std::string&& data, std::string pars, std::string ext );
    class Impl;
    class Impl2;
    COWPimpl<Impl> m_impl;
    COWPimpl<Impl2> m_impl2;
    OptionalTextDataSP m_textDataSP;
    friend class MatInfoCfg;
  };

  NCRYSTAL_API std::ostream& operator<< (std::ostream&, const MatCfg&);

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  inline void MatCfg::set_temp( double t ) { set_temp(Temperature{t}); }
  inline MatCfg MatCfg::clone() const { return MatCfg(*this); }
  inline void MatCfg::scale_density( double val )
  {
    if ( val!= 1.0 ) {
      DensityState ds;
      ds.type = DensityState::Type::SCALEFACTOR;
      ds.value = val;
      set_density(ds);
    }
  }
  inline std::ostream& operator<< (std::ostream& os, const MatCfg& cfg) { cfg.dump(os,false); return os; }

  static_assert(std::is_move_constructible<MatCfg>::value, "");
  static_assert(std::is_move_assignable<MatCfg>::value, "");
  static_assert(std::is_copy_constructible<MatCfg>::value, "");
  static_assert(std::is_copy_assignable<MatCfg>::value, "");
}

#endif
