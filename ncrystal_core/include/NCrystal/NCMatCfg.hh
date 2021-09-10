#ifndef NCrystal_MatCfg_hh
#define NCrystal_MatCfg_hh

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

#include "NCrystal/NCTypes.hh"
#include "NCrystal/NCVariant.hh"
#include "NCrystal/NCTextData.hh"
#include <ostream>

namespace NCrystal {

  class SCOrientation;

  class MatInfoCfg;

  class NCRYSTAL_API MatCfg {
  public:

    /////////////////////////////////////////////////////////////////////////////
    // Class which is used to define the configuration of a given material,    //
    // which usually includes not only the path to a material data file in any //
    // supported format, but also other parameters needed for a particular     //
    // use-case, can be set up such as temperature, packing factor or          //
    // orientation.                                                            //
    /////////////////////////////////////////////////////////////////////////////


    /////////////////////////////////////////////////////////////////////////////
    // Constructor:                                                            //
    /////////////////////////////////////////////////////////////////////////////
    //
    // Construct material configuration by supplying at least the path to a data
    // file (in any supported format). Optionally, it is possible to set
    // configuration parameters by appending one or more occurances of
    // "<name>=<value>" to the file name, all separated with semi-colons (;).
    //
    // Thus the single line:
    //
    //    MatCfg cfg("myfile.ncmat;temp=77.0;packfact=0.8");
    //
    // Has the same effect as the three lines:
    //
    //    MatCfg cfg("myfile.ncmat");
    //    cfg.set_temp(77.0);
    //    cfg.set_packfact(0.8);
    //
    // The supported variable names and the set methods which they invoke
    // are documented below.
    //

    MatCfg( const std::string& datafile_and_parameters );
    MatCfg( const char* datafile_and_parameters );
    //
    // Note that it is also possible to set the (initial) values of parameters,
    // by embedding a statement like NCRYSTALMATCFG[temp=500.0;dcutoff=0.2] into
    // the data itself (in a suitable location in the file of course, i.e. in a
    // comment). These parameters can of course still be overridden, and such
    // embedded configuration can be ignored entirely by appending
    // ";ignorefilecfg" to the filename:
    //
    //    MatCfg cfg("myfile.ncmat;ignorefilecfg").
    //
    // Also note that for convenience some parameters optionally support the
    // specification of units when setting their values via strings. For
    // example, the following ways of setting the temperature all correspond to
    // the freezing point of water:
    //
    //     ";temp=273.15" ";temp=273.15K" ";temp=0C" "temp=32F"
    //
    // No matter how the temperature was set, get_temp() will return the value
    // in kelvin (273.15) afterwards.
    //
    // Important notice: The MatCfg constructor immediately opens the specified
    // input file and caches the content in memory (without needless copies in
    // case of in-memory files), and makes it available via the textData and
    // textDataSP methods. Thus, keeping a very large number of MatCfg objects
    // around can potentially incur a sizeable overhead. However, the usage
    // should be rather extreme before this can become a problem (e.g. keeping
    // 1000 MatCfg objects around which all refer to unusually large 1MB files
    // could potentially mean 1GB of memory retained). As most input files are
    // typically 100 times smaller than this and since users typically only keep
    // a few MatCfg alive at once, this is not expected to be a problem in
    // practice. Notice also the cloneThinned() method below.

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


    /////////////////////////////////////////////////////////////////////////////
    // Possible parameters and their meaning:                                  //
    /////////////////////////////////////////////////////////////////////////////
    //
    // temp........: [ double, fallback value is -1.0 ]
    //               Temperature of material in Kelvin. The special value of
    //               -1.0 implies that the info factory should determine an
    //               appropriate value based on the input data, falling back to
    //               293.15K when the input has not special preference (example:
    //               an S(alpha,beta) scattering kernel might be valid at a
    //               particular temperature value only).
    //               [ Recognised units: "K", "C", "F" ]
    //
    // coh_elas....: [ bool, fallback value is true ]
    //               If enabled, scatter factories will include coherent elastic
    //               (i.e. Bragg diffraction) components for crystalline
    //               materials. See also pseudo-parameters "elas" and "bragg"
    //               below which can both be used to change the value of the
    //               coh_elas parameter.
    //
    // incoh_elas....: [ bool, fallback value is true ]
    //               If enabled, scatter factories will include incoherent elastic
    //               components for crystalline materials. See also pseudo-parameters
    //               "elas" and "bkgd" below which can both be used to change the value
    //               of the incoh_elas parameter.
    //
    // inelas........: [ string, fallback value is "auto" ]
    //               Influence inelastic scattering models chosen by scatter
    //               factories. The default value of "auto" leaves the choice up
    //               to scatter factories, and a value of "none", "0", "false",
    //               or "sterile", all disables inelastic scattering.  The
    //               standard NCrystal scatter factory currently supports
    //               additional values: "external", "dyninfo", "vdosdebye", and
    //               "freegas", and internally the "auto" mode will simply
    //               select the first possible of those in the listed order
    //               (falling back to "none" when nothing is possible). Note
    //               that "external" is not possible for .ncmat files (but is
    //               for .nxs files), and when using .laz/.lau only "none" is
    //               possible.  The "dyninfo" mode will simply base modelling on
    //               whatever dynamic information is available for each element
    //               in the input data. The "vdosdebye" and "freegas" modes
    //               overrides this, and force those models for all elements if
    //               possible.  The "external" mode implies usage of an
    //               externally provided cross-section curve and an
    //               isotropic-elastic scattering model. See also the parameter
    //               "vdoslux" and pseudo-parameter "bkgd" below.
    //
    // dcutoff.....: [ double, fallback value is 0 ]
    //               D-spacing cutoff in Angstrom. Crystal planes with spacing
    //               below this value will not be created. The special setting
    //               dcutoff=0 causes the code to attempt to select an
    //               appropriate threshold automatically, and the special
    //               setting dcutoff=-1 means that HKL lists should not be
    //               created at all (often used with bragg=false).
    //               [ Recognised units: "Aa", "nm", "mm", "cm", "m" ]
    //
    // packfact....: [ double, fallback value is 1.0 ]
    //               Packing factor which can be less than 1.0 for powders,
    //               which can thus be modelled as polycrystals with reduced
    //               density (not to be confused with the *atomic* packing
    //               factor).
    //
    // mos.........: [ double, no fallback value ]
    //               Mosaic FWHM spread in mosaic single crystals, in radians.
    //               [ Recognised units: "rad", "deg", "arcmin", "arcsec" ]
    //
    // dir1........: [ special, no fallback value ]
    //               Used to specify orientation of single crystals, by
    //               providing both a vector in the crystal frame and the lab
    //               frame, as explained in more detail in the file
    //               NCSCOrientation.hh. If the six numbers in the two vectors
    //               are respectively (c1,c2,c3) and (l1,l2,l3), this is
    //               specified as: "dir1=@crys:c1,c2,c3@lab:l1,l2,l3" If
    //               (c1,c2,c3) are points in hkl space, simply use "@crys_hkl:"
    //               instead, as in: "dir1=@crys_hkl:c1,c2,c3@lab:l1,l2,l3"
    //
    // dir2........: [ special, no fallback value ]
    //               Similar to dir1, but the direction might be modified
    //               slightly in case of imprecise input, up to the value of
    //               dirtol).  See the file NCSCOrientation.hh for more details.
    //
    // dirtol......: [ double, fallback value is 0.0001 ]
    //               Tolerance parameter for the secondary direction of the
    //               single crystal orientation, here in radians. See the file
    //               NCSCOrientation.hh for more details.
    //               [ Recognised units: "rad", "deg", "arcmin", "arcsec" ]
    //
    // lcaxis......: [ vector, no fallback value ]
    //               Used to specify symmetry axis of anisotropic layered
    //               crystals with a layout similar to pyrolytic graphite, by
    //               providing the axis in lattice coordinates using a format
    //               like "0,0,1" (does not need to be normalised). Specifying
    //               this parameter along with an orientation (see dir1, dir2
    //               and dirtol parameters) will result in a specialised single
    //               crystal scatter model being used.
    //
    //
    /////////////////////////////////////////////////////////////////////////////
    // Options mainly of interests to experts and NCrystal developers:
    //
    // dcutoffup...: [ double, fallback value is infinity ]
    //               Like dcutoff, but representing an upper cutoff instead
    //               [ Recognised units: "Aa", "nm", "mm", "cm", "m" ]
    //
    // infofactory.: [ string, fallback value is "" ]
    //               By supplying the name of an NCrystal factory, this
    //               parameter can be used by experts to circumvent the usual
    //               factory selection algorithms and instead choose the factory
    //               for creating NCrystal::Info instances directly.
    //               (TODO: Mention how to set flags like expandhkl)
    //
    // scatfactory.: [ string, fallback value is "" ]
    //               Similar to infofactory, this parameter can be used to
    //               directly select factory with which to create
    //               NCrystal::Scatter instances. As a special feature (needed
    //               for plugin development), factories can be excluded by
    //               adding them with a "!" in front of their name. Multiple
    //               entries can be added by separating them with an "@" sign
    //               (but at most one non-excluded entry can appear).
    //
    // absnfactory.: [ string, fallback value is "" ]
    //               Similar to infofactory, this parameter can be used to
    //               directly select factory with which to create
    //               NCrystal::Absorption instances. As a special feature (needed
    //               for plugin development), factories can be excluded by
    //               adding them with a "!" in front of their name. Multiple
    //               entries can be added by separating them with an "@" sign
    //               (but at most one non-excluded entry can appear).
    //
    // mosprec.....: [ double, fallback value is 1.0e-3 ]
    //               Approximate relative precision in implementation of mosaic
    //               model in single crystals. Affects both approximations used
    //               and truncation range of Gaussian. Values must be in the
    //               range [1e-7,1e-1].
    //
    // lcmode......: [ int, fallback value is 0 ]
    //               Choose which modelling is used for layered crystals (has no
    //               effect unless lcaxis is also set). The default value
    //               indicates the recommended model, which is both fast and
    //               accurate. A positive value triggers a very slow but simple
    //               reference model, in which n=lcmode crystallite orientations
    //               are sampled internally (the model is accurate only when n
    //               is very high). A negative value triggers a different (and
    //               multi-thread unsafe!) model in which each crossSection call
    //               triggers a new selection of n=-lcmode randomly oriented
    //               crystallites.
    //
    // sccutoff....: [ double, fallback value is 0.4Aa ]
    //               Single-crystal d-spacing cutoff in Angstrom. When creating
    //               single-crystal scatterers, crystal planes with spacing
    //               below this value will be modelled as having an isotropic
    //               mosaicity distribution. This usually results in very great
    //               computational speedups for neutrons at wavelengths below
    //               2*sccutof. The tradeoff is in principle an incorrect angular
    //               dependency when calculating the cross-section for scattering
    //               on *these* individual planes, but in practice the net effect
    //               is usually not particularly significant due to the very large
    //               number of very weak planes affected. Setting sccutoff=0
    //               naturally disables this approximation.
    //               [ Recognised units: "Aa", "nm", "mm", "cm", "m" ]
    //
    // vdoslux.....: [ int, fallback value is 3 ]
    //               Setting affecting "luxury" level when expanding phonon
    //               spectrums (VDOS) into scattering kernels, affecting things
    //               like number of (alpha,beta) grid points in the resulting
    //               kernel and what energy range is covered by the kernel. In
    //               very rough terms, the levels have the following approximate
    //               impact (exact impact depends on both the given VDOS as well
    //               as values of other configuration parameters):
    //                  0 : Extremely crude, 100x50 grid, Emax=0.5eV
    //                      Costs 0.1MB mem, 0.02s init time
    //                  1 : Crude, 200x100 grid, Emax=1eV
    //                      Costs 0.5MB mem, 0.04s init time
    //                  2 : Decent, 400x200 grid, Emax=3eV
    //                      Costs 2MB mem, 0.08s init time
    //                  3 : Good, 800x400 grid, Emax=5eV
    //                      Costs 8MB mem, 0.2s init time
    //                  4 : Very good, 1600x800 grid, Emax=8eV
    //                      Costs 30MB mem, 0.8s init time
    //                  5 : Extremely good, 3200x1600 grid, Emax=12eV
    //                      Costs 125MB mem, 5s init time
    //               Levels 2-4 are intended for normal usage, level 5 as a
    //               validation reference, while levels 0 and 1 are intended for
    //               the case where the input VDOS data is anyway just a crude
    //               estimate. For cases where NCrystal has no actual VDOS input
    //               data provided and instead generates an idealised Debye
    //               spectrum on the fly based on the Debye temperature, the
    //               vdoslux level actually used will be 3 less than the one
    //               specified in this variable (but at least 0).
    //
    // atomdb......: [ string, fallback value is "" ]
    //               Modify atomic definitions if supported by the info factory
    //               (in practice this is unlikely to be supported by anything
    //               except NCMAT data). The string must follow a syntax
    //               identical to that used in @ATOMDB sections of NCMAT file
    //               (cf. ncmat_doc.md), with a few exceptions explained here:
    //               First of all it is allowed to use semicolons (':') to
    //               divide words (in fact all semicolons will simply be
    //               replaced with spaces during evaluation), which makes it
    //               possible to write cfg-strings without spaces (this might
    //               occasionally be useful, e.g. on the command line). Next,
    //               '@' characters play the role of line separators. Finally,
    //               when used with an NCMAT file that already includes an
    //               internal @ATOMDB section, the effect will essentially be to
    //               combine the two sections by appending the atomdb lines from
    //               the cfg parameter to the lines already present in the input
    //               data. The exception is the case where the cfg parameter
    //               contains an initial line with the single word "nodefaults",
    //               the effect of which will always be the same as if it was
    //               placed on the very first line in the @ATOMDB section
    //               (i.e. NCrystal's internal database of elements and isotopes
    //               will be ignored).


    /////////////////////////////////////////////////////////////////////////////
    // Special pseudo-parameters available for usage in configuration strings
    // (and only there), for convenience and backwards compatiblity:
    //
    // bragg..: This is simply an alias for the "coh_elas" parameter.
    // elas...: Assigning a boolean value to this pseudo-parameter will change
    //          both the "coh_elas" and "incoh_elas" parameters at once.
    // bkgd...: Assigning "0" or "none" to this pseudo-parameter will result in
    //          coh_elas and inelas parameters being set to false and "none"
    //          respectively. No other values are accepted as this parameter
    //          exists purely for backwards compatiblity reasons. Users should
    //          now use the "inelas" parameter to affect the choice of inelastic
    //          model, and the "incoh_elas" parameter to toggle incoherent
    //          elastic scattering.

    /////////////////////////////////////////////////////////////////////////////
    // Methods for setting parameters:                                         //
    /////////////////////////////////////////////////////////////////////////////
    //
    //Directly set from C++ code:
    void set_temp( Temperature );
    void set_temp( double );
    void set_dcutoff( double );
    void set_dcutoffup( double );
    void set_packfact( double );
    void set_mos( MosaicityFWHM );
    void set_mosprec( double );
    void set_sccutoff( double );
    void set_dirtol( double );
    void set_coh_elas( bool );
    void set_incoh_elas( bool );
    void set_inelas( const std::string& );
    void set_infofactory( const std::string& );
    void set_scatfactory( const std::string& );
    void set_absnfactory( const std::string& );
    void set_lcmode( int );
    void set_vdoslux( int );
    void set_atomdb( const std::string& );
    void set_lcaxis( const LCAxis& );
    //
    //Special setter method, which will set all orientation parameters based on
    //an SCOrientation object:
    void setOrientation( const SCOrientation& );
    void set_dir1( const HKLPoint&, const LabAxis& );
    void set_dir1( const CrystalAxis&, const LabAxis& );
    void set_dir2( const HKLPoint&, const LabAxis& );
    void set_dir2( const CrystalAxis&, const LabAxis& );
    std::pair<Variant<CrystalAxis,HKLPoint>,LabAxis> get_dir1() const;
    std::pair<Variant<CrystalAxis,HKLPoint>,LabAxis> get_dir2() const;
    //
    // Set parameters from a string, using the same format as that supported by
    // the constructor, e.g. "par1=val1;...;parn=valn":
    void applyStrCfg( const std::string&  );

    /////////////////////////////////////////////////////////////////////////////
    // Methods for accessing parameters or derived information:                //
    /////////////////////////////////////////////////////////////////////////////
    //
    // Directly access from C++
    Temperature get_temp() const;
    double get_dcutoff() const;
    double get_dcutoffup() const;
    double get_packfact() const;
    MosaicityFWHM get_mos() const;
    double get_mosprec() const;
    double get_sccutoff() const;
    double get_dirtol() const;
    LCAxis get_lcaxis() const;
    const std::string& get_scatfactory() const;
    const std::string& get_absnfactory() const;
    int  get_lcmode() const;
    int  get_vdoslux() const;
    const std::string& get_atomdb() const;
    const std::vector<VectS>& get_atomdb_parsed() const;

    bool get_coh_elas() const;
    bool get_incoh_elas() const;
    const std::string& get_inelas() const;

    //infofactory option decoded:
    std::string get_infofact_name() const;
    bool get_infofactopt_flag(const std::string& name) const;
    double get_infofactopt_dbl(const std::string& name, double defval) const;
    int get_infofactopt_int(const std::string& name, int defval) const;

    // Specialised getters for derived information:
    bool isSingleCrystal() const;//true if mos or orientation parameters are set
    bool isPolyCrystal() const;//same as !isSingleCrystal()
    SCOrientation createSCOrientation() const;//Create and return a new SCOrientation object based cfg.
    bool isLayeredCrystal() const;//true if lcaxis parameter is set

    //Test if was constructed with ";ignorefilecfg" keyword:
    bool ignoredEmbeddedConfig() const;

    //Validate infofactory flags and options to prevent silently ignoring
    //unused options. Call only from *selected* factory, to throw BadInput in
    //case of unknown options:
    void infofactopt_validate(const std::set<std::string>& allowed_opts) const;

    //Serialise in various forms (note that if the MatCfg object was not
    //constructed from a text string, the string returned by toStrCfg(true) will
    //of course not be useful for creating a new equivalent MatCfg object).
    std::string toStrCfg( bool include_datafile = true ) const;
    std::string toEmbeddableCfg() const;//Produces a string like "NCRYSTALMATCFG[temp=500.0;dcutoff=0.2]"
                                        //which can be embedded in data files.
    void dump( std::ostream &out, bool add_endl=true ) const;

    /////////////////////////////////////////////////////////////////////////////
    // Associated text data (content and meta data of data file)               //
    /////////////////////////////////////////////////////////////////////////////
    //
    // Both content and meta data are available, unless the MatCfg object was
    // "thinned" (see below), in which case the textData()/textDataSP() methods
    // should not be used (but the textDataUID()/getDataType methods are always
    // OK). However, normally MatCfg objects are not thinned and most code
    // should simply use the methods.

    const TextData& textData() const;
    TextDataSP textDataSP() const;

    //UID of text data:
    const TextDataUID textDataUID() const;

    //Data type of text data (usually the extension of the data file). Calling
    //this method should be the only manner of determining the format of the
    //associated textData:
    const std::string& getDataType() const;

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

    //Create a "thinned" clone which does not keep a strong reference to the
    //associated TextData object, but which still retains all other information
    //(including the textDataUID) needed for comparisons. Thus, such a clone is
    //suitable for e.g. cache keys which must be kept around for comparison but
    //which should not retain strong references to potentially large text
    //data. Accessing the textData()/textDataSP() methods on a thinned object
    //will result in an exception.
    MatCfg cloneThinned() const;
    bool isThinned() const;

    /////////////////////////////////////////////////////////////////////////////
    // Interface for NCrystal factory infrastructure:                          //
    /////////////////////////////////////////////////////////////////////////////
    //
    // Advanced interface, which is intended solely for use by the NCrystal
    // factory infrastructure, for validation, factory selection, and to avoid
    // unnecessarily re-initialisation of expensive objects.

    //Create const-view wrapper of MatCfg with access to the restricted set of
    //parameters which is allowed to use in info factories:
    MatInfoCfg createInfoCfg() const;

    //Verify that the parameter values are not inconsistent or
    //incomplete. Throws BadInput exception if they are. This will automatically
    //be invoked by the factory infrastructure:
    void checkConsistency() const;

    //Convenience interface for setting/decoding scatfactory+absnfactory parameters:
    struct FactRequested {
      std::string specific;
      std::set<std::string> excluded;
    };
    FactRequested get_scatfactory_parsed() const;
    FactRequested get_absnfactory_parsed() const;
    void set_scatfactory(const FactRequested& );
    void set_absnfactory(const FactRequested& );

    //Unspecified ordering (for usage as map keys):
    bool operator<(const MatCfg&) const;

  private:
    struct constructor_args : private MoveOnly {
      OptionalTextDataSP td; std::string pars, origfn;
    };
    explicit MatCfg( constructor_args&& );
    struct from_raw_t {};
    explicit MatCfg( from_raw_t, std::string&& data, std::string pars, std::string ext );
    const std::string& get_infofactory() const;//undecoded, internal usage only
    struct Impl;
    COWPimpl<Impl> m_impl;
    OptionalTextDataSP m_textDataSP;
    friend class MatInfoCfg;
  };

  class MatInfoCfg {
  public:
    //Reduced const-view wrapper of MatCfg, restricting access to the parameters
    //that are allowed only in Info-factories. The comparison operator likewise
    //ignores all other parameters.
    const TextDataUID textDataUID() const;
    const std::string& getDataType() const;
    const TextData& textData() const;
    TextDataSP textDataSP() const;
    Temperature get_temp() const;
    double get_dcutoff() const;
    double get_dcutoffup() const;
    const std::string& get_atomdb() const;
    const std::vector<VectS>& get_atomdb_parsed() const;
    std::string get_infofact_name() const;
    bool get_infofactopt_flag(const std::string& name) const;
    double get_infofactopt_dbl(const std::string& name, double defval) const;
    int get_infofactopt_int(const std::string& name, int defval) const;
    void checkConsistency() const;
    void dump( std::ostream &out, bool add_endl=true ) const;
    void infofactopt_validate(const std::set<std::string>& ao) const;
    std::string toStrCfg( bool include_datafile = true ) const;
    MatInfoCfg clone() const;
    MatInfoCfg cloneThinned() const;
    bool isThinned() const;
    MatInfoCfg(const MatCfg&);
    MatInfoCfg(MatCfg&&);
    MatInfoCfg(const MatInfoCfg&) = default;
    MatInfoCfg& operator=(const MatInfoCfg&) = default;
    MatInfoCfg( MatInfoCfg&& ) = default;
    MatInfoCfg& operator=(MatInfoCfg&&) = default;
    bool operator<(const MatInfoCfg&) const;
  private:
    MatCfg m_cfg;
  };

  inline std::ostream& operator<< (std::ostream&, const MatCfg&);
  inline std::ostream& operator<< (std::ostream&, const MatInfoCfg&);

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCrystal {
  inline void MatCfg::set_temp( double t ) { set_temp(Temperature{t}); }
  inline MatCfg MatCfg::clone() const { return MatCfg(*this); }
  inline bool MatCfg::isThinned() const { return m_textDataSP == nullptr; }
  inline const TextData& MatCfg::textData() const { return textDataSP(); }
  inline const TextDataUID MatInfoCfg::textDataUID() const { return m_cfg.textDataUID(); }
  inline const std::string& MatInfoCfg::getDataType() const { return m_cfg.getDataType(); }
  inline const TextData& MatInfoCfg::textData() const { return m_cfg.textData(); }
  inline TextDataSP MatInfoCfg::textDataSP() const { return m_cfg.textDataSP(); }
  inline Temperature MatInfoCfg::get_temp() const { return m_cfg.get_temp(); }
  inline double MatInfoCfg::get_dcutoff() const { return m_cfg.get_dcutoff(); }
  inline double MatInfoCfg::get_dcutoffup() const { return m_cfg.get_dcutoffup(); }
  inline const std::string& MatInfoCfg::get_atomdb() const { return m_cfg.get_atomdb(); }
  inline const std::vector<VectS>& MatInfoCfg::get_atomdb_parsed() const { return m_cfg.get_atomdb_parsed(); }
  inline std::string MatInfoCfg::get_infofact_name() const { return m_cfg.get_infofact_name(); }
  inline bool MatInfoCfg::get_infofactopt_flag(const std::string& name) const { return m_cfg.get_infofactopt_flag(name); }
  inline double MatInfoCfg::get_infofactopt_dbl(const std::string& name, double defval) const { return m_cfg.get_infofactopt_dbl(name,defval); }
  inline int MatInfoCfg::get_infofactopt_int(const std::string& name, int defval) const { return m_cfg.get_infofactopt_int(name,defval); }
  inline void MatInfoCfg::checkConsistency() const { m_cfg.checkConsistency(); }
  inline void MatInfoCfg::infofactopt_validate(const std::set<std::string>& ao) const { m_cfg.infofactopt_validate(ao); }
  inline MatInfoCfg MatInfoCfg::clone() const { return m_cfg.clone(); }
  inline MatInfoCfg MatInfoCfg::cloneThinned() const { return m_cfg.cloneThinned(); }
  inline MatInfoCfg::MatInfoCfg(const MatCfg& cfg) : m_cfg(cfg) {}
  inline MatInfoCfg::MatInfoCfg(MatCfg&& cfg) : m_cfg(std::move(cfg)) {}
  inline bool MatInfoCfg::isThinned() const { return m_cfg.isThinned(); }
  inline std::ostream& operator<< (std::ostream& os, const MatCfg& cfg) { cfg.dump(os,false); return os; }
  inline std::ostream& operator<< (std::ostream& os, const MatInfoCfg& cfg) { cfg.dump(os,false); return os; }
}

#endif
