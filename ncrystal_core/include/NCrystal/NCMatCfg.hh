#ifndef NCrystal_MatCfg_hh
#define NCrystal_MatCfg_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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
#include <set>
#include <ostream>

namespace NCrystal {

  class Info;
  class SCOrientation;

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
    //
    // Note that it is also possible to set the (initial) values of parameters,
    // by embedding a statement like NCRYSTALMATCFG[temp=500.0;dcutoff=0.2] into
    // the datafile itself. These parameters can of course still be overridden,
    // and such embedded configuration can be ignored entirely by appending
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

    /////////////////////////////////////////////////////////////////////////////
    // Possible parameters and their meaning:                                  //
    /////////////////////////////////////////////////////////////////////////////
    //
    // temp........: [ double, fallback value is 293.15 ]
    //               Temperature of material in Kelvin.
    //               [ Recognised units: "K", "C", "F" ]
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
    //               Mosaic spread in mosaic single crystals, in radians.
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
    // bragg.......: [ bool, fallback value is true ]
    //               If enabled, scatter factories will include Bragg
    //               diffraction components.
    //
    // bkgd........: [ string, fallback value is "best" ]
    //               Influence background model chosen by scatter factories. The
    //               default value of "best" implies that they should pick the
    //               most realistic one available, while a value of "none"
    //               prevents background components from being added. Depending
    //               on the model in question, it is possible to specify
    //               additional options by appending them using ':' and with '@'
    //               signs for assignments (see examples below).
    //
    //               The default NCrystal factories currently support "external"
    //               and "phonondebye" models, both of which will model energy
    //               transfers as fully thermalising if temperature is known and
    //               otherwise elastic. This can be overridden using
    //               "thermalise" or "elastic" flags, and the phonon expansion
    //               order for the phonondebye model can be controlled with
    //               "nphonon". Examples: "external:thermalising",
    //               "phonondebye:elastic", "phonondebye:thermalising:nphonon@20"
    //
    // lcaxis......: [ vector, no fallback value ]
    //               Used to specify symmetry axis of anisotropic layered
    //               crystals with a layout similar to pyrolythic graphite, by
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
    //
    // scatfactory.: [ string, fallback value is "" ]
    //               Similar to infofactory, this parameter can be used to
    //               directly select factory with which to create
    //               NCrystal::Scatter instances.
    //
    // absnfactory.: [ string, fallback value is "" ]
    //               Similar to infofactory, this parameter can be used to
    //               directly select factory with which to create
    //               NCrystal::Absorption instances.
    //
    // mosprec...: [ double, fallback value is 1.0e-3 ]
    //               Approximate relative precision in implementation of mosaic
    //               model in single crystals. Affects both approximations used
    //               and truncation range of Gaussian.
    //
    // lcmode....: [ int, fallback value is 0 ]
    //               Choose which modelling is used for layered crystals (has no
    //               effect unless lcaxis is also set). The default value
    //               indicates the recommended model, which is both fast and
    //               accurate. A positive value triggers a very slow but simple
    //               reference model, in which n=lcmode crystallite orientations
    //               are sampled internally (the model is accurate only when n
    //               is very high). A negative value triggers a different model
    //               in which each crossSection call triggers a new selection of
    //               n=-lcmode randomly oriented crystallites.
    //
    // sccutoff.....: [ double, fallback value is 0.5Aa ]
    //               Single-crystal d-spacing cutoff in Angstrom. When creating
    //               single-crystal scatterers, crystal planes with spacing
    //               below this value and Fsquared less than 10% of the maximum
    //               Fsquared in the crystal will be modelled as having a
    //               isotropic mosaicity distribution. This usually results in
    //               very great computational speedups for neutrons at
    //               wavelengths below 2*sccutof. The tradeoff is in principle
    //               incorrect angular dependency when calculating the
    //               cross-section for scattering on *these* individual planes,
    //               but in practice the net effect is usually not particularly
    //               significant due to the very large number of very weak
    //               planes affected. Setting sccutoff=0 naturally disables this
    //               approximation.
    //               [ Recognised units: "Aa", "nm", "mm", "cm", "m" ]


    /////////////////////////////////////////////////////////////////////////////
    // Methods for setting parameters:                                         //
    /////////////////////////////////////////////////////////////////////////////
    //
    //Directly set from C++ code:
    void set_temp( double );
    void set_dcutoff( double );
    void set_dcutoffup( double );
    void set_packfact( double );
    void set_mos( double );
    void set_mosprec( double );
    void set_sccutoff( double );
    void set_dirtol( double );
    void set_overridefileext( const std::string& );
    void set_bragg( bool );
    void set_bkgd( const std::string& );
    void set_infofactory( const std::string& );
    void set_scatfactory( const std::string& );
    void set_absnfactory( const std::string& );
    void set_lcmode( int );
    //
    //Special setter method, which will set all orientation parameters based on
    //an SCOrientation object:
    void setOrientation( const SCOrientation& );

    void set_dir1( bool crystal_dir_is_point_in_hkl_space,
                   const double (&crystal_direction)[3],
                   const double (&lab_direction)[3] );
    void set_dir2( bool crystal_dir_is_point_in_hkl_space,
                   const double (&crystal_direction)[3],
                   const double (&lab_direction)[3] );
    void get_dir1( bool& crystal_dir_is_point_in_hkl_space,
                   double (&crystal_direction)[3],
                   double (&lab_direction)[3] );
    void get_dir2( bool& crystal_dir_is_point_in_hkl_space,
                   double (&crystal_direction)[3],
                   double (&lab_direction)[3] );
    void set_lcaxis( const double (&axis)[3] );
    //
    // Set parameters from a string, using the same format as that supported by
    // the constructor, e.g. "par1=val1;...;parn=valn":
    void applyStrCfg( const std::string&  );

    /////////////////////////////////////////////////////////////////////////////
    // Methods for accessing parameters or derived information:                //
    /////////////////////////////////////////////////////////////////////////////
    //
    // Directly access from C++
    double get_temp() const;
    double get_dcutoff() const;
    double get_dcutoffup() const;
    double get_packfact() const;
    double get_mos() const;
    double get_mosprec() const;
    double get_sccutoff() const;
    double get_dirtol() const;
    bool get_bragg() const;
    void get_lcaxis( double (&axis)[3] ) const;
    const std::string& get_overridefileext() const;
    const std::string& get_scatfactory() const;
    const std::string& get_absnfactory() const;
    int  get_lcmode() const;

    //Bkgd option decoded:
    std::string get_bkgd_name() const;
    bool get_bkgdopt_flag(const std::string& name) const;
    double get_bkgdopt_dbl(const std::string& name, double defval) const;
    int get_bkgdopt_int(const std::string& name, int defval) const;

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

    //Validate bkgd/infofactory flags and options to prevent silently ignoring
    //unused options. Call only from *selected* factory, to throw BadInput in
    //case of unknown options:
    void bkgdopt_validate(const std::set<std::string>& allowed_opts) const;
    void infofactopt_validate(const std::set<std::string>& allowed_opts) const;

    //Datafile (never decode extension by hand):
    const std::string& getDataFile() const;//Resolved path to the datafile
    const std::string& getDataFileAsSpecified() const;//Path to the datafile as specified in the constructor
    const std::string& getDataFileExtension() const;//Extension of datafile (actual unless overridefileext is set)

    //Serialise in various forms:
    std::string toStrCfg( bool include_datafile = true, const std::set<std::string> * only_parnames = 0 ) const;
    std::string toEmbeddableCfg() const;//Produces a string like "NCRYSTALMATCFG[temp=500.0;dcutoff=0.2]"
                                        //which can be embedded in data files.
    void dump( std::ostream &out, bool add_endl=true ) const;

    //Test if was constructed with ";ignorefilecfg" keyword:
    bool ignoredEmbeddedConfig() const;


    /////////////////////////////////////////////////////////////////////////////
    // Copy/assign/clone/destruct                                              //
    /////////////////////////////////////////////////////////////////////////////
    //
    // Copy/assignment/cloning is allowed and is a priori cheap, since internal
    // data structures are shared until modified (aka copy-on-write):
    MatCfg(const MatCfg&);
    MatCfg& operator=(const MatCfg&);
    MatCfg clone() const { return MatCfg(*this); }
    //A completely "clean" clone, unconnected to the original and no internal caches set is obtainable with:
    MatCfg cloneUnshared() const { MatCfg c(*this); c.cow(); return c; }
    //
    //Destructor:
    ~MatCfg();


    /////////////////////////////////////////////////////////////////////////////
    // Interface for NCrystal factory infrastructure:                          //
    /////////////////////////////////////////////////////////////////////////////
    //
    // Advanced interface, which is intended solely for use by the NCrystal
    // factory infrastructure, in order to avoid unnecessarily re-initialising
    // expensive objects.
    //
    // Monitor parameter accesses (does not take ownership of the spy and
    // parameter modifications will be forbidden while such monitoring is in
    // place):
    struct AccessSpy {
      virtual ~AccessSpy(){};
      virtual void parAccessed(const std::string& parname) = 0;
    };
    bool hasAccessSpy(AccessSpy*) const;
    void addAccessSpy(AccessSpy*) const;
    void removeAccessSpy(AccessSpy*) const;
    //
    //Verify that the parameter values are not inconsistent or
    //incomplete. Throws BadInput exception if they are. This will automatically
    //be invoked by the factory infrastructure:
    void checkConsistency() const;
    //
    //returns a string like "par1=val1;...;parn=valn" for the specified
    //parameters only (using "<>" as value for unset parameters):
    void getCacheSignature( std::string& signature,
                            const std::set<std::string>& parameters ) const;

  private:
    const std::string& get_bkgd() const;//undecoded, internal usage only
    const std::string& get_infofactory() const;//undecoded, internal usage only
    struct Impl;
    Impl* m_impl;
    void cow();

  };

  inline std::ostream& operator<< (std::ostream& s, const MatCfg& cfg) { cfg.dump(s,false); return s; }

}

#endif
