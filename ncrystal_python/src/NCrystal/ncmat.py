
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2025 NCrystal developers                                   ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################

"""

Utilities for NCMAT creation, most notably the NCMATComposer helper class which
provides a pythonic and flexible object oriented interface for creating NCMAT
data representing a given material, and for verifying, writing, registering, or
loading said data.

"""

class NCMATComposer:
    """Helper class for composing an NCrystal material.

    One can either start with an empty instance, or bootstrap the process by
    starting the material from an existing NCMAT file, an NCrystal cfg-string,
    an NCrystal.Info object, a CIF file, or from a variety of files (anything
    that can be loaded with ASE and contains a suitable crystal
    structure). Depending on the specifics of the problem, one will then use
    various methods on the object to add or modify material information as
    desired, before ultimately the create_ncmat() method can be used to create
    the resulting NCMAT data. As a convenience the .write() and .register_as()
    methods can be used to write this data as a physical or virtual file
    respectively, or the .load() method can be used to directly load the
    data. Finally, the .plot_xsect() and .inspect() methods can also be used to
    quickly load and investigate the material. It might also be useful to print
    the output of .create_ncmat(), to understand what material has been created.

    For a more in-depth discussion and usage examples, please refer to the
    example Jupyter-Lab notebooks in the
    https://github.com/mctools/ncrystal-notebooks repository.

    Note that the "labels" used in the various methods, to identify different
    atoms in the materials are custom user labels. They are simply a way to
    provide a handle associated with a particular atom role or species in the
    material. For many simple materials, all labels will simply have the name of
    the element or isotope ("Al", "Li6", "D", ...), and each type of element and
    isotope will only appear once in the material. That is certainly OK, for
    such simple materials, and in these cases there is no need to concern
    one-self with the atomic composition of each label, since it can be inferred
    from the label itself. However, sometimes different atoms of the same type
    will appear in more than one role in the material (i.e. there might be one
    group of loosely bound aluminium atoms, while a different group of aluminium
    atoms are more tightly bound). Or, it might be that a given label is
    actually associated with an actual mix of elements of isotopes (i.e. a label
    might be associated with a particular isotopic enrichment, or a
    multi-element mixture). In these cases, each label should have it's precise
    composition defined with a call to the .set_composition(..) method.

    In other words: Unless the label happens to coincide with the name of an
    element or isotope ("Al", "Li6", "D", ...), its atomic composition MUST be
    defined by a call to the .set_composition(..) method. And if a given element
    has more than one role in the material, it must of course also have more
    than one label - one for each role (e.g. "tightAl" and "looseAl").

    To ensure high data quality and reliability, all crystalline materials will
    by default have their crystal structure verified for consistency with the
    indicated spacegroup, before it is converted to NCMAT data. For that reason,
    the third-party module "spglib" must be installed in order for this
    verification to take place (note that the .refine_crystal_structure() method
    can also be invoked before the material is used, to detect the spacegroup in
    if it was not provided). Fortunately, "spglib" is available both on PyPI and
    in the conda-forge conda channel.

    """

    def __init__(self, data = None, fmt = None, quiet = False, plotlabel = None ):
        """Initialise. Either an empty instance, or if data (and possibly fmt)
        is provided, from a variety of data. Typically, the fmt parameter can be
        left out, but occasionally it might be needed to specify it
        explicitly. Supported fmt's are "cif", "cfgstr", "ncmat", (for NCMAT
        data), "info" (for NCrystal.Info or LoadedMaterial objects), ase (for
        ASE structures). Finally, the format "via_ase" is special, since it will
        always process the input with ASE, even if another more obvious format
        exists (i.e. CIF files will be processed via ASE first, not passed along
        directly to the NCrystal.CIFLoader). As an alternative to the "fmt"
        parameter, one can also invoke the various .from_xxx(..)  static
        methods. This has the potential advantage of allowing the specification
        of parameters dedicated only to that input data type.

        If quiet=True, no informative messages will be emitted, and the plotlabel
        will be passed along to .set_plotlabel().

        """
        from ._ncmatimpl import NCMATComposerImpl as Impl
        self.__impl = ( data if isinstance( data, Impl )
                        else ( data.__impl if isinstance( data, NCMATComposer )
                               else Impl( data = data, fmt = fmt, quiet = quiet ) ) )
        if plotlabel:
            self.set_plotlabel(plotlabel)

    @staticmethod
    def from_cif( cifsrc, quiet=False, mp_apikey = None,
                  uiso_temperature = None,
                  override_spacegroup = None, **kwargs ):
        """Initialise from CIF data (anything suitable for the
        NCrystal.cifutils.CIFSource class). Internally, the loading is carried
        out by instantiating a CIFLoader instance (from the NCrystal.cifutils
        module), and calling it's .create_ncmat_composer() method. Any provided
        arguments will be passsed along to the call to .create_ncmat_composer(),
        except for the override_spacegroup argument which is used when the
        CIFLoader is instantiated.

        For advanced work, it might be better to use the CIFLoader class
        directly.
        """
        from ._ncmatimpl import NCMATComposerImpl as Impl
        return NCMATComposer( Impl.from_cif( cifsrc=cifsrc, quiet=quiet,
                                             uiso_temperature = uiso_temperature,
                                             override_spacegroup = override_spacegroup,
                                             mp_apikey = mp_apikey, **kwargs ) )

    @staticmethod
    def from_cfgstr( cfgstr ):
        """Loads an Info object from the given cfg-string, and pass it on to
        .from_info(..)."""
        from ._ncmatimpl import NCMATComposerImpl as Impl
        return NCMATComposer( Impl.from_cfgstr( cfgstr = cfgstr ) )

    @staticmethod
    def from_ase( ase_obj, ase_format = None, *, quiet = False, **kwargs ):
        """Initialise from an ASE object. This is actually done internally by
        first converting the ASE object to CIF data, and then subsequently
        passing it on to the .from_cif(..) method (any **kwargs will be
        passed along to the .from_cif call). If the provided ase_obj is not
        actually already an ASE object, an attempt is done to load it into one
        via the ase.io.read module (mimicking the fmt="via_ase" behaviour of the
        __init__ method). In that case, any value in the ase_format parameter
        will be passed along to ase.io.read's format parameter.
        """
        from ._ncmatimpl import NCMATComposerImpl as Impl
        return NCMATComposer( Impl.from_ase( ase_obj = ase_obj, quiet = quiet,
                                             ase_format = ase_format, **kwargs ) )

    @staticmethod
    def from_info( info_obj ):
        """Initialise from provided Info (or LoadedMaterial) object.

        In case the material is a gas, the temperature of the resulting
        NCMATComposer will for safety be locked (due to the strong temperature
        dependency of densities for gasses).
        """
        from ._ncmatimpl import NCMATComposerImpl as Impl
        return NCMATComposer( Impl.from_info( info_obj = info_obj ) )

    @staticmethod
    def from_ncmat( data ):
        """Loads an Info object from the given NCMAT data, and passes it on to
        .from_info(..). If the NCMAT data has no newlines and doesn't start with
        'NCMAT', it will actually be passed directly on to .from_cfgstr(..)
        instead.
        """
        from ._ncmatimpl import NCMATComposerImpl as Impl
        return NCMATComposer( Impl.from_ncmat( data ) )

    @staticmethod
    def from_hfg( spec,
                  formula, *,
                  density,
                  title,
                  debyetemp = 400.0,
                  verbose = True,
                  notrim = False ):
        """Constructs an amorphous hydrogen-rich material via the hfg2ncmat
        function from the NCrystal.hfg2ncmat module. Refer to that function for
        usage instructions.
        """
        from .hfg2ncmat import _default_debye_temp, hfg2ncmat
        from ._ncmatimpl import NCMATComposerImpl as Impl
        assert _default_debye_temp()==400.0, "from_hfg default must be updated"
        ncmat = hfg2ncmat( spec = spec,
                           formula = formula,
                           density = density,
                           title = title,
                           debyetemp = debyetemp,
                           verbose = verbose,
                           notrim = notrim )
        c = NCMATComposer( Impl.from_ncmat( ncmat,
                                            keep_header=True ) )
        c.set_plotlabel(title)
        return c

    def __call__( self, cfg_params = None, **kwargs ):
        """Convenience short-cut for the .create_ncmat() method"""
        return self.create_ncmat( cfg_params = cfg_params, **kwargs )

    def update_atomdb( self, element_or_isotope, data = None, *,  mass = None, coh_scat_len=None, incoh_xs=None, abs_xs=None ):
        """Specify atom data for particular element or isotope. Either specify
        the data parameter as a string (in @ATOMDB format like "28.97649466525u
        4.7fm 0.001b 0.101b") or AtomData object, or put all four physics values
        (mass, coh_scat_len, incoh_xs, and abs_xs) directly as floating point
        numbers. More specifically these numbers are mass value (amu), coherent
        scattering length (fm), incoherent cross section (barn), and absorption
        cross section @v_n=2200m/s (barn).
        """
        return self.__impl.update_atomdb( element_or_isotope,data=data,mass=mass,
                                          coh_scat_len=coh_scat_len,incoh_xs=incoh_xs,abs_xs=abs_xs )

    def find_label( self, element, allow_multi=False ):
        """Attempt to locate label which has a composition containing provided
        element (which can be a string name or Z value like Al, 13, B10,
        ...). If exactly one such label exists, it will be returned, otherwise
        None is returned. If allow_multi==True, the result(s) will instead be
        returned as a list, and all matching labels will be returned.
        """
        return self.__impl.find_label( element = element, allow_multi = allow_multi )

    def set_composition( self, label, *composition ):
        """Use like .set_composition('H','D') or .set_composition('H is D'). Mixtures
        can also be defined, for instance like .set_composition('H','0.2 H1 0.8
        H2'). For a more in-depth discussion of the syntax, consult the @ATOMDB
        section on https://github.com/mctools/ncrystal/wiki/NCMAT-format

        Set atomic composition associated with a given label. Usage examples:

        .set_composition('mylbl','Al') #mylbl now has a composition of only (natural) Al atoms.
        .set_composition('Al','Al') #Al atoms (notice this call is not needed
                                    #since the label is already an element name).
        .set_composition('Al','Cr') #Well if you want to use confusing labels you can
        .set_composition('mylbl','Li6') #specific isotope
        .set_composition('mylbl','D') #deuterium is supported by either 'D' or 'H2'

        More complicated mixtures are supported as well:

        .set_composition('mylbl',0.95,'B10',0.05,'B11')#enriched boron
        .set_composition('mylbl',[0.95,'B10',0.05,'B11'])#same but defined via sequence
        .set_composition('mylbl',[(0.95,'B10'),(0.05,'B11')])#same but defined via sequence of pairs
        .set_composition('mylbl','0.95 B10 0.05 B11')#same but via a string in @ATOMDB-like form.
        .set_composition('mylbl is 0.95 B10 0.05 B11')#same but even more @ATOMDB like.
        .set_composition('Al is 0.99 Al 0.005 Cr 0.005 B10')#mix and match elements and isotopes
        """
        return self.__impl.set_composition( label, *composition )

    def remap_atom( self,  element_or_isotope, *composition ):
        """Calling this method updates all usage of the denoted
        element_or_isotope marker ('C', 'Al', 'D', 'B10', ...) and replaces it
        with the particular composition. For instance, .remap_atom('H','D') can
        be used to fully deuterate a material. For a description of the allowed
        compositon syntax, refer to the .set_composition method.

        The remapping will only affect labels and compositions that have already
        been added to the NCMATComposer object.
        """
        return self.__impl.remap_atom( element_or_isotope, *composition )

    def clear_comments( self ):
        """Clear comments added with the add_comments method."""
        return self.__impl.clear_comments()

    def add_comments( self, comments, add_empty_line_divider=False ):
        """Add comments (string or list of strings) to the top of NCMAT data."""
        return self.__impl.add_comments( comments = comments,
                                         add_empty_line_divider = add_empty_line_divider )

    def clone(self):
        """Return a new independent NCMATComposer object which is a clone of the
        current one."""
        return NCMATComposer( self.__impl.clone() )

    def get_cache_key( self ):
        """
        Returns tuple of integers which will be unchanged if generated NCMAT
        content is unchanged. Note that this is partially based on memory
        addresses of python objects, so the keys are NOT reproducible outside
        the current process.
        """
        return self.__impl.get_cache_key()

    def to_dict(self):
        return self.__impl.to_dict()

    def set_cellsg( self, *, a,b,c, alpha,beta,gamma, spacegroup=None ):
        """Assumes alpha,beta,gamma in degrees, a,b,c in Aangstrom, spacegroup
        an int in range 1..230."""
        return self.__impl.set_cellsg( a=a,b=b,c=c,
                                       alpha=alpha,beta=beta,gamma=gamma,
                                       spacegroup=spacegroup )

    def set_cellsg_cubic( self, a, *, spacegroup=None ):
        """Like .set_cellsg but will set b=c=a and alpha=beta=gamma=90."""
        return self.__impl.set_cellsg_cubic( a = a, spacegroup = spacegroup )

    def set_atompos( self, atompos ):
        """The atompos parameter must be list of entries
        (label,x,y,z) or (label,x,y,z,site_occupancy).

        If provided, the site_occupancy value (a number in (0,1] must be the
        same for all identical labels.

        The site_occupancy parameter is for now highly experimental and will
        result in generated NCMAT files with fake sterile atoms inserted for
        technical reasons (resulting in a higher density but hopefully correct
        physics as long as the user does not subsequently try to override the
        density directly via the cfg-string parameters). This is perfectly fine
        for usage in a MC simulation context where NCrystal provides the entire
        physics of the material (e.g. McStas), but in a context where a base
        material must be created using the composition from NCrystal
        (e.g. Geant4 or OpenMC), it most likely will result in an error.
        """
        return self.__impl.set_atompos( atompos = atompos )

    def refine_crystal_structure( self, symprec = 0.01, quiet = False ):
        """Attempt to refine the crystal structure (does nothing if not a
        crystalline material). This ignores the spacegroup number (if any), and
        uses spglib to refine and standardise the crystal structure based on the
        provided unit cell parameters and atomic positions. The unit cell
        parameters, spacegroup number, and atomic positions, are all updated
        accordingly.
        """
        return self.__impl.refine_crystal_structure( symprec = symprec,
                                                     quiet = quiet )

    def verify_crystal_structure( self, symprec = 0.01, quiet = False ):
        """Attempt to verify the crystal structure (does nothing if not a
        crystalline material). This requires a spacegroup number to be
        available, and uses spglib to verify that the crystal structure (unit
        cell and atom positions) is consistent with that spacegroup
        number. Raises an exception in case the structure could not be verified,
        otherwise does nothing.
        """
        return self.__impl.verify_crystal_structure( symprec = symprec,
                                                     quiet = quiet )

    def set_density( self, value, unit = 'g/cm3' ):
        """
        Specify material density for non-crystalline materials. If provided,
        the unit must be a string which is one of "g/cm3" (default) , "kg/m3",
        or "atoms/Aa3". Note that any material can have its density modified
        with the cfg-string parameter "density", while the present method is
        instead used to provide the basic density embedded in NCMAT format for
        non-crystalline materials.
        """
        return self.__impl.set_density( value = value, unit = unit )

    def set_fraction( self, label, value ):
        """
        Specify fraction of component associated with a given label. This is
        not needed for crystalline materials, where the fractions can be
        inferred from the unit cell contents. Note that the various
        set_dyninfo_xxx methods have an optional fraction parameter, which can
        also be used to set the fraction directly.

        """
        return self.__impl.set_fraction( label = label, value = value )

    def allow_fallback_dyninfo( self, debye_temp = 300.0 ):
        """For crystalline materials only, allow a fall-back modelling of any
        components with absent dynamic information via a VDOS-Debye model using
        the specified Debye temperature.
        """
        return self.__impl.allow_fallback_dyninfo( debye_temp = debye_temp )

    def set_dyninfo_vdos( self, label, vdos_egrid, vdos, *, fraction = None, comment = None ):
        """Set dynamics of component to be modelled by a 1D phonon density of
        state (DOS) curve (VDOS="Vibrational DOS"). This is the preferred way to
        provide dynamics for components of solid materials, as it not only
        gives enough information to estimate the Debye-Waller factors needed for
        elastic scattering in all such materials, but also allows realistic
        inelastic scattering via expansion into full 2D scattering kernels.
        This model is not appropriate for gaseous or liquid materials.
        """
        return self.__impl.set_dyninfo_vdos( label = label,
                                             vdos_egrid = vdos_egrid,
                                             vdos = vdos,
                                             comment = comment,
                                             fraction = fraction )

    def set_dyninfo_vdosdebye( self, label, debye_temp, *, fraction = None, comment = None ):
        """Set dynamics of component to be modelled by a Debye temperature. This
        allows temperature-dependent atomic displacements and Debye-Waller
        factors to be estimated, and inelastic scattering will be modelled based
        on an idealised VDOS curve (a parabola up to a cutoff frequency
        corresponding to E=k*T_debye). This model is only appropriate for solid
        materials.
        """
        return self.__impl.set_dyninfo_vdosdebye( label = label,
                                                  debye_temp = debye_temp,
                                                  comment = comment,
                                                  fraction = fraction )

    def set_dyninfo_debyetemp( self, label, debye_temp, *, fraction = None, comment = None ):
        """Alias for set_dyninfo_vdosdebye."""
        return self.set_dyninfo_vdosdebye( label=label, debye_temp=debye_temp,
                                           fraction=fraction, comment=comment )

    def set_dyninfo_msd( self, label, msd, *, temperature, fraction = None, comment = None,  ):
        """Calculate and set Debye temperature based on
        mean-squared-displacement (msd) value (in Aa^2). This also needs the
        temperature value for which the msd value is associated (in
        kelvin). Note that this temperature value is NOT necessarily the
        temperature at which the material will later be used. See
        .set_dyninfo_vdosdebye(..) for further information about the resulting
        modelling.
        """
        return self.__impl.set_dyninfo_msd( label = label,
                                            msd = msd,
                                            temperature = temperature,
                                            comment = comment,
                                            fraction = fraction )

    def set_dyninfo_uiso( self, label, uiso, temperature, *, fraction = None, comment = None ):
        """Alias for set_dyninfo_msd. The uiso value is just another name for
        "msd" and both have a unit of Aa^2."""
        return self.set_dyninfo_msd( label=label, msd=uiso, temperature=temperature,
                                     fraction=fraction, comment=comment )

    def set_dyninfo_scatknl( self, label, *, alphagrid, betagrid, temperature,
                             sab = None, sab_scaled = None, egrid = None,
                             fraction = None, comment = None ):
        """Set dynamics of component to be modelled by the provided 2D
        S(alpha,beta) ("sab") scattering kernel. Such kernels are
        temperature-dependent, so one must also specify the temperature for
        which the kernel is valid, which of necessity will lock the temperature
        of the entire material at that temperature. As NCrystal currently is not
        able to estimate atomic displacements and Debye-Waller factors from such
        kernels, they can not be used with crystalline materials at all, and are
        in principle not suitable for amorphous solids either, as even
        incoherent-elastic contributions will be absent. Thus, this modelling is
        mainly intended for liquids or gasses.
        """
        return self.__impl.set_dyninfo_scatknl( label = label,
                                                alphagrid = alphagrid,
                                                betagrid = betagrid,
                                                temperature =  temperature,
                                                sab = sab,
                                                sab_scaled = sab_scaled,
                                                egrid = egrid,
                                                comment = comment,
                                                fraction = fraction )

    def set_dyninfo_freegas( self, label, *, fraction = None, comment = None ):
        """Set dynamics of component to be modelled by a free gas description (a gas
        of non-interacting atoms). This model is not suitable for solid
        materials, except as a crude approximation. It is specifically
        disallowed for crystalline materials, since it does not provide atomic
        displacements (Debye-Waller factors). See .set_dyninfo_debyetemp(..),
        .set_dyninfo_msd(..) and .allow_fallback_dyninfo() for easy
        (low-realism) alternatives for crystalline materials.
        """
        return self.__impl.set_dyninfo_freegas( label = label,
                                                comment = comment,
                                                fraction = fraction )

    def set_dyninfo_sterile( self, label, *, fraction = None, comment = None ):
        """Set dynamics of component to be modelled as being without any scattering
        cross section. This model is obviously unrealistic and exists for
        debugging purposes only. It can not be used for crystalline materials
        (for these one could use .update_atomdb instead to introduce atoms with
        no scattering cross section into the crystal).
        """
        return self.__impl.set_dyninfo_sterile( label = label,
                                                comment = comment,
                                                fraction = fraction )

    def set_dyninfo_from_object( self, label, source_dyninfo, comment = None, fraction=None ):
        """Set dyninfo of given label based on existing NCrystal.DynamicInfo
        object.

        Note that in principle the flow NCMAT-data -> Info -> NCMATComposer ->
        NCMAT-data is not guaranteed to be 100% lossless, but in practice the
        resulting physics is unlikely to be unaffected by this. For a guaranteed
        loss-less operation, one must manually edit the NCMAT data and copy the
        exact lines of the relevant @DYNINFO section into the new NCMAT data
        (this is beyond the scope of the NCMATComposer).
        """
        return self.__impl.set_dyninfo_from_object( label = label,
                                                    source_dyninfo = source_dyninfo,
                                                    comment = comment,
                                                    fraction = fraction )

    def transfer_dyninfo_objects( self, source, mapping = None, allow_none = False ):
        """Set dyninfo from source. Source can be a list of
        NCrystal.DynamicInfo, a single NCrystal.DynamicInfo info object, or any
        sort of material source like an Info object, a cfg-str, NCMAT data, etc.

        If a mapping dict is provided, it must be a mapping between labels used
        on the NCMATComposer object and display labels in the source DynamicInfo
        objects (example: {"H":"H","mylbl":"Al","O16":"O"}). If a mapping is NOT
        provided, an automatic mapping will be attempted based on element Z
        values between source entries and labels already registered in the
        composer (in this case it must be possible to determine the composition
        of all labels). Any excess entries in the source will be ignored.

        An NCBadInput error will be raised if no dynamic information is actually
        transferred, unless allow_none=True is provided.

        """
        return self.__impl.transfer_dyninfo_objects( source = source,
                                                     mapping = mapping,
                                                     allow_none = allow_none )

    def set_state_of_matter( self, state_of_matter ):
        """Explicitly specify the state of matter, which must be a string with
        value of "solid", "liquid", or "gas". Call it with an empty string or
        None to clear the information.

        Note that NCrystal will automatically classify materials as "solid" if
        they are crystalline materials and/or have dynamic information of the
        types vdos, vdosdebye, or msd. Trying to change the state of matter
        value for such materials will result in an error at a later stage.
        """
        return self.__impl.set_state_of_matter( state_of_matter )


    def add_secondary_phase(self, fraction, cfgstr, normalise = True ):
        """Add a secondary phase to the material, by specifying the phase volume
        fraction and a cfg-strings defining the content of the secondary
        phase. For instance a material with air bubbles inside, might be
        emulated by calling .add_secondary_phase(0.05,"gasmix::air"), which
        would add a secondary phase to the material, occupying 5% of the volume
        and containing air. For details, please see the discussion under the
        heading "The @OTHERPHASES section" in the document at
        https://github.com/mctools/ncrystal/wiki/NCMAT-format.

        Unless normalise=False, the cfg-string will be normalised with a call to
        normaliseCfg(..) from the NCrystal.cfgstr module.
        """
        return self.__impl.add_secondary_phase( fraction = fraction,
                                                cfgstr = cfgstr,
                                                normalise = normalise )

    def add_hard_sphere_sans_model( self, sphere_radius ):
        """As a technology preview of future more complete SANS support, one can
        for now use this method to add hard-sphere SANS scattering between the
        primary phase, and the first phase added with .add_secondary_phase(..)
        Note that as this results in a @CUSTOM_HARDSPHERESANS section in the
        resulting NCMAT data (because it is a tech-preview), warnings will be
        emitted when the data is subsequently loaded by NCrystal.

        The implemented hard-sphere model only has a single free parameter, the
        sphere_radius in angstrom.
        """
        return self.__impl.add_hard_sphere_sans_model(sphere_radius=sphere_radius)

    def add_raw_content( self,  content ):
        """This is an experts-only method for adding raw text data to be
        appended to the generated NCMAT data. Note that multiple calls to this
        method will simply append more content, if you wish to remove content
        again you must call clear_raw_content.

        Note that this method does not necessarily add a newline to your content
        if it is missing.

        """
        return self.__impl.add_raw_content( content )

    def clear_raw_content( self ):
        """Remove any content added by calls to .add_raw_content(..)"""
        return self.__impl.clear_raw_content()

    def get_raw_content( self ):
        """Return any content added by calls to .add_raw_content(..). Returns an
        empty string if no such content was added."""
        return self.__impl.get_raw_content()

    def set_custom_section_data( self, section_name, content ):
        """Add a @CUSTOM_<sectionname> section with the provided content to the
        generated NCMAT data. Note that multiple calls to this method with the
        same section_name will simply override the content of that custom
        section. To completely remove a previously added custom section, use the
        .clear_custom_section_data(..) method.
        """
        return self.__impl.set_custom_section_data( section_name, content )

    def get_custom_section_data( self, section_name = None ):
        """Access any @CUSTOM_<sectionname> contents which was previously added
        with the .set_custom_section_data(..) method. Returns None if data for
        that section was not added. If called without parameters, a dictionary
        of all such data in the form { section_name : content, ... } is
        returned instead.
        """
        return self.__impl.get_custom_section_data( section_name )

    def clear_custom_section_data( self, section_name = None ):
        """Remove any @CUSTOM_<sectionname> section previously added by calling
        set_custom_section_data. Does nothing if no such section was previously
        added. Calling with no arguments clears all custom section data.
        """
        return self.__impl.clear_custom_section_data( section_name )

    def lock_temperature( self, value ):
        """
        Lock the temperature of the material to the given value. This not
        only changes the default temperature of the material, but also "locks"
        it in the sense that any attempts at using the cfg-level variable "temp"
        to modify the temperature will result in an error. See also the
        .set_default_temperature() method. Call with value=None to clear any
        previous effects of .lock_temperature() or .set_default_temperature().
        """
        return self.__impl.lock_temperature( value )

    def set_default_temperature( self, value ):
        """
        Modify the default temperature of the material from the usual
        293.15K. The temperature value can still be modified subsequently using
        the cfg parameter "temp" as usual. To prevent this, use the
        .lock_temperature() method instead. Call with value=None to clear any
        previous effects of .lock_temperature() or .set_default_temperature().
        """
        return self.__impl.set_default_temperature( value )

    def get_temperature_setting( self ):
        """Returns the temperature setting resulting from calls to
        .lock_temperature() or .set_default_temperature(). Returns a tuple of
        two values: the temperature value in kelvin, and a boolean indicating
        whether the value has been locked or is merely a changed
        default. Returns (None,None) in the absence of such a setting.
        """
        return self.__impl.get_temperature_setting()

    @property
    def state_of_matter( self ):
        """
        Access the state of matter which will be one of "solid", "liquid",
        or "gas". However, this is only available if explicitly set, otherwise
        None is returned.
        """
        return self.__impl.get_state_of_matter()

    def get_labels( self ):
        """List of all labels registered so far."""
        return self.__impl.get_labels()

    def write( self, path, cfg_params = None ):
        """Produce NCMAT data and write it to the provided file path. If any
        cfg_params are provided, they will be embedded in the NCMAT data using
        the NCRYSTALMATCFG[..] syntax.
        """
        return self.__impl.write( path = path,
                                  cfg_params = cfg_params )

    def register_as( self, virtual_filename, cfg_params = None ):
        """Produce NCMAT data and register it in memory with the provided
        virtual file-name, using the registerInMemoryFileData(..) function from
        the NCrystal.datasrc module. If any cfg_params are provided, they will
        be embedded in the NCMAT data using the NCRYSTALMATCFG[..] syntax.
        """
        return self.__impl.register_as( virtual_filename = virtual_filename,
                                        cfg_params = cfg_params )

    def load(self, cfg_params = None, *, force = False ):
        """Will create NCMAT and load it with the directLoad(..) function of the
        NCrystal.core module. The cfg_params parameter can be used to apply cfg
        parameters (e.g. cfg_params="temp=200K;dcutoff=0.1")., while force=True
        will prevent the NCMATComposer from simply returning a previously loaded
        material (usually one should just leave the default force=False value
        untouched).
        """
        return self.__impl.load( cfg_params = cfg_params, force = force )

    def plot_xsect( self, cfg_params = None, **kwargs_plot_xsect ):
        """
        Quick plot (with matplotlib) showing the cross sections produced by
        current material, possibly after appending certain cfg_params. This is
        using the plot_xsect function from the NCrystal.plot module, so refer to
        that function for available arguments.
        """
        return self.__impl.plot_xsect( self, cfg_params, kwargs_plot_xsect )

    def inspect( self, cfg_params = None, **kwargs_plot_xsect ):
        """Similar to the plot_xsect(..) function, but will also print out
        information (e.g. .dump()) about the Info, Scatter, and Absorption
        objects.
        """
        return self.__impl.inspect( self, cfg_params, kwargs_plot_xsect )

    def get_chemical_composition( self, as_str = False ):
        """Chemical composition as list of [(elemiso,count),...] where elemiso
        is a marker like "Al", "D", "B10", .... For crystals, the counts will
        represent the number of elements/isotopes per unit cell. Returns None in
        case of incomplete information.

        If as_str = True, the result will instead be returned encoded in a
        string (e.g. "H2O").
        """
        return self.__impl.get_chemical_composition( as_str=as_str )

    def create_ncmat( self, cfg_params = None, *,
                      meta_data = False,
                      verify_crystal_structure = True ):
        """Creates and returns NCMAT data based on the material settings
        provided so far. Any cfg_params provided here will be embedded into the
        NCMAT data itself, using the NCRYSTALMATCFG[..] syntax.

        If meta_data=True, the return value will be a tuple whose first element
        is the created NCMAT data, and whose second element is a dictionary with
        a high-level meta-data concerning the material (like chemical
        composition, spacegroup, etc.).

        Unless verify_crystal_structure=False (not recommended!) any crystalline
        material will have its structure verified by calling
        .verify_crystal_structure(), possibly raising an exception in case of
        issues found. Depending on the scenario, it might be possible to prevent
        this by calling .refine_crystal_structure(), which will hopefully fix it
        up before verification (but be vigilant and double-check that the
        refinement did not simply hide some fundamental flaw of the input data).
        """
        return self.__impl.create_ncmat( cfg_params = cfg_params,
                                         meta_data = meta_data,
                                         verify_crystal_structure = verify_crystal_structure )


    def set_plotlabel( self, lbl ):
        """Sets the plotlabel (cf. .plotlabel)."""
        return self.__impl.set_plotlabel(lbl)

    @property
    def plotlabel( self ):
        """An optional label for the material, which can be used for plotting
        purposes in legends, etc.
        """
        return self.__impl.plotlabel

    def as_spglib_cell( self ):
        """For a crystalline material, return the unit cell definition in a
        format suitable for usage in spglib calls."""
        return self.__impl.as_spglib_cell()

    def _unofficial_vdos2sab_ignore( self, *, order_low, order_high = None, mode = None ):
        """This expert-only method can be used to exclude certain phonon orders
        from the scattering kernel during vdos2sab expansion. It is primarily
        intended as a debugging aid for visualisation purposes, or for
        developers intending to develop coherent single-phonon physics.

        Be aware that, for now, trying to use the Info.DI_VDOS[Debye].loadKernel
        or .plot_knl(), will NOT reflect this hack!
        """
        return self.__impl._unofficial_vdos2sab_ignore( order_low = order_low,
                                                        order_high = order_high,
                                                        mode = mode )

    def _get_impl_obj( self ):
        return self.__impl

def formatVectorForNCMAT(name,values,indent='  '):
    """Utility function for help in python scripts composing .ncmat files,
       transforming an array of of values into a properly formatted text string,
       with word-wrapping, usage of <val>r<n> syntax, etc.
    """
    from ._ncmatimpl import formatVectorForNCMAT as _
    return _( name = name,
              values = values,
              indent = indent )

def _rawParseNCMAT(text_data_name,*,asJSONStr=False):
    """Parses NCMAT content and returns as Python data structure (a dictionary). The
       format of this data structure should be mostly self-evident by
       inspection, and is not guaranteed to stay the same across NCrystal
       versions. If asJSONStr=true, the data structure will be returned as a
       JSON-encoded string, instead of a Python dictionary.

       WARNING: This function is considered experimental and is currently NOT
       feature complete. It only returns data from a few select NCMAT sections."""
    from ._chooks import _get_raw_cfcts
    _js = _get_raw_cfcts()['nc_ncmat2json'](text_data_name if not hasattr(text_data_name,'rawData') else text_data_name.rawData )
    if asJSONStr:
        return _js
    import json
    return json.loads(_js)
