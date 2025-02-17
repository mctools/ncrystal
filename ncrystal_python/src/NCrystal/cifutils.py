
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

"""Utilities for converting local or online CIF data to NCMAT data.

Two utility classes (CIFSource and CIFLoader) are included, and along with the
NCMATComposer they are used to implement the related cmd-line scripts (in
particular ncrystal_cif2ncmat), but they can of course also be used directly via
the Python API below by expert users.

"""

__all__ = ['CIFSource','CIFLoader',
           'produce_validation_plot',
           'produce_validation_plots']

from . import _common as _nc_common
from . import ncmat as _nc_ncmat
from . import _ncmatimpl as _nc_ncmatimpl
from . import core as _nc_core

class CIFSource:

    """
    Generalised CIF source, either in the form of text data, a file path, a
    database ID (in the form of a string like "codid::9008460" or "mpid::87)",
    for either the Crystallography Open Database
    (https://www.crystallography.net/cod/) or the Materials Project
    (https://www.materialsproject.org/), an existing CIFSource object, or any
    object like CIFLoader which has a cifsrc property which is a CIFSource
    object.
    """

    def __init__( self, data, *, allow_fail = False, name = None ):
        """Initialise from data in various formats (see class
        description). Unless allow_fail=True, an unrecognised input format will
        result in an NCBadInput exception being thrown. If allow_fail=True the
        .invalid property can be used to check if the loading failed.

        The name argument can be used to assign a name to otherwise anonymous
        text data. It can also be used to override the name.

        """
        self.__codid = None
        self.__mpid = None
        self.__fp = None
        self.__textdata = None
        self.__name_override = name

        def check_status():
            if ( not allow_fail ) and self.invalid:
                raise _nc_core.NCBadInput('Could not detect CIF source type')

        #Check if data is CIFSource object or has one as a .cifsrc property (e.g. LoadedCIF, CIFAnalyser)
        o = getattr( data, 'cifsrc', None ) or data
        if isinstance(o,CIFSource):
            self.__codid = o.__codid
            self.__mpid = o.__mpid
            self.__fp = o.__fp
            self.__textdata = o.__textdata
            self.__name_override = o.__name_override
            return
        def _setfp( pth ):
            pth = pth.decode() if hasattr(pth,'decode') else pth
            self.__fp = _nc_common._lookup_existing_file( pth )
            if self.__fp is None:
                #Try to look up via NCrystal's TextData infrastructure:
                from .core import createTextData, NCFileNotFound
                try:
                    td = createTextData(pth)
                    tdname = td.dataSourceName
                except NCFileNotFound:
                    td = None
                    tdname = None
                if td is not None:
                    self.__textdata = td.rawData
                    if tdname and not self.__name_override:
                        self.__name_override = tdname
                    self.__fp = None
            if not allow_fail and self.__fp is None:
                raise _nc_core.NCBadInput('Could not detect CIF source type'
                                          ' (tried to load from path, but'
                                          ' could not locate "%s")'%pth)

        if hasattr(data,'__fspath__'):
            _setfp(data)
            return check_status()
        if hasattr( data, 'startswith' ):
            str_to_sb = ( lambda s : s ) if not hasattr(data,'decode') else ( lambda s : s.encode() )

            if data.startswith(str_to_sb('codid::')) and data[7:].isdigit():
                self.__codid = int(data[7:])
                return check_status()
            if data.startswith(str_to_sb('mpid::')) and data[6:].isdigit():
                self.__mpid = int(data[6:])
                return check_status()
        if isinstance( data, str ):
            if '\n' in data.strip():
                self.__textdata = data
            else:
                _setfp(data)
            return check_status()
        if isinstance( data, bytes ):
            if b'\n' in data.strip():
                self.__textdata = data.decode()
            else:
                _setfp(data)
            return check_status()
        check_status()

    @property
    def invalid( self ):
        """Check if instance is invalid (only possible if constructed with
        allow_fail=True). """
        return all(e is None for e in (self.__codid,self.__mpid,
                                       self.__fp,self.__textdata))

    @property
    def codid( self ):
        """
        None or a database ID (integer) for the Crystallography Open Database
        (https://www.crystallography.net/cod/).
        """
        return self.__codid

    @property
    def mpid( self ):
        """
        None or a database ID (integer) for the Materials Project
        (https://www.materialsproject.org/).
        """
        return self.__mpid

    @property
    def textdata( self ):
        """
        None or a string containing raw CIF data.
        """
        return self.__textdata

    @property
    def name( self ):
        """A short string describing the data. This might be a file name
        (without directory part), or a special string like 'mpid::xyz'
        'codid::xyz'. Returns None if no name is available.
        """
        if self.__name_override:
            return self.__name_override
        if self.__mpid:
            return 'mpid::%i'%self.__mpid
        if self.__codid:
            return 'codid::%i'%self.__codid
        if self.__fp:
            return self.__fp.name

    @property
    def filepath( self ):
        """
        None or a pathlib.Path to an on-disk CIF file.
        """
        return self.__fp

    @property
    def is_remote( self):
        """True if either .codid or .mpid is not None"""
        return any( ( e is not None ) for e in (self.__codid,self.__mpid) )

    def load_data( self, quiet = False, mp_apikey = None ):
        """Try to load the data and return it as a string. Depending on the
        source, this might result in reading a file or querying an online
        database. In case the data is defined via a materials project ID (mpid),
        a materials project API key must be provided either via the mp_apikey
        parameter, or in the MATERIALSPROJECT_USER_API_KEY environment variable.
        """
        if self.invalid:
            raise _nc_core.NCDataLoadError( 'Can not load_data'
                                            ' from invalid CIF source' )
        if self.__textdata is not None:
            return self.__textdata
        if self.__fp is not None:
            if not quiet:
                _nc_common.print(f"Loading data from file {self.__fp.name}")
            return self.__fp.read_text()
        if self.__mpid is not None:
            return _mp_get_cifdata( self.__mpid, quiet = quiet, apikey = mp_apikey )
        assert self.__codid
        return _cod_get_cifdata( self.__codid, quiet = quiet )

class CIFLoader:

    """Class which is used to load and analyse CIF data, primarily so it can be
    used to transfer the crystal structure into an NCMATComposer object, with
    the intention of producing NCMAT data and new materials for usage with
    NCrystal.

    Behind the scenes the class uses third-party projects gemmi and spglib to
    parse the CIF data and process the symmetries within.

    The perhaps biggest issue to understand when basing NCrystal materials on
    CIF data, is that usually CIF data does not contain any or only incomplete
    information about the dynamics of the atoms. Thus, for a complete
    description of a high-quality material, one might want to pair up the usage
    of a CIFLoader, with a PhononDOSAnalyser from the NCrystal.vdos
    module.

    However, some CIF data contain information about atomic displacements
    (either "Uiso" or the anisotropic equivalent), which NCrystal is able to
    convert into a Debye temperature in order to provide a temperature-dependent
    model of dynamics and atomic displacement. However, such "Uiso" information
    is only useful if one knows at which temperature it was obtained, and that
    information is almost always absent in most CIF data out there. Thus, in the
    last step - when using the .create_ncmat(..) or .create_ncmat_composer(..)
    methods of the CIFLoader object, one should ideally use the uiso_temperature
    parameter to provide this information to NCrystal (usually after a call to
    print(cifloader.raw_cif_textdata) to manually inspect the meta-data, read
    any linked journal papers, etc.).

    TODO: Add link here to online tutorial once it is ready.

    """

    def __init__( self, cifsrc, quiet = False, mp_apikey = None,
                  refine_with_spglib = True, merge_equiv = True, override_spacegroup = None ):
        """Initialise from cifsrc which can either be a CIFSource object, or any
        sort of data which can be used to initialise a CIFSource object (cf. the
        CIFSource class documentation)

        If refine_with_spglib=False, the structure will not be verified and
        refined with spglib. It is strongly recommended to not disable this
        refinement, as it is needed to ensure correctness and consistency of the
        resulting material.

        The override_spacegroup parameter can in some cases be used to select a
        particular setting of the detected spacegroup. For instance, spacegroup
        number 70 is available in two settings, and loading a CIF file with
        those might result in a warning that both "F d d d:1" and "F d d d:2"
        are available. Trying again with override_spacegroup="F d d d:1" can be
        used to explicitly select one of them (it might require some
        investigation of course to determine which one is right). A list of
        space group settings can be found at:

        https://cci.lbl.gov/sginfo/hall_symbols.html

        Setting merge_equiv=False will prevent usage of the same NCMAT label for
        multiple lists of equivalent atomic positions. However, note that this
        might affect the symmetry and space group of the structure, and is in
        general not recommended.

        If quiet=True, no informative messages will be emitted.

        Finally, the mp_apikey can be used to specify an API key in case the
        cifsrc is an entry for the materialsproject.org. To avoid having to
        specify the key here, one can instead provide the key in the environment
        variable MATERIALSPROJECT_USER_API_KEY.
        """
        with _nc_common.WarningSpy( block = quiet ) as warnlist:
            load, gemmi, struct = _actual_init_gemmicif( cifsrc,
                                                         quiet = quiet,
                                                         mp_apikey = mp_apikey,
                                                         refine_with_spglib = refine_with_spglib,
                                                         merge_equiv = merge_equiv,
                                                         override_spacegroup = override_spacegroup )
        assert len(load)==9
        self.__cifsrc = load['cifsrc']
        self.__actual_codid = load['actual_codid']
        self.__actual_mpid = load['actual_mpid']
        self.__cif_chemformula = load['cif_chemformula']
        self.__cifdata = load['cifdata']
        self.__cif_raw = load['cif_raw']
        self.__cellsg = load['cellsg']
        self.__atoms = load['atoms']
        self.__cifdescr = load['cifdescr']
        self.__warnlist = tuple( warnlist )

        #Hidden support for additional processing in derived class:
        _ = getattr( self, '_process_raw_gemmi_objects', None )
        if _:
            _( gemmi, struct )

    @property
    def cifsrc( self ):
        """The CIFSource object on which everything is based."""
        return self.__cifsrc

    @property
    def actual_codid( self ):
        """None or the database ID of the Crystallography Open Database. This
        might be present even if .cifsrc.codid is None, since the CIF data might
        have been downloaded manually, but provide the ID in its metadata."""
        return self.__actual_codid

    @property
    def actual_mpid( self ):
        """None or the database ID of the Material Project. This might be
        present even if .cifsrc.mpid is None, since the CIF data might have
        been downloaded manually, but provide the ID in its metadata.
        """
        return self.__actual_mpid

    @property
    def raw_cif_textdata( self ):
        """The raw CIF input data, as a string."""
        return self.__cifdata

    @property
    def raw_cif_chemformula( self ):
        """The chemical formula encoded in the CIF data (usually in the
        _chemical_formula_sum field) as a string."""
        return self.__cif_chemformula

    @property
    def raw_cif_dict( self ):
        """The CIF data parsed into a dictionary."""
        return self.__cif_raw

    @property
    def atoms( self ):
        """A tuple of all atoms in the crystal structure, including positions in
        the unit cell. The format of each atom is essentially a dictionary
        similar to:
               {'aniso': None,
                'cif_labels': ['Al'],
                'composition': ((1.0, 'Al'),),
                'equivalent_positions': [(0.0, 0.0, 0.0),
                                         (0.0, 0.5, 0.5),
                                         (0.5, 0.0, 0.5),
                                         (0.5, 0.5, 0.0)],
                'occupancy': 1.0,

                'uiso': None})

        Here 'aniso' might provide information about
        anisotropic displacements, 'uiso' (Uiso) about isotropic displacements,
        and occupancy is the site occupancy. The 'cif_labels' are the
        corresponding labels in the CIF data, and the 'composition' +
        'equivalent_positions' should be self-explanatory.
        """
        return self.__atoms

    @property
    def cellsg( self ):
        """
        Get information about unit cell parameters and space group, as a
        dictionary with a format like:

           {'a': 4.0496, 'b': 4.0496,'c': 4.0496,
            'alpha': 90,'beta': 90,'gamma': 90,
            'spacegroup': {'hm': 'Fm-3m', 'number': 225} }

        """
        return self.__cellsg

    @property
    def extracted_description( self ):
        """
        Description of the material based on the meta-data found in the CIF
        data. Returned as a list of strings.
        """
        return self.__cifdescr

    @property
    def warnings( self ):
        """
        List of warnings emitted during processing of the CIF data, as a
        sequence of tuples of two strings: (warning_message, warning_category).
        """
        return self.__warnlist

    def create_ncmat( self, *args, **kwargs ):
        """
        Return NCMAT data corresponding to the currently held crystal
        parameters. All parameters are forwarded (same as calling
        create_ncmat_composer(*args,**kwargs).create_ncmat())
        """
        c = self.create_ncmat_composer(*args,**kwargs)
        return c.create_ncmat()

    def create_ncmat_composer( self, *,
                               uiso_temperature = None,
                               skip_dyninfo = False,
                               top_comments = None,
                               fallback_debye_temp = 300.0,
                               quiet = False,
                               remap = None,
                               no_formula_check = False ):
        """Setup and return an NCMATComposer object, based on the currently held
        crystal parameters.

        The most important thing to be aware of is that, if you wish to take
        advantage of any Uiso/biso values found in the CIF data, you must
        provide the uiso_temperature parameter value (kelvin), to indicate the
        temperature at which these are valid (most likely after inspecting the
        .raw_cif_textdata manually for hints about the correct value).

        If quiet=True, no informative messages will be emitted.

        The fallback_debye_temp parameter can be used to change the default
        Debye temperature assigned to atoms, for which Uiso/biso information was
        not available.

        The other parameters are more esoteric:

        If skip_dyninfo=True, no ncmat.set_dyninfo_... calls will be performed,
        presumably because the caller will perform such calls subsequently
        (i.e. with a PhononDOSAnalyser).

        The remap parameter can be used to remap elements and isotopes found in
        the CIF data to any desired composition (cf. the
        NCMATComposer.set_composition for the allowed composition syntax. For
        instance, if one wishes to take a CIF file containing hydrogen atoms and
        deuterate it partly, one could used remap = { 'H' : '0.95 D 0.05 H', }.

        Finally, unless no_formula_check=True, an exception is raised if the
        chemical formula of the resulting material is incompatible with any
        chemical formula indicated in the input CIF meta data. This check is
        normally nice to keep, since it can catch many cases where the wrong
        setting of a spacegroup was used.
        """
        return _impl_create_ncmat_composer( self,
                                            uiso_temperature = uiso_temperature,
                                            skip_dyninfo = skip_dyninfo,
                                            top_comments = top_comments,
                                            fallback_debye_temp = fallback_debye_temp,
                                            #TODO? #merge_equiv = merge_equiv,
                                            quiet = quiet,
                                            remap = remap,
                                            no_formula_check = no_formula_check )


def produce_validation_plots( files, verbose_lbls = True, pdf_target = None,
                              **plot_kwargs ):
    """Function which can produce several of the validation plots resulting
    from usage of the produce_validation_plot(..) function, and possibly even
    embed them in a PDF file. Any plot_kwargs will be passed along to the
    produce_validation_plot(..) function.
    """
    from .plot import ( _import_matplotlib_plt,
                        _import_matplotlib_pdfpages )

    if pdf_target:
        pdfpages = _import_matplotlib_pdfpages()
        the_real_pdfpages = ( pdfpages.the_real_inspected_object
                              if hasattr(pdfpages,'the_real_inspected_object')
                              else None )
        already_pdfpages = isinstance(pdf_target,the_real_pdfpages) if the_real_pdfpages else False
        pdf = pdf_target if already_pdfpages else pdfpages(pdf_target)
        assert pdf
    else:
        pdf, already_pdfpages = None, False

    if pdf:
        plt = _import_matplotlib_plt()
    for f in files:
        produce_validation_plot( f, verbose_lbls = verbose_lbls, do_show = not pdf,
                                 line_width_scale = 0.5, **plot_kwargs )
        if pdf:
            pdf.savefig()
            plt.close()

    if pdf and not already_pdfpages:
        pdf.close()
        _nc_common.print("created %s"%pdf_target)

def produce_validation_plot( data_or_file, verbose_lbls = True, line_width_scale = 1,
                             quiet = False, xlabel = None, legend_args = None, do_newfig = True,
                             do_show = True, do_legend = True, do_grid=True, do_tight_layout=True ):

    """Function which can produce validation plots for NCMAT data defining
    crystalline materials like the ones on
    https://github.com/mctools/ncrystal/wiki/Data-library. The plot will compare
    the powder Bragg diffraction cross section curve of the NCMAT file, with
    those obtained by replacing the crystal structure with a crystal structure
    loaded directly from online database entries that might be mentioned in the
    comments of the NCMAT data. This makes it easy to verify that the crystal
    structure is still compatible (or not) with those crystal structures
    mentioned in the comments.

    Currently two online databases are supported, with recognised URLs in a
    format like either of the following:

       https://www.crystallography.net/cod/9008460.html
       https://www.materialsproject.org/materials/mp-87

    """

    def lookupAndProduce( cifsrc, *,dynamics, remap):
        cifloader = CIFLoader( cifsrc, quiet = quiet )
        ncmat = cifloader.create_ncmat_composer( uiso_temperature = None,#dynamics provided separately here
                                                 remap = remap,
                                                 skip_dyninfo = True )
        ncmat.transfer_dyninfo_objects( dynamics )
        return ncmat.create_ncmat()

    from .plot import _import_matplotlib_plt
    plt = _import_matplotlib_plt()
    if do_newfig:
        plt.figure()

    def multcreate( data ):
        return _nc_core.directMultiCreate(data,cfg_params='comp=bragg')
    _file = None
    fn = None
    if hasattr( data_or_file, '__fspath__' ):
        pass
    elif isinstance( data_or_file, _nc_core.TextData ):
        mc = multcreate( data_or_file )
        contentiterable = data_or_file
        fn = data_or_file.dataSourceName
    elif isinstance( data_or_file, bytes ) or isinstance( data_or_file, str ):
        _ = data_or_file.decode() if isinstance( data_or_file, bytes ) else data_or_file
        if '\n' in _ or _.startswith('NCMAT'):
            contentiterable = _.splitlines()
            mc = multcreate(_)
            fn = 'Anonymous NCMAT data'
    if not fn:
        fn = data_or_file
        contentiterable = _nc_core.createTextData(data_or_file)
        mc = multcreate(contentiterable)

    if hasattr(fn,'__fspath__'):
        import pathlib
        fn = pathlib.Path(fn).name
    elif '/' in fn:
        fn = fn.split('/')[-1]
    if not quiet:
        _nc_common.print(f'Attempting to validate {fn}')

    if not mc.info.isCrystalline():
        if not quiet:
            _nc_common.print(f"Ignoring non-crystalline material {fn}")
        return

    def _extractID(s,pattern):
        if pattern not in s:
            return None
        ll=s.split(pattern)[1:]
        while ll:
            e,ll = ll[0],ll[1:]
            d=''
            while e and e[0].isdigit():
                d+=e[0]
                e=e[1:]
            if d and int(d)>0:
                yield int(d)

    import re as _re
    _re_atomdbspecs = _re.compile(r"\[ *with ([a-zA-Z ]+)->([ a-zA-Z0-9+-\.]+) *\]")
    def _extractAtomDBSpec(s):
        #Look for remapping specs like "[with H->D]" and return in @ATOMDB format
        #(i.e. "H is D").
        if '->' not in s:
            return
        m = _re_atomdbspecs.search(s)
        return ' '.join(('%s is %s'%m.groups()).split()) if m else None

    ids = []
    for ll in contentiterable:
        atomdb = _extractAtomDBSpec(ll)
        for e in _extractID(ll,'materialsproject.org/materials/mp-'):
            ids.append( dict(dbtype='mp',entryid=e,atomdb=atomdb))
        for e in _extractID(ll,'crystallography.net/cod/'):
            ids.append( dict(dbtype='cod',entryid=e,atomdb=atomdb))

    #order-preserving remove duplicates:
    _seen = set()
    newids=[]
    for d in ids:
        key = tuple(sorted((k,v) for k,v in d.items()))
        if key in _seen:
            continue
        newids.append(d)
        _seen.add(key)
    ids = newids

    def _atomdb_to_remap( atomdb):
        _atomdb = list( ' '.join(e.strip().split())
                        for e in (atomdb or '').replace(':',' ').split('@') )
        _atomdb = list( e for e in _atomdb if e )
        ll = []
        for c in _atomdb:
            p=c.replace(':',' ').split()
            if not len(p)>=3 or p[1]!='is':
                raise _nc_core.NCBadInput('invalid atomdb remap syntax in "%s"'%c)
            ll.append( (p[0],' '.join(p[2:]) ) )
        return ll

    dynamics = mc.info
    cmps = [ (fn, mc ) ]
    for _ in ids:
        dbname = _['dbtype']
        entryid = _['entryid']
        atomdb = _['atomdb']
        lpargs=dict(remap=_atomdb_to_remap(atomdb),dynamics=dynamics)
        if dbname=='mp':
            lpargs['cifsrc']='mpid::%i'%entryid
        else:
            assert dbname=='cod'
            lpargs['cifsrc']='codid::%i'%entryid
        out = lookupAndProduce( **lpargs )
        lbl=f'{dbname}-{entryid}'
        if atomdb:
            lbl = f'{lbl} with replacement "{atomdb}"'
        cmps.append( ( lbl, multcreate(out) ) )
    from .constants import ekin2wl
    wlmax = max(1.05*ekin2wl(mc.scatter.domain()[0]) for _,mc in cmps)
    from ._numpy import _np_linspace
    wls = _np_linspace(0.001,wlmax,10000)
    if len(cmps)<=1:
        _nc_common.warn(f"No online DB IDs found in {fn}")
        return

    sgno_all = set()
    natoms_all = set()

    col_ordered = [_nc_common._palette_Few.get(k,k) for k in
                   ('black',
                    'blue',
                    'orange',
                    'green',
                    'red',
                    'brown',
                    'purple',
                    'yellow',
                    'pink',
                    'gray')]

    def fix_lbl_for_plt(lbl):
        return lbl.replace('_',"$"+'\\'+"mathrm{\\_}$")

    for i,(lbl,mc) in enumerate(cmps):
        lbl = str(lbl)
        if '/' in lbl:
            lbl=lbl.split('/')[-1]
        elif lbl.startswith('mp-'):
            lbl = 'Materials Project entry '+lbl[3:]
        elif lbl.startswith('cod-'):
            lbl = 'Crystallography Open Database entry '+lbl[4:]
        lw=4 if i>0 else 2
        lw *= line_width_scale

        si = mc.info.structure_info if mc.info.hasStructureInfo() else None
        if not quiet:
            _nc_common.print('Structure[%s]:'%lbl,si)
        if si and verbose_lbls:
            sgno=si.get('spacegroup',None)
            lbl += ' (SG-%s, %i atoms/cell, Vcell=%g)'%(sgno or 'unspecified',si['n_atoms'],si['volume'])
            sgno_all.add(sgno)
            natoms_all.add(si['n_atoms'])

        plt.plot(wls,mc.scatter.xsect(wl=wls),
                 label=fix_lbl_for_plt(lbl),
                 linewidth=lw,
                 alpha=0.5 if i==0 else 0.5,
                 color = col_ordered[i] if i<len(col_ordered) else None,
                 dashes = [] if i==0 else [4 if len(cmps)>2 else 2,2]+[2,2]*i)
    if len(natoms_all)>1:
        _nc_common.warn('WARNING validation detected different number of atoms/cell in tested curves')
    if len(sgno_all)>1:
        _nc_common.warn('WARNING validation detected different space group numbers in tested curves')

    plt.xlim(0.0)
    plt.ylim(0.0)
    plt.title('NB: This compares the crystal structure (space group, lattice, atom positions). Phonons/dynamics always taken from .ncmat file',fontsize=6)
    plt.xlabel(xlabel or 'Neutron wavelength (%s)'%(b'\xc3\x85'.decode()))
    plt.ylabel('Coherent elastic cross section (barn)')
    if do_legend:
        plt.legend(loc='best',handlelength=5,**(legend_args or {}))
    if do_grid:
        plt.grid()
    if do_tight_layout:
        plt.tight_layout()
    if do_show:
        plt.show()

def _extract_descr_from_cif( raw_cif_dict, cifsrc, ciftextdata ):

    def normalise_result( _l ):
        #split any newlines and trim whitespace before returning:
        res = []
        for e in _l:
            for ee in e.splitlines():
                res.append( ee.rstrip() )
        return res

    if not raw_cif_dict:
        _nc_common.warn('Could not properly decode CIF metadata due to inability to access raw CIF data as dictionary.')
        return [], None, None, None

    #cif is cif section with publication info:
    ll=[]
    def extract(key,*altkeys,expectlist = False):
        if not any( key in s for _,s in sorted(raw_cif_dict.items())):
            if altkeys:
                return extract( altkeys[0], *altkeys[1:], expectlist = expectlist )
            return [] if expectlist else ''
        s = list( d[key] for _,d in sorted(raw_cif_dict.items()) if key in d)[0]
        if not expectlist:
            if s is None:
                return ''
            if not hasattr(s,'strip') and hasattr(s,'__len__') and len(s)==1:
                s=s[0]#convert single element list to first element
            s=str(s).strip()
            return '' if s in ('?','.') else s
        else:
            _l=list( ('' if (e is None or (hasattr(e,'strip') and e.strip()=='?')) else str(e).strip()) for e in s )
            return list( e for e in _l if e )

    title = extract('_publ_section_title','_citation_title')
    authors = extract('_publ_author_name','_citation_author_name',expectlist=True)
    journalname = extract('_journal_name_full','_citation_journal_full')
    year = extract('_journal_year','_citation_journal_year','_citation_year')
    doi = extract('_journal_paper_doi')
    chemformsum = extract('_chemical_formula_sum','_cod_original_formula_sum')

    extracted_codid = extract('_cod_database_code')
    if extracted_codid:
        extracted_codid = int(extracted_codid) if extracted_codid.isdigit() else None

    extracted_mpid = None
    _sp = 'Note from NCrystal.cifutils: Data from materialsproject.org / mp-'# NB this exact format is produced elsewhere in this file.
    for e in ciftextdata.splitlines():
        if _sp in e:
            _ = e.split(_sp,1)
            if len(_) >= 2:
                _ = _[1].split()[0].split('#',1)[0]
                if _.isdigit():
                    extracted_mpid = int(_)

    if extracted_codid is not None and cifsrc.codid is not None and cifsrc.codid != extracted_codid:
        _nc_common.warn(f'Embedded CODID {extracted_codid} is different than requested {cifsrc.codid}!')

    if extracted_mpid is not None and cifsrc.mpid is not None and cifsrc.mpid != extracted_mpid:
        _nc_common.warn(f'Embedded MPID {extracted_mpid} is different than requested {cifsrc.mpid}!')

    _thecodid = cifsrc.codid or extracted_codid
    _thempid = cifsrc.mpid or extracted_mpid

    audit_creation_method = extract('_audit_creation_method')
    audit_creation_date = extract('_audit_creation_date')
    audit_author_name = extract('_audit_author_name')

    nauthors = len(authors or [])
    if authors:
        if len(authors)==1:
            authors = authors[0]
        elif len(authors)==2:
            authors = f'{authors[0]} and {authors[1]}'
        else:
            authors = f'{authors[0]}, et al.'
    if title:
        ll.append(f'"{title}"')
    if  (journalname and year) and not authors:
        authors = '<unknown authors>'
    if journalname and year:
        _jy = f'{journalname}, {year}'
    elif journalname or year:
        _jy = journalname or year
    else:
        _jy = ''
    if authors:
        ll.append(f'{authors} [{_jy}]' if _jy else f'Author{"" if nauthors==1 else "s"}: {authors}')

    if doi:
        ll.append(f'DOI: https://dx.doi.org/{doi}')

    if audit_creation_method:
        ll.append(f'CIF creation method: {audit_creation_method}')
    if audit_creation_date:
        ll.append(f'CIF creation date: {audit_creation_date}')
    if audit_author_name:
        ll.append(f'CIF created by: {audit_author_name}')

    if _thecodid:
        ll += [ f'Crystallography Open Database entry {_thecodid}', _codid2url(_thecodid) ]
    if _thempid:
        ll += [ 'The Materials Project', _mpid2url(_thempid) ]

    return normalise_result( ll ), _thecodid, _thempid, chemformsum

def _codid2url( codid ):
    return f'https://www.crystallography.net/cod/{codid}.html'

def _mpid2url( mpid ):
    return f'https://www.materialsproject.org/materials/mp-{mpid}'

def _impl_create_ncmat_composer( cifloader, *,
                                 uiso_temperature,
                                 skip_dyninfo,
                                 top_comments,
                                 fallback_debye_temp,
                                 #merge_equiv,
                                 quiet,
                                 remap,
                                 no_formula_check ):
        with _nc_common.WarningSpy( block = quiet ) as extra_ana_warnings:
            composer = _impl_create_ncmat_composer_internal( cifloader = cifloader,
                                                             #merge_equiv = merge_equiv,
                                                             uiso_temperature = uiso_temperature,
                                                             skip_dyninfo = skip_dyninfo,
                                                             top_comments = top_comments,
                                                             remap = remap,
                                                             no_formula_check = no_formula_check )
            if not skip_dyninfo and fallback_debye_temp and fallback_debye_temp > 0.0:
                composer.allow_fallback_dyninfo( fallback_debye_temp )

        def fmtwarnings( warnings, descrtxt, nwmax ):
            ll=[]
            if not warnings:
                return ll
            ll.append(f'Notice: The following WARNINGS were emitted {descrtxt}:')
            ll.append('')

            wleft=warnings[::-1]
            while wleft:
                w = wleft.pop()
                ncount = 1 + wleft.count(w)
                if ncount > 1:
                    wleft = [e for e in wleft if e!=w]
                wmsg,wcat = w
                wmsg = list(e.strip() for e in wmsg.splitlines() if e.strip())
                s0 = '  %s :'%wcat if ncount==1 else '  (%ix) %s : '%(ncount,wcat)
                for i,m in enumerate(wmsg):
                    if i==nwmax:
                        break
                    pref=s0 if i==0 else ' '*(len(s0))
                    if i+1==nwmax and len(wmsg)>nwmax:
                        m='<%i lines of output hidden>'%(len(wmsg)-nwmax)
                    ll.append('%s %s'%(pref,m))
            return ll

        composer.add_comments( fmtwarnings( list(cifloader.warnings) + extra_ana_warnings,
                                            'when loading the CIF data',nwmax=10),
                               add_empty_line_divider = True )
        return composer

def _aniso_nontrivial( aniso_dict ):
    return aniso_dict and ( any( abs(aniso_dict[e])>1e-13 for e in ('u12','u13','u23') )
                            or ( aniso_dict['u11'] != aniso_dict['u22'] )
                            or ( aniso_dict['u11'] != aniso_dict['u33'] ) )


def _impl_create_ncmat_composer_internal( cifloader, *, uiso_temperature, skip_dyninfo, top_comments, remap, no_formula_check ):

    src = cifloader

    src_atoms = src.atoms

    #Apply remap:
    def _stdatomname( name ):
        return {'D':'H2','T':'H3'}.get(name,name)

    remap_decoded = {}
    for elemiso, compos in (remap or []):
        elemiso, compos = _nc_ncmatimpl._decode_composition(elemiso,compos)
        elemiso = _stdatomname( elemiso )
        remap_decoded[ elemiso ] = compos
    remap = remap_decoded
    remap_str = ''

    if remap:
        new_src_atoms = []
        for atom in src_atoms:
            if not any( e in remap for fr,e in atom['composition'] ):
                new_src_atoms.append( atom )
                continue
            newcompos = []
            for fr,e in atom['composition']:
                remap_compos = remap.get(e,None)
                if remap_compos is None:
                    newcompos.append( (fr,e) )
                else:
                    for remap_fr, remap_elem in remap_compos:
                        newcompos.append( ( fr*remap_fr, remap_elem ) )
                d = dict( (k,v) for k,v in atom.items() if k!='composition' )
                d['composition'] = newcompos
                new_src_atoms.append( d )
        src_atoms = tuple( new_src_atoms )
        #For the description:
        remap_strs = []
        for elemiso, compos in sorted(remap.items()):
            compos_str = compos[0][1] if len(compos)==1 else ' '.join( f'{_f:g} {_n}' for _f,_n in compos)
            remap_strs.append(f'{elemiso} -> {compos_str}')
        remap_str = ', '.join(remap_strs)

    #DYNINFO analysis can only be done now, where we know the skip_dyninfo and uiso_temperature values:

    ncmat = _nc_ncmat.NCMATComposer()
    all_atompos = []
    for idx, a in enumerate( src_atoms ):
        lbl = f'cif_species_{idx}'
        compos = [ (fr,e) for fr,e in a['composition'] ]
        if len(compos)==1 and compos[0][0]==1.0:
            composstr = compos[0][1]
        else:
            composstr = ' '.join(f'{f:g} {n}' for f,n in compos)
        ncmat.set_composition( lbl, compos )
        if not skip_dyninfo:
            uiso = a.get('uiso',None)
            if uiso is not None:
                assert uiso>0.0
                if uiso_temperature is None:
                    _nc_common.warn(f'ignoring uiso info present in CIF input for "{composstr}" since uiso_temperature parameter value is not provided\n')
                else:
                    ncmat.set_dyninfo_msd( lbl, msd=uiso, temperature=uiso_temperature )
            _aniso = a.get('aniso',None)

            if _aniso_nontrivial(_aniso):
                _nc_common.warn('Anisotropic displacement for %s ignored (%s)'%(composstr,', '.join (f'{k}={v:g}' for k,v in sorted(_aniso.items()))))

        occu = a.get('occupancy',1.0)
        for c in a['equivalent_positions']:
            all_atompos.append( ( lbl, c[0], c[1], c[2], occu ) )

    ncmat.set_atompos( all_atompos )

    st = src.cellsg
    ncmat.set_cellsg( a=st['a'], b=st['b'], c=st['c'],
                      alpha=st['alpha'], beta=st['beta'], gamma=st['gamma'],
                      spacegroup=st['spacegroup']['number'] )

    if top_comments is not None:
        ncmat.add_comments(top_comments)

    _extracted_description = src.extracted_description
    if not _extracted_description:
        if src.cifsrc.name:
            _ds = 'CIF data' if src.cifsrc.is_remote else 'CIF file'
            _extracted_description = [f'{_ds}: {src.cifsrc.name}']
        else:
            _extracted_description = ['Anonymous CIF data']
    if _extracted_description:
        ncmat.add_comments(['Structure converted (with NCrystal'
                            '.cifutils module) from:',''])

        _thedescr = _extracted_description
        if remap_str:
            src.actual_mpid
            _codidurl = _codid2url(src.actual_codid) if src.actual_codid else None
            _mpidurl = _mpid2url(src.actual_mpid) if src.actual_mpid else None
            if _codidurl and _codidurl in _thedescr:
                _=f'{_codidurl} [with {remap_str}]'
                _thedescr = [e.replace(_codidurl,_) for e in _thedescr]
            elif _mpidurl and _mpidurl in _thedescr:
                _=f'{_mpidurl} [with {remap_str}]'
                _thedescr = [e.replace(_mpidurl,_) for e in _thedescr]
            else:
                _thedescr+=['',f'Note: With custom remappings: {remap_str}']

        ncmat.add_comments(['  '+e for e in _thedescr])
        ncmat.add_comments([''])

    ncmat.add_comments('IMPORTANT NOTICE: This is a mostly automatic conversion which has not been',add_empty_line_divider=True)
    ncmat.add_comments('                  verified!  In particular the @DYNINFO sections might need')
    ncmat.add_comments('                  post-editing. Before distributing this file to other people,')
    ncmat.add_comments('                  please review this, amend the comments here to document,')
    ncmat.add_comments('                  anything done, and remove this notice.')

   #formula
    def _extract_atomcount( atom ):
        occu, npos = atom['occupancy'], len(atom['equivalent_positions'])
        return npos * ( int(occu) if occu == int(occu) else float(occu) )
    total_composition = []
    for a in src_atoms:
        atomcount = _extract_atomcount( a )
        for elemfrac, eleminfo in a['composition']:
            total_composition.append( ( eleminfo, elemfrac * atomcount ) )
    formula = _nc_common.format_chemform( total_composition )


    expected_formula = cifloader.raw_cif_chemformula

    if expected_formula:
        def formula_to_dict( f ):
            if not f:
                return None
            if isinstance(f,str):
                return formula_to_dict( _decode_formula(f) )
            if isinstance(f,dict):
                return f
            res = {}
            for lbl,count in f:
                assert count>=0.0
                if count==0:
                    continue
                if lbl in res:
                    res[lbl] += count
                else:
                    res[lbl] = count
            return res or None
        from .atomdata import allElementNames
        _all_elements = allElementNames()
        def _decode_formula( s ):
            l1 = [e for e in _all_elements if len(e)==1]+['D','T']
            l2 = [e for e in _all_elements if len(e)==2]
            res = []
            while s.strip():
                s=s.strip()
                c = [ e for e in l2 if s.startswith(e) ]
                if not c:
                    c = [ e for e in l1 if s.startswith(e) ]
                if len(c)!=1:
                    return None
                c = c[0]
                s=s[len(c):].strip()
                leaddigits = ''
                while s and s[0].isdigit():
                    leaddigits += s[0]
                    s=s[1:]#nb: no strip here, we want to stop on a space!
                n = int(leaddigits) if leaddigits else 1
                res.append( (c,n) )
            return formula_to_dict(res)
        def formulas_incompatible( formula1, formula2 ):
            f1 = formula_to_dict(formula1)
            f2 = formula_to_dict(formula2)
            if not f1 or not f2:
                return False
            #map D,T,H2,H3 -> H to ensure fewer false positives:
            def _tmp(d):
                return dict( (('H' if k in ('D','T','H2','H3') else k),v)
                             for k,v in d.items() )
            f1,f2 = _tmp(f1),_tmp(f2)
            if f1 == f2:
                return False
            if not set(f1.keys())==set(f2.keys()):
                return True
            kref = list(f1.keys())[0]
            assert f1[kref] > 0.0 and f2[kref]>0.0
            f2scale = f1[kref]/f2[kref]
            tol = 1e-2#not tighter to avoid false positives!
            for k,v1 in f1.items():
                v2 = f2[k]*f2scale
                if abs(v1-v2)>tol and abs(v1-v2)/(1e-300+abs(v1)+abs(v2)) > tol:
                    return True
            return False

        if not no_formula_check and expected_formula:
            actual_expected_formula_dict = formula_to_dict( expected_formula )
            if remap:
                newf = {}
                def _collect( _name, _frac ):
                    _name = _stdatomname(_name)
                    if _name not in newf:
                        newf[_name] = _frac
                    else:
                        newf[_name] += _frac
                for name, frac in actual_expected_formula_dict.items():
                    _remapped = remap.get(_stdatomname(name),None)
                    if _remapped is not None:
                        for remap_frac, remap_name in _remapped:
                            _collect( remap_name, remap_frac * frac )
                    else:
                        _collect( name, frac )
                actual_expected_formula_dict = newf
            if formulas_incompatible(total_composition,actual_expected_formula_dict):
                s = f'"{expected_formula}"'
                if remap:
                    _ = _nc_common.format_chemform( list( sorted(actual_expected_formula_dict.items() ) ) )
                    s += f' remapped to "{_}"'
                raise _nc_core.NCBadInput(f'Formula encoded in CIF data ({s}) is not compatible with formula of loaded structure ("{formula}")')

    return ncmat


def _impl_merge_atoms( atoms ):
    ll = []
    for a in atoms:
        pos = list(a['equivalent_positions'])
        cif_labels = list( a['cif_labels'] )
        other_metadata = list( sorted( (k,v) for k,v in a.items()
                                       if k not in ('equivalent_positions',
                                                    'cif_labels') ) )
        found = False
        for k,v in ll:
            if k == other_metadata:
                v[0] += list( pos )
                v[1] += list( cif_labels )
                found = True
                break
        if not found:
            ll.append( (other_metadata,[list(pos),list(cif_labels)]) )
    res = []
    for other_metadata, ( pos, cif_labels ) in ll:
        d = dict( (k,v) for k,v in sorted(other_metadata) )
        d['equivalent_positions'] = list( sorted( pos ) )
        d['cif_labels'] = list( sorted( cif_labels ) )
        res.append( d )
    return res

def _suggest_filename( ncmat_metadata, cifloader ):
    ll = ['autogen']
    ll.append( ncmat_metadata['chemform'] )
    sgnum = ncmat_metadata.get('cellsg',{}).get('spacegroup',None)
    if sgnum:
        ll.append( 'sg%i'%sgnum )
    if cifloader.actual_codid:
        ll.append( 'cod%i'%cifloader.actual_codid )
    if cifloader.actual_mpid:
        ll.append( 'mp%i'%cifloader.actual_mpid )
    return '_'.join(ll) + '.ncmat'

def _format_spglib_cell( cellsg, atoms ):
    lattice = _nc_ncmatimpl._cellparams_to_spglib_lattice(cellsg)
    atomic_points = []
    atomic_types = []
    for i,a in enumerate(atoms):
        for p in a['equivalent_positions']:
            atomic_points.append( p )
            atomic_types.append( i )
    return ( lattice, atomic_points, atomic_types )

def _impl_refine_cell( cellsg, atoms ):
    orig_cell = _format_spglib_cell( cellsg, atoms )

    d = _nc_ncmatimpl._spglib_refine_cell( orig_cell )#, symprec = 0.01, allow_axis_swap = True )
    assert len(d)==7

    refined_cell = d['refined_cell']
    new_atom_pos = dict( (i,[]) for i in range(max(refined_cell[2])+1) )
    for pos, atomidx in zip(refined_cell[1],refined_cell[2]):
        new_atom_pos[ atomidx ] += [ (pos[0],pos[1],pos[2]) ]

    new_atoms = [ dict( (k,v if k!='equivalent_positions' else new_atom_pos[idx])
                        for k,v in e.items()) for idx,e in enumerate(atoms) ]
    if not d['can_keep_anisotropic_properties']:
        new_atoms = [ dict( (k,v if k!='aniso' else None)
                            for k,v in e.items()) for idx,e in enumerate(new_atoms) ]




    new_cellsg = dict( d['cellparams_snapped'].items() )
    new_cellsg['spacegroup'] = dict( number = d['sgno'],
                                     hm = d['sgsymb_hm'] )

    return new_cellsg, new_atoms, d['warnings'], d['msgs']

def _gemmi_wrap_to_unit( gemmi_fractional ):
    #gemmi's wrap_to_unit() can return 1.0, but we want those to be 0.0 (so
    #(a,b,1.0) and (a,b,0.0) does not appear to be two different
    #points). Update: I am not 100% sure about this, but keeping this function
    #for added robustness.
    c = gemmi_fractional.wrap_to_unit()
    for i in range(3):
        if c[i]==1.0:
            c[i]=0.0
    return c

def _getMaterialsProjectAPIKEY():
    import os
    apikey = os.environ.get('MATERIALSPROJECT_USER_API_KEY',None)
    if len(apikey.strip()) < 25:
        raise _nc_core.NCException('ERROR: Legacy API key found in MATERIALSPROJECT_USER_API_KEY. Please instead provide a new API key from https://materialsproject.org/api (Keys for the new API are ~32 characters, whereas keys for the legacy API are ~16 characters).')
    if not apikey:
        raise _nc_core.NCException('ERROR: Missing API key for materialsproject.org access. To fix, '
                           +'please make sure the environment variable MATERIALSPROJECT_USER_API_KEY'
                           +' contains your personal access key that you see on'
                           +' https://www.materialsproject.org/dashboard (after logging in).')
    return apikey

def _use_local_cif_cache( fn, text_data = None, quiet = False ):
    #NB: This simple implementation use no locking to guard against race
    #conditions!!! But we do perform write+move instead of simply write, which
    #is a bit more "atomic".
    if _nc_common.ncgetenv_bool('ONLINEDB_FORBID_NETWORK'):
        def notfound():
            n = _nc_common.expand_envname('ONLINEDB_FORBID_NETWORK')
            raise RuntimeError('Error: Trying to access remote DB but'
                               f' {n} is set')
        time_limit_hours = 24*7*365*1000#revisit this in 3023
    else:
        def notfound():
            pass
        time_limit_hours = 24*7
    d = _nc_common.ncgetenv('ONLINEDB_CACHEDIR')
    if not d:
        return notfound()
    import os.path
    import pathlib
    p = pathlib.Path(os.path.expanduser(d)).resolve().absolute() if d else None
    if not p:
        return notfound()
    pfn = ( p / fn )
    if text_data:
        #store data into local cache:
        if not quiet:
            _nc_common.print(f"Adding {fn} to local file cache in $NCRYSTAL_ONLINEDB_CACHEDIR")
        p.mkdir(parents=True, exist_ok=True)
        pfn_tmp = p / f'{fn}_tmp_{os.getpid()}'
        _nc_common.write_text(pfn_tmp,text_data)
        pfn_tmp.replace(pfn)
        return notfound()
    #retrieve data from local cache:
    if not pfn.exists():
        return notfound()
    #Since online data can change, entries expire eventually:
    import datetime
    timelim = datetime.timedelta( hours = time_limit_hours )
    now = datetime.datetime.now(tz=datetime.timezone.utc)
    mtime = datetime.datetime.fromtimestamp(pfn.stat().st_mtime, tz=datetime.timezone.utc)
    if ( now-mtime ) > timelim:
        #entry expired:
        if not quiet:
            _nc_common.print(f"Ignoring (and removing) expired {fn} from local file cache in $NCRYSTAL_ONLINEDB_CACHEDIR")
        pfn.unlink()
        return notfound()
    if not quiet:
        _nc_common.print(f"Getting {fn} from local file cache in $NCRYSTAL_ONLINEDB_CACHEDIR")
    return pfn.read_text()

_cod_cache = []
def _cod_get_cifdata( codid, quiet = False ):
    for _codid, _result in _cod_cache:
        if codid == _codid:
            if not quiet:
                _nc_common.print(f"Using cached Crystallography Open Database result for entry {codid}")
            return _result
    cache_fn = 'cod_%i.cif'%codid
    #check file cache:
    c = _use_local_cif_cache( cache_fn, quiet = quiet )
    if c:
        return c
    if not quiet:
        _nc_common.print(f"Querying the Crystallography Open Database for entry {codid}")

    mirror_urls = [
        "https://www.crystallography.net/cod/%i.cif",#canonical first
        'https://qiserver.ugr.es/cod/%i.cif',
        'http://cod.ibt.lt/cod/%i.cif',
    ]

    result = None
    while result is None:
        #Fail gently and try next mirror - except for the last attempt.
        url = mirror_urls.pop(0)%codid
        result = _nc_common.download_url( url,
                                          timeout = 10.0,
                                          quiet_network_fail = bool(mirror_urls) )
        if result or not mirror_urls:
            break
        nextdescr = '/'.join(mirror_urls[0].split('/')[0:-1])+'/'
        _nc_common.print(f'Retrival failed. Trying mirror at: {nextdescr}')

    assert result is not None

    if len(_cod_cache)==10:
        _cod_cache.pop(0)
    _cod_cache.append( (codid, result ) )
    #update file cache:
    _use_local_cif_cache( cache_fn, quiet = quiet, text_data = result )
    return result

_mp_cache = []
def _mp_get_cifdata( mpid, quiet = False, apikey = None ):
    if hasattr(mpid,'startswith') and mpid.startswith('mp-'):
        mpid = mpid[3:]
    if hasattr(mpid,'startswith') and mpid.startswith('mpid::'):
        mpid = mpid[6:]
    mpid = int(mpid)
    for _mpid, _result in _mp_cache:
        if mpid == _mpid:
            if not quiet:
                _nc_common.print(f"Using cached materialsproject.org result for entry mp-{mpid}")
            return _result
    #check file cache:
    cache_fn = 'mp_%i.cif'%mpid
    c = _use_local_cif_cache( cache_fn, quiet = quiet )
    if c:
        return c

    if not quiet:
        _nc_common.print(f"Querying materialsproject.org for entry mp-{mpid}")

    if apikey is None:
        apikey = _getMaterialsProjectAPIKEY()

    with _nc_common.WarningSpy(blockfct = lambda msg, cat : cat in ('PendingDeprecationWarning','DeprecationWarning') ):
        try:
            import mp_api.client#NB: This might trigger a spurious FPE
        except ImportError:
            raise ImportError('Could not import mp_api.client. Installing the mp-api package will most likely'
                              ' fix this (perhaps with a command like "conda install -c conda-forge mp-api"'
                              ' or "python3 -mpip install mp-api").')

    with _nc_common.WarningSpy(blockfct = lambda msg,cat : msg.lower().startswith('mpcontribs-client not installed') ):
        with mp_api.client.MPRester(apikey) as mpr:
            s = mpr.get_structure_by_material_id( f'mp-{mpid}', conventional_unit_cell=True )
            #If we do no use the symprec argument in the next line, we will get a an unrefined P1 structure:
            result = s.to(fmt='cif',symprec=1e-4, significant_figures=15, angle_tolerance=5.0, refine_struct=True)
            #But we want to sanity check that this refinement gives the same spacegroup result as listed in the MP database:
            mp_expected_sg_number = s.get_space_group_info()[1]
    sg_checked = False
    errmsg = f'Unable to reliably determine spacegroup when trying to retrieve structure for mp-{mpid} from materialsproject.org'
    for ll in result.splitlines():
        p = ll.split('#',1)[0].split()
        if p and p[0]=='_symmetry_Int_Tables_number':
            if not len(p)>=2 or not p[1].isdigit():
                raise _nc_core.NCBadInput(errmsg+f' (unexpected format of line: "{ll}")')
            _sgnum = int(p[1])
            if _sgnum == mp_expected_sg_number:
                sg_checked = True
            else:
                raise _nc_core.NCBadInput(errmsg+f' (expected SG-{mp_expected_sg_number} but got SG-{_sgnum})')

    if not sg_checked:
        raise _nc_core.NCBadInput(errmsg+' (no line with "_symmetry_Int_Tables_number" produced)')

    #Finally embed origin as a note (so we can extract it later for the ncmat
    #header comments). Keep the format below synchronised with the reader
    #elsewhere in this file!!
    result += f'\n# Note from NCrystal.cifutils: Data from materialsproject.org / mp-{mpid}\n'

    #Update cache and return:
    if len( _mp_cache ) == 10:
        _mp_cache.pop(0)
    _mp_cache.append( (mpid, result ) )

    #update file cache:
    _use_local_cif_cache( cache_fn, quiet = quiet, text_data = result )

    return result

def _codid2url( codid ):
    return f'https://www.crystallography.net/cod/{codid}.html'
def _mpid2url( mpid ):
    return f'https://www.materialsproject.org/materials/mp-{mpid}'

####################### NEW GEMMI STUFF ############################################
_import_gemmi_cache = [None,None]
def _import_gemmi( *, sysexit = False ):
    if _import_gemmi_cache[0] is not None:
        return _import_gemmi_cache[0], _import_gemmi_cache[1]
    try:
        import gemmi#both available on pypi and conda-forge
        import gemmi.cif
    except ImportError:
        m = ( 'Could not import gemmi modules needed to process CIF files.'
              +' The gemmi package is available on both PyPI ("python3 -mpip install'
              +' gemmi") and conda ("conda install -c conda-forge gemmi")' )
        if sysexit:
            raise SystemExit(m)
        else:
            raise ImportError(m)
    _import_gemmi_cache[0] = gemmi
    _import_gemmi_cache[1] = gemmi.cif
    return gemmi, gemmi.cif

_guessmap = [None]
def _guess_spacegroup_name( gemmi, s ):
    def _guess_keys( s ):
        s=''.join(s.split())
        return [ s, s.replace('/',''),s.replace('-',''),s.replace('-','').replace('/','' ) ]
    def _init_guess_map():
        guess = {}
        def _ag( i, s ):
            if s not in guess:
                guess[s] = set([i])
            else:
                guess[s].add(i)
        def add_guess( i, s ):
            for k in _guess_keys(s):
                _ag( i, k )
        for i in range(1,230+1):
            sg = gemmi.find_spacegroup_by_number(i)
            add_guess( i, sg.xhm() )
            add_guess( i, sg.short_name() )
            add_guess( i, sg.hall )
            add_guess( i, sg.hm )
        return guess
    if not _guessmap[0]:
        _guessmap[0] = _init_guess_map()
    gm = _guessmap[0]
    possible_sgnos = set()
    for k in _guess_keys(s):
        for i in gm.get(k,[]):
            possible_sgnos.add( i )
    return list(sorted(possible_sgnos))

def _load_with_gemmi( cifblock, allow_fixup = True ):
    #Load cifblock into gemmi struct. We might perform in-place editing of the
    #cifblock, if gemmi does not immediately recognise the space group.
    gemmi, gemmi_cif = _import_gemmi()

    struct = gemmi.make_small_structure_from_block( cifblock )

    assert struct
    if hasattr(struct,'find_spacegroup'):
        sg = struct.find_spacegroup()
    else:
        sg = struct.spacegroup

    if not sg:
        #Doing what the old struct.find_spacegroup() was doing:
        sg = gemmi.find_spacegroup_by_name( hm=struct.spacegroup_hm,
                                            alpha=struct.cell.alpha,
                                            gamma=struct.cell.gamma)

    if sg or not allow_fixup:
        return struct, sg

    #Let us see if we can find the spacegroup with a bit of manual
    #intervention. If we can, we update the cifblock and do a full reload
    #(otherwise the gemmi struct object won't have the spacegroup and images
    #properly initialised):

    _ = cifblock.find(['_space_group_IT_number'])
    if not _:
        _ = cifblock.find(['_symmetry_Int_Tables_number'])

    if _ and len(_)==1 and len(_[0])==1 and str(_[0][0]).isdigit() and (1<=int(str(_[0][0]))<=230):
        _sgnum = int(_[0][0])
        #NB: We can do this for some groups only (cubic?) Or perhaps we could use it as a double-check only?
        sg = gemmi.find_spacegroup_by_number(_sgnum)
        assert sg and sg.number == _sgnum

    sg_hm = str(struct.spacegroup_hm).strip()
    if sg_hm and not sg:
        def sg_searchname(x):
            return gemmi.find_spacegroup_by_name(x,
                                                 struct.cell.alpha,
                                                 struct.cell.gamma )
        attempts = [ sg_hm ]
        #TODO: We used to do this, but currently it does not make a difference:
        #for e in ('H','R'):#NB H+R is not enough, see table 6 at http://cci.lbl.gov/sginfo/hall_symbols.html
        #    if sg_hm.endswith(e):
        #        attempts.append( sg_hm[:-1] + ' : ' + sg_hm[-1] )
        for a in attempts:
            a = ' '.join(a.split())
            sg = sg_searchname( a )
            if sg:
                _nc_common.warn(f'Had to interpret spacegroup "{sg_hm}" as "{a}" before Gemmi could recognise it.')
                break

    if sg_hm and not sg:
        possible = _guess_spacegroup_name( gemmi, sg_hm )
        if len(possible) == 1:
            sgnum = possible[0]
            sg = gemmi.find_spacegroup_by_number(sgnum)
            _nc_common.warn(f'Had to interpret spacegroup "{sg_hm}" as "{sg.hm}" before Gemmi could recognise it.')
        elif len(possible) > 1:
            poshm = list( (sgnum,gemmi.find_spacegroup_by_number(sgnum).hm) for sgnum in possible )
            poshm = ', '.join( '"%s"(number %i)'  for sgnum,sgstr in poshm )
            _nc_common.warn(f'Failed to interpret spacegroup interpret spacegroup "{sg_hm}". Could be any of: {poshm}.')
    if not sg:
        #manual intervention did not help:
        return struct, None

    #patch up the cifblock with info about manually found spacegroup, and try again:
    _update_spacegroup_in_cifblock( cifblock, sg )

    return _load_with_gemmi( cifblock, allow_fixup = False )

def _update_spacegroup_in_cifblock_from_name( gemmi, cifblock, name ):
    sg = gemmi.find_spacegroup_by_name( name )
    if not sg:
        raise _nc_core.NCBadInput(f'Unknown spacegroup: "{name}"')
    _update_spacegroup_in_cifblock( cifblock, sg )

def _update_spacegroup_in_cifblock( cifblock, sg, use_xhm = True ):
    cifblock.set_pair('_space_group_name_H-M_alt', sg.xhm() if use_xhm else sg.hm )
    cifblock.set_pair('_space_group_IT_number', str(sg.number) )
    cifblock.set_pair('_space_group_name_Hall', sg.hall )
    cifblock.set_pair('_space_group_crystal_system', sg.crystal_system_str() )


def _actual_init_gemmicif( cifsrc, *, quiet, mp_apikey, refine_with_spglib, merge_equiv, override_spacegroup = None ):

    from math import fsum as _math_fsum

    result = {}
    pos_tolerance = 0.0001#NB: Best if value matches the one in NCInfoBuilder.cc

    from types import MappingProxyType
    dict_ro, list_ro = MappingProxyType, tuple

    gemmi, gemmi_cif = _import_gemmi()

    if not isinstance( cifsrc, CIFSource ):
        cifsrc = CIFSource( cifsrc )

    result['cifsrc'] = cifsrc

    _cifdata = cifsrc.load_data( quiet = quiet, mp_apikey = mp_apikey )
    if not quiet:
        _nc_common.print("Attempting to load CIF data with gemmi")

    result['cifdata'] = _cifdata
    try:
        cif_doc = gemmi_cif.read_string( _cifdata )
    except ValueError as e:
        errmsg = 'CIF parsing error (from Gemmi): "%s"'%e
        cif_doc = None
    if cif_doc is None:
        raise _nc_core.NCBadInput(errmsg or 'Unknown CIF parsing error from Gemmi')

    cif_doc_as_json = cif_doc.as_json()
    import json
    try:
        cif_doc_as_dict = json.loads( cif_doc_as_json )
    except json.decoder.JSONDecodeError as e:
        cif_doc_as_dict = None
        _nc_common.warn('Could not decode raw CIF data to dictionary (Gemmi bug?): JSONDecodeError("%s")'%e)

    result['cif_raw'] = cif_doc_as_dict
    if len(cif_doc) == 0:
        raise _nc_core.NCBadInput('CIF data had no blocks!')
    if len(cif_doc) == 1:
        cif_block_with_structure = cif_doc.sole_block()
    else:
        cif_block_with_structure = None
        for e in cif_doc:
            if not len(e.find(['_atom_site_fract_x'])):
                continue
            if cif_block_with_structure:
                cif_block_with_structure = None
                break
            cif_block_with_structure = e
        if not cif_block_with_structure:
            raise _nc_core.NCBadInput('Could not automatically determine block in CIF data with structure info')

    assert cif_block_with_structure

    if override_spacegroup:
        assert isinstance(override_spacegroup,str)
        _nc_common.warn(f'Overriding spacegroup to "{override_spacegroup}" due to explicit request')
        _update_spacegroup_in_cifblock_from_name( gemmi, cif_block_with_structure, override_spacegroup )

    #record original h-m-alt entry:
    orig_hm_alt = None
    _ = cif_block_with_structure.find(['_space_group_name_H-M_alt'])
    if _ and len(_)==1 and len(_[0])==1:
        orig_hm_alt = str(_[0][0]).strip()

    struct, sg = _load_with_gemmi( cif_block_with_structure )

    if not struct:
        raise _nc_core.NCBadInput('Could not load structure with Gemmi')

    if not sg:
        raise _nc_core.NCBadInput('Could not determine space group from CIF data')

    #warn if spacegroup setting might be ambiguous:
    if sg.number != 1:
        _ = list( sorted( (e.number,e.xhm(),e.is_reference_setting()) for e in gemmi.spacegroup_table() if ( e.number==sg.number and ':' in e.xhm())))
        #_ = [ (no,xhm,isref) for no,xhm,isref in _ if isref ]
        if len(_)>1:
            _str = '", "'.join( xhm for no,xhm,isref in _)
            if not orig_hm_alt or not any( (orig_hm_alt==e[1] or orig_hm_alt.replace(' ','')==e[1].replace(' ','')) for e in _ ):
                _nc_common.warn(f'SG-{sg.number} available in multiple'
                                f' choices ("{_str}") and which one was not'
                                ' encoded explicitly in the '
                                '_space_group_name_H-M_alt CIF field. Consider'
                                ' overriding the space group explicitly when'
                                ' loading this file.')

    _ = cif_block_with_structure.find(['_atom_type_number_in_cell'])
    _ = sum((sum(([e] for e in ll),[]) for ll in _),[]) if _ else []
    if _ and not any( e is None for e in _ ):
        expected_tot_atom_in_orig_cell = sum(float(e) for e in _)
    else:
        expected_tot_atom_in_orig_cell = None
    if not ( isinstance(sg.number,int) or sg.number.isdigit() ) or not ( 1<=int(sg.number)<=230 ):
        raise _nc_core.NCBadInput(f'Could not determine space group from CIF data (it loaded with invalid SG number: {sg.number}')
    if not struct.cell.is_crystal():
        raise _nc_core.NCBadInput('Loaded structure is non-crystalline (according to Gemmi).')

    sgnumber = int(sg.number)

    if not struct.cell.is_compatible_with_spacegroup( sg ):
        raise _nc_core.NCBadInput('Loaded unit cell is not compatible with deduced space group (according to Gemmi).')

    cellsg = dict ( a = struct.cell.a, b = struct.cell.b, c = struct.cell.c,
                    alpha = struct.cell.alpha,
                    beta = struct.cell.beta,
                    gamma = struct.cell.gamma )

    fractcoord_approx_1angstrom = 1.0 / ( (struct.cell.a+struct.cell.b+struct.cell.c)/3.0 )


    #Although the sg object also provides sg.ccp4 and sg.hall, we record here
    #just sgnumber and xhm, since that is what spglib refinement also provides:
    cellsg['spacegroup'] = dict_ro( dict( number = sgnumber,
                                          hm = sg.xhm(),
                                         ) )

    collected_atoms = []

    def _expand_coord_to_all_other_images( coord, allow_finetune = 3 ):
        pos0 = _gemmi_wrap_to_unit(coord)
        #First find all brute-force expanded coords (ignore those within machine
        #precision of each other!):
        ll = [ pos0 ]
        for candidate in ( _gemmi_wrap_to_unit(img.apply( pos0 )) for img in struct.cell.images):
            use = True
            for c in ll:
                if _nc_ncmatimpl._unit_cell_point_dist(candidate,c) < 1e-10:
                    use = False
                    break
            if use:
                ll.append( candidate )
        if not allow_finetune:
            return ll
        #Now, check how many of these are very close to the initial point:
        _ucpdist = _nc_ncmatimpl._unit_cell_point_dist
        lclose = [ _nc_ncmatimpl._remap_fract_pos_pt(e) for e in ll
                   if _ucpdist(pos0,e) < 0.01*fractcoord_approx_1angstrom ]
        assert len(lclose) > 0
        if len(lclose) == 1:
            #no issues, just return:
            return ll
        #Input might have had inexact coordinates for atoms at special
        #positions. Fine-tune pos0 as average over the close points, and rerun:
        def _fract_delta( x1, x0 ):
            #remember that distance of (0,0,eps) and (0,0,1-eps) is 2eps, not 1-2eps.
            dx = x1 - x0
            if dx > 0.5:
                dx -= 1.0
            if dx < -0.5:
                dx += 1.0
            return dx
        def _fract_delta_pt( xyz1, xyz0 ):
            return tuple( _fract_delta( xyz1[i], xyz0[i] ) for i in range(3) )

        pos0_new = [ _math_fsum( [lclose[0][i]] + [ _fract_delta(c[i],lclose[0][i])/len(lclose) for c in lclose ] ) for i in range(3) ]
        #snap to (0,0,0):
        for i in range(3):
            if abs(pos0_new[i])< 1e-16:
                pos0_new[i] = 0.0
        def fmt( c ):
            return f'({c[0]:.15g},{c[1]:.15g},{c[2]:.15g})'
        _nc_common.warn('Fractional coordinate %s interpreted as special position %s to avoid numerical precision issues'%(fmt(coord),fmt(pos0_new)))
        return _expand_coord_to_all_other_images( gemmi.Fractional(pos0_new[0],pos0_new[1],pos0_new[2]), allow_finetune = (allow_finetune-1) )


    for site in struct.sites:
        if not float( site.occ ) > 0.0:
            continue
        pos0 = _gemmi_wrap_to_unit( site.fract )
        expanded_coords = _expand_coord_to_all_other_images( pos0 )

        aniso = None
        if site.aniso.nonzero:
            aniso = dict( (k, float(getattr(site.aniso,k))) for k in ('u11','u22','u33','u12','u13','u23') )
            if all( v==0.0 for v in aniso.values() ):
                aniso = None
        u_iso = float( site.u_iso )
        if not u_iso > 0.0:
            u_iso = None

        if aniso and not u_iso:
            #this logic might already exist on the gemmi-side, but to be safe we add it here as well:
            u_iso = ( aniso['u11'] + aniso['u22'] + aniso['u33'] ) / 3

        site_en = str(site.element.name).strip()
        if site_en == 'X' or not site_en:
            #Issue seen in codid::9005777, which had no atomic symbols but
            #labels like "CaM1", "CaM2", "OA1".."OA3", "OB1", "OB2", "OC1",
            #"OC3",... Correct if beginning of string starts with exactly one
            #element name, followed by either a digit or an upper case string
            #followed by digits.
            from .atomdata import allElementNames
            def guess_elem_name( site_label ):
                sl, traildigit = _nc_common._split_trailing_digit( site_label )
                if traildigit is None:
                    return
                #first check element names with 2 chars:
                candidates = [ e for e in allElementNames() if ( sl.startswith(e) and len(e)==2 ) ]
                if len(candidates) > 1:
                    return
                #Then those with 1 char if no hit already (also check for 'D' and 'T'):
                if not candidates:
                    candidates = [ e for e in allElementNames() if ( sl.startswith(e) and len(e)==1 ) ]
                    candidates += [ e for e in ('D','T') if ( sl.startswith(e) and len(e)==1 ) ]
                if len(candidates)!=1:
                    return
                c = candidates[0]
                assert site_label.startswith(c)
                leftover = site_label[len(c):].strip()
                while leftover and leftover[0].isupper()  and leftover[0].isalpha():
                    leftover = leftover[1:].strip()
                if not leftover or leftover.isdigit():
                    return c
            elem_name = guess_elem_name( site.label )
            if elem_name:
                _nc_common.warn(f'Assuming atomic symbol "{elem_name}" based on CIF label "{site.label}"')
            else:
                raise _nc_core.NCBadInput(f'Neither Gemmi nor NCrystal\'s custom'
                                          ' code could deduce an atomic symbol'
                                          f' for the CIF label "{site.label}"')
            elem_Z = None
        else:
            elem_name, elem_Z = site_en, int(site.element.atomic_number)

        elem_A = None
        if elem_name == 'D':
            elem_name, elem_Z, elem_A = 'H',1,2
        if elem_name == 'T':
            elem_name, elem_Z, elem_A = 'T',1,3
        elem_marker = elem_name + ( str(elem_A) if elem_A else '' )

        from .atomdata import elementNameToZValue
        z2name = elementNameToZValue( elem_name, allow_isotopes = False )
        if z2name is None:
            raise _nc_core.NCBadInput(f'Not an element name: "{elem_name}"')

        if elem_Z is None:
            elem_Z = z2name
        if elem_Z != z2name:
            raise _nc_core.NCBadInput('Wrong atomic number (%i) provided for element "%s"'%(elem_Z,elem_name))

        collected_atoms.append( dict( expanded_coords = expanded_coords,
                                      element_marker = elem_marker,
                                      occupancy = float(site.occ),
                                      cif_label = str(site.label),
                                      uiso = u_iso,
                                      aniso = aniso ) )

    #Now merge collected atoms which occupy the same sites:
    def has_same_sites( atom1, atom2 ):
        l1,l2 = atom1['expanded_coords'],atom2['expanded_coords']
        assert len(l1)>0 and len(l2)>0
        return len(l1) == len(l2) and any( ( _nc_ncmatimpl._unit_cell_point_dist(c2,l1[0])<pos_tolerance) for c2 in l2 )

    atoms_to_process, final_atoms = collected_atoms, []
    while atoms_to_process:
        atom1 = atoms_to_process[0]
        atoms_with_same_pos = [ atom1 ]
        atoms_with_different_pos = []
        for atom2 in atoms_to_process[1:]:
            if has_same_sites(atom1,atom2):
                if atom1['uiso'] != atom2['uiso'] or atom1['aniso'] != atom2['aniso']:
                    raise _nc_core.NCBadInput('atoms labelled "%s" and "%s" have same cell '%(atom1['cif_label'],atom2['cif_label'])
                                              +'positions but different uiso information. This is not supported.')
                atoms_with_same_pos.append( atom2 )
            else:
                atoms_with_different_pos.append( atom2 )

        cif_labels = list( a['cif_label'] for a in atoms_with_same_pos )
        cif_labels_str = '"%s"'%('", "'.join(cl for cl in cif_labels))

        composition = list( ( a['occupancy'], a['element_marker'] ) for a in atoms_with_same_pos )
        occupancy = _math_fsum( fr for fr, elem in composition )
        assert occupancy > 0.0
        if occupancy != 1.0:
            composition = list( (fr/occupancy, elem) for fr, elem in composition )
            if occupancy > (1.0 + 1e-6) :
                m=f'Error: Too high total occupancy ({occupancy} which is >1) for atoms with CIF labels: {cif_labels_str}'
                if abs(int(occupancy+0.5)-occupancy)<1e-10:
                    m+='. A possible cause could be multiple entries in the CIF file using atomic positions that are actually symmetrically equivalent!'
                raise _nc_core.NCBadInput(m)
            occupancy = min( 1.0, occupancy )

        final_atom_info = dict( cif_labels = list_ro( cif_labels ),
                                equivalent_positions = list_ro( (c.x,c.y,c.z) for c in atoms_with_same_pos[0]['expanded_coords'] ),
                                composition = list_ro( sorted( composition ) ),
                                occupancy = occupancy,
                                uiso = atoms_with_same_pos[0]['uiso'],
                                aniso = atoms_with_same_pos[0]['aniso'] )
        final_atoms.append( dict_ro( final_atom_info )  )
        atoms_to_process = atoms_with_different_pos

    if merge_equiv:
        final_atoms = list_ro( dict_ro(a) for a in _impl_merge_atoms( final_atoms ) )

    if expected_tot_atom_in_orig_cell:
        n_actual = sum( len(a['equivalent_positions']) for a in final_atoms )
        def rd( x,y ):
            return abs(x-y)/(max(1e-300,abs(x)+abs(y)))
        if rd( n_actual, expected_tot_atom_in_orig_cell ) > 0.01:
            raise _nc_core.NCBadInput(f'Expanded number of atoms per cell ({n_actual}) is different '
                                      'from what is stated explicitly in the CIF data'
                                      f' ({expected_tot_atom_in_orig_cell}). Wrong spacegroup (or'
                                      ' choice/setting improperly specified)?')

    if refine_with_spglib:
        cellsg, final_atoms, warnings, msgs = _impl_refine_cell( cellsg, final_atoms )
        for w in warnings:
            _nc_common.warn(w)
        if not quiet:
            for m in msgs:
                _nc_common.print(m)
    else:
        _nc_common.warn('Structure refinement and verification with spglib was disabled due to explicit request.')

    if int(cellsg['spacegroup']['number']) == 1:
        #The following warning might be a bit obsolete now that we usually
        #refine cells, but we keep it here for now as an extra safeguard:
        _nc_common.warn('Space-group number 1 ("P1") detected in CIF input. A lot (but not all) of CIF data'
                       ' with this space group listed is in a non-conventional format not suitable for NCrystal. Be wary!')

    result['cellsg'] = dict_ro( cellsg )
    def _atom_sort( a ):
        #well-defined sort order for reproducible output
        return ( tuple(sorted(a.get('composition',[]))), len(a.get('equivalent_positions')),a.get('uiso',None),a.get('cif_labels',None), tuple(a.items()) )
    result['atoms'] = list_ro( dict_ro(a) for a in sorted(final_atoms,key = _atom_sort) )

    result['cifdescr'], result['actual_codid'],result['actual_mpid'],result['cif_chemformula'] =_extract_descr_from_cif( result['cif_raw'], result['cifsrc'], result['cifdata'] )
    return result, gemmi, struct
