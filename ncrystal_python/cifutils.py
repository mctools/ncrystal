"""
Utilities for converting local or online CIF files to NCMAT data.

Various utility functions are included as needed to support the related cmd-line
scripts (in particular ncrystal_cif2ncmat and ncrystal_onlinedb2ncmat), but they
can of course also be used directly via the Python API below by expert users.

"""


################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2022 NCrystal developers                                   ##
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


#FIXME: Doc string on top and functions below... __all__ as well.
#FIXME: check if error messages says "script"
#FIXME: E.g. cif_willend/FECO_B2.cif has some info about origins which we don't user in cif2descr
#FIXME: cif_uiso_temperature should be user-accessible in cmdline scripts
#FIXME: Add unit tests and examples for the cifutils python module! Also documentation online?
#FIXME: Do not raise SystemExit anywhere in a module!!

value_fallback_debye_temp = 300.0

import math
import functools
import pathlib
import os
import warnings
import numpy as np

#NB: cmdline scripts should have more helpful error messages for these:
import pymatgen
import pymatgen.symmetry.analyzer

#Main NCrystal module
import importlib
__name = 'NCrystal.cifutils' if __name__=='__main__' else __name__
NC = importlib.import_module( '..', __name )#importing __init__.py
nc_common = importlib.import_module( '.._common', __name )
nc_experimental = importlib.import_module( '..experimental', __name )
print = nc_common.print

class CifSource:
    """
    Generalised CIF source, either in the form of text data, a file path, or
    a database ID for either the Crystallography Open Database
    (https://www.crystallography.net/cod/) or the Materials Project
    (https://www.materialsproject.org/)
    """

    def __init__(self,*,codid=None,mpid=None,filepath=None,textdata=None):
        self.__codid,self.__mpid,self.__fp,self.__textdata = None,None,None,None
        if codid is not None:
            self.__codid = int(codid)
        if mpid is not None:
            self.__mpid = int(mpid)
        if filepath is not None:
            self.__fp = pathlib.Path(filepath)
        if textdata is not None:
            #textdata, decoding bytes to str assuming ascii/utf8:
            self.__textdata = ( textdata if isinstance(textdata,str) else textdata.decode() )
        if ( int(not self.__textdata is None)
             + int(not self.__codid is None)
             + int(not self.__mpid is None)
             + int(not self.__fp is None) ) != 1:
            raise ValueError('CifSource objects must be initialised with exactly one argument')

    @staticmethod
    def fromAnySrc( anysrc ):
        if isinstance(anysrc,CifSource):
            return anysrc
        if isinstance(anysrc,LoadedCIF) or isinstance(anysrc,CIFAnalyser):
            return anysrc.cifsrc
        if hasattr(anysrc,'__fspath__'):
            return CifSource( filepath = anysrc )
        if hasattr( anysrc, 'startswith' ):
            if anysrc.startswith('codid::') and anysrc[7:].isdigit():
                return CifSource(codid=int(anysrc[7:]))
            if anysrc.startswith('mpid::') and anysrc[6:].isdigit():
                return CifSource(mpid=int(anysrc[6:]))
        if isinstance( anysrc, str ):
            return CifSource(textdata=anysrc) if '\n' in anysrc else CifSource(filepath=pathlib.Path(anysrc))
        if isinstance( anysrc, bytes ):
            return CifSource(textdata=anysrc.decode()) if b'\n' in anysrc else CifSource(filepath=pathlib.Path(anysrc))
        raise ValueError('Could not detect CIF source type')

    @property
    def codid( self ): return self.__codid
    @property
    def mpid( self ): return self.__mpid
    @property
    def textdata( self ): return self.__textdata
    @property
    def filepath( self ): return self.__fp

def getMaterialsProjectAPIKEY():
    import os
    apikey=os.environ.get('MATERIALSPROJECT_USER_API_KEY',None)
    if not apikey:
        raise RuntimeError('ERROR: Missing API key for materialsproject.org access. To fix, '
                           +'please make sure the environment variable MATERIALSPROJECT_USER_API_KEY'
                           +' contains your personal access key that you see on'
                           +' https://www.materialsproject.org/dashboard (after logging in).')
    return apikey

def _classifySG(sgno):
    assert 1<=sgno<=230
    l=[(195,'cubic'),(168,'hexagonal'),(143,'trigonal'),
       (75,'tetragonal'),(16,'orthorombic'),(3,'monoclinic'),(1,'triclinic')]
    for thr,nme in l:
        if sgno>=thr:
            return nme
    assert False

def _cif2descr(cif):
    #cif is cif section with publication info

    #FIXME also support:
    # ('_citation_journal_full', ['Chemistry of Materials']),
    #             ('_citation_year', ['2019']),
    #             ('_citation_journal_volume', ['31']),
    #             ('_citation_page_first', ['7203']),
    #             ('_citation_page_last', ['7211']),
    #             ('_citation_journal_id_ASTM', ['CMATEX']),
    #             ('_citation_author_citation_id',

    l=[]
    def extract(key,*,expectlist = False):
        if not key in cif:
            return [] if expectlist else ''
        s=cif[key]
        if not expectlist:
            s=s.strip()
            return '' if s=='?' else s
        else:
            l=list( ('' if e.strip()=='?' else e.strip()) for e in s )
            return list( e for e in l if e )

    title = extract('_publ_section_title')
    authors = extract('_publ_author_name',expectlist=True)
    year = extract('_journal_year')
    doi = extract('_journal_paper_doi')
    codid = extract('_cod_database_code')
    audit_creation_method = extract('_audit_creation_method')
    audit_creation_date = extract('_audit_creation_date')
    audit_author_name = extract('_audit_author_name')

    if authors:
        if len(authors)==1:
            authors = authors[0]
        elif len(authors)==2:
            authors = f'{authors[0]} and {authors[1]}'
        else:
            authors = f'{authors[0]}, et al.'
    if title:
        l.append(f'"{title}"')
    if authors and year:
       l.append(f'{authors} ({year})')
    if doi:
        l.append(f'DOI: https://dx.doi.org/{doi}')
    if codid:
        l += [ f'Crystallography Open Database entry {codid}',
               f'https://www.crystallography.net/cod/{codid}.html' ]
    if audit_creation_method:
        l.append(f'Creation method: {audit_creation_method}')
    if audit_creation_date:
        l.append(f'Creation date: {audit_creation_date}')
    if audit_author_name:
        l.append(f'Created by: {audit_author_name}')

    return l

class LoadedCIF:
    """Information loaded from a CIF source.

    For reference it also contains information about the original CIF source, as
    well as a list of warnings emitted while loading (these will have been
    silenced if quiet=True, otherwise they will have been raised as usual to the
    calling code).

    """


    def __init__( self, cifsrc, *, quiet ):
        self.__cifsrc = CifSource.fromAnySrc( cifsrc )
        with nc_common.WarningSpy() as warnlist:
            d = _loadCIF_actual( cifsrc, quiet )
        assert len(d)==5
        self.__w,self.__s,self.__d = warnlist,d['structure'],d['descr']
        self.__cr,self.__cs,self.__cp = d['cif_raw'],d['cif_block_sites'],d['cif_block_publ']

    @staticmethod
    def fromAnySrc( anysrc, quiet = False ):
        if isinstance(anysrc,LoadedCIF):
            return anysrc
        if isinstance(anysrc,CIFAnalyser):
            return anysrc.loaded_cif
        return LoadedCIF(anysrc,quiet=quiet)

    @property
    def cifsrc( self ):
        return self.__cifsrc

    @property
    def cif_raw( self ):
        return self.__cr

    @property
    def cif_publ( self ):
        return self.__cp

    @property
    def cif_sites( self ):
        return self.__cs

    @property
    def structure( self ):
        return self.__s

    @property
    def descr( self ):
        #FIXME: To CIFAnalyser?
        return self.__d

    @property
    def warnings( self ):
        return self.__w

    def _append_warnlist(self,warnlist):
        #FIXME: Get rid of this mutability!!
        for wmsg,wcat in warnlist:
            self.__w.append( (wmsg,wcat) )

def _rawLoad_CifAndStructure(cifdata,quiet):
    import pymatgen.core.structure
    import pymatgen.io.cif

    #First just load:
    if not quiet:
        print("Attempting to load CIF data (using pymatgen module).")
    cif_obj = pymatgen.io.cif.CifFile.from_string( cifdata )
    cif_raw = cif_obj.data

    #Dig out publ + site blocks from cif:
    if not quiet:
        print("Attempting to locate relevant sections of CIF data.")
    cifblocks_publ,cifblocks_sites = [],[]
    for blockname,blockdata in cif_raw.items():
        has_key_starting_with = lambda x : any(k.startswith(x) for k,v in blockdata.data.items())
        is_publ = ( has_key_starting_with('_publ_')
                    or has_key_starting_with('_journal_')
                    or has_key_starting_with('_audit_')
                    or has_key_starting_with('_cod_database_code') )
        is_sites = has_key_starting_with('_atom_site_')
        if is_publ:
            cifblocks_publ.append( blockdata.data )
        if is_sites:
            cifblocks_sites.append( blockdata.data )

    if len(cifblocks_sites) > 1:
        raise ValueError('CIF file contains multiple sections with atom site info (multiphase CIF files are not supported).')
    if len(cifblocks_publ) > 1:
        warnings.warn('CIF file contains multiple sections with publication info. All but the first will be ignored.')
    cif_sites = cifblocks_sites[0] if cifblocks_sites else None
    cif_publ = cifblocks_publ[0] if cifblocks_publ else None

    #Extract the structure:
    if not quiet:
        print("Attempting to retrieve crystal structure from CIF data (using pymatgen module).")
    structure = pymatgen.core.structure.Structure.from_str( cifdata, fmt="cif" )

    #Make sure it is a conventional, rather than primitive, cell:
    sga = pymatgen.symmetry.analyzer.SpacegroupAnalyzer(structure)
    structure_conventional = sga.get_conventional_standard_structure()

    #Done:
    return cif_raw, cif_sites, cif_publ, structure_conventional

def _loadCIF_actual( cifsrc, quiet ):
    cifsrc = CifSource.fromAnySrc( cifsrc )
    textdata = cifsrc.textdata
    if textdata is None and cifsrc.filepath is not None:
        if not quiet:
            print(f"Loading data from file {cifsrc.filepath}")
        textdata = cifsrc.filepath.read_text()
    if textdata is not None:
        cif_raw,cif_sites,cif_publ,structure = _rawLoad_CifAndStructure(textdata,quiet)
        if cifsrc.filepath is not None:
            descr = [f'File: {cifsrc.filepath.name}']
        else:
            descr = [f'Anonymous CIF data']
        cifdescr = _cif2descr(cif_publ) if cif_publ else None
        if cifdescr:
            descr.append('')
            descr += cifdescr

    elif cifsrc.codid is not None:
        if not quiet:
            print(f"Querying the Crystallography Open Database for entry {cifsrc.codid}")
        import requests
        r=requests.get("https://www.crystallography.net/cod/%i.cif"%(cifsrc.codid))
        r.raise_for_status()#throw exception in case of e.g. 404
        cif_raw,cif_sites,cif_publ,structure = _rawLoad_CifAndStructure(r.text,quiet)
        descr = _cif2descr(cif_publ) if cif_publ else ['']
        url = f'https://www.crystallography.net/cod/{cifsrc.codid}.html'
        if not any(url in e for e in descr):
            descr += [ f'Crystallography Open Database entry {cifsrc.codid}',
                       f'https://www.crystallography.net/cod/{cifsrc.codid}.html' ]
    else:
        assert cifsrc.mpid is not None
        from pymatgen.ext.matproj import MPRester
        if not quiet:
            print(f"Querying the Materials Project for entry {cifsrc.mpid}")
        with MPRester(getMaterialsProjectAPIKEY()) as m:
            structure = m.get_structure_by_material_id("mp-%i"%cifsrc.mpid,
                                                       conventional_unit_cell=True)
        #No actual CIF data provided (at least not in the API we use):
        cif_raw = {}
        cif_publ = {}
        cif_sites = {}
        descr=['The Materials Project',
               f'https://www.materialsproject.org/materials/mp-{cifsrc.mpid}']

    #always return in same format:
    return dict( structure = structure, descr = descr,
                 cif_raw=cif_raw,
                 cif_block_publ=cif_publ,
                 cif_block_sites=cif_sites )

def extract_uiso( cif ):
    """Try to dig out atomic displacements (Uiso) information from the cif data.
    """
    #FIXME: Emit warnings as appropriate! + do this extraction in the CIFAnalyser

    cif = LoadedCIF.fromAnySrc(cif)
    c = cif.cif_sites

    if not c:
        #e.g. the case for materialsprojects input.
        return {}
    if not '_atom_site_label' in c:
        nc_common.warn('Could not extract uiso values from CIF data: No _atom_site_label entries in data!')
        return {}

    if not('_atom_site_U_iso_or_equiv' in c or '_atom_site_B_iso_or_equiv' in c):
        #This return statement will trigger for most files, a warning here would be too verbose.
        return {}
    albls = c['_atom_site_label']
    if '_atom_site_type_symbol' in c and len(c['_atom_site_type_symbol'])==len(albls):
        albls = c['_atom_site_type_symbol']
    elif '_atom_type_symbol' in c and len(c['_atom_type_symbol'])==len(albls):
        albls = c['_atom_type_symbol']

    auiso = c.get('_atom_site_U_iso_or_equiv',None)
    expected_adptype = 'Uiso'
    auiso_factor = 1.0
    if auiso is None:
        auiso = c['_atom_site_B_iso_or_equiv']
        auiso_factor = 1.0 / ( 8 * math.pi**2 )
        expected_adptype = 'Biso'
    if len(albls)!=len(auiso):
        nc_common.warn('Could not extract uiso values from CIF data: Number of uiso/biso fields does not match number of sites.')
        return {}

    adptypes = c.get('_atom_site_adp_type',None)
    if adptypes is None:
        adptypes = c.get('_atom_site_thermal_displace_type',None)
    if adptypes is not None:
        if len(albls)!=len(adptypes):
            nc_common.warn('Could not extract uiso values from CIF data: Number of _atom_site_thermal_displace_type fields does not match number of sites.')
            return {}
        if not all(e==expected_adptype for e in adptypes):
            nc_common.warn('Could not extract uiso values from CIF data: Not all sites has same type of uiso/biso information.')
            return {}

    #Ok, all labels have Uiso/Biso. Proceed to decode them into elem2uiso map,
    #with various checks along the way.

    def extractElementName(lbl):
        if lbl.endswith('+') or lbl.endswith('-'):
            #Remove oxidation state, e.g. "Fe2+" -> "Fe"
            l = lbl[:-1]
            while l and l[-1].isdigit():
                l = l[:-1]
            lbl = l
        for e in nc_experimental.knownElementNames():
            #Case insensitive cmp, support e.g. "AL" instead of "Al":
            if lbl.lower()==e.lower():
                return e
        return None

    #Now check all the names are proper element names and with unique values
    #(otherwise it is too complicated to support):
    d={}
    lbls_badnames = set()
    for lbl,uiso in zip(albls,auiso):
        #Map to proper element name, or warn+skip if not possible:
        if lbl in lbls_badnames:
            continue
        elementName = extractElementName(lbl)
        if elementName is None:
            nc_common.warn(f'Ignoring Uiso/Biso info for label {lbl} since it is not a recognised element name')
            lbls_badnames.add( lbl )
            continue
        if not elementName in d:
            d[elementName] = set()
        d[elementName].add( uiso )

    elem2uiso = {}

    def decode_value(s):
        #support values with uncertainties like '0.080(10)' (simply discards the
        #uncertainty).
        try:
            if s.count('(')==1 and s.count(')')==1:
                s=s.split('(',1)[0]
            return float(s)
        except ValueError:
            pass

    for elem,uiso_values in d.items():
        if len(uiso_values) > 1:
            nc_common.warn( f'Ignoring Uiso/Biso info for element "{elem}" since more than one value were provided: '
                            + '"%s"'%('", "'.join(sorted(uiso_values))))
        else:
            val = decode_value( uiso_values.pop() )
            if val is None or not (0.0<val<1e5):
                nc_common.warn( f'Ignoring Uiso/Biso info for element "%s" since value could not be decoded or is out of range'%elem )
            else:
                elem2uiso[ elem ] = val * auiso_factor

    return elem2uiso

#dyninfos can either be a dyninfo object (originating in an Info object), or an
#instance of either AtomMSD or AtomDebyeTemp.

class AtomMSD:

    def __init__(self,*,name,msd,temperature):
        """Takes name of atom (e.g. "Al"), mean-squared-displacement (msd), and
        temperature at which the msd is valid. Units are respectively
        squared-Angstrom and kelvin.
        """
        assert float( msd ) > 0.0
        assert float( temperature ) > 0.0
        self.__name,self.__m,self.__t = str(name),float(msd),float(temperature)

    @property
    def name(self):
        """Name of atom."""
        return self.__name

    @property
    def msd(self):
        """Mean-squared-displacement value (Aa^2)"""
        return self.__m

    @property
    def temperature(self):
        """Temperature value  (kelvin) associated with the msd value"""
        return self.__t

class AtomDebyeTemp:

    def __init__(self,*,name,debye_temperature):
        """Takes name of atom (e.g. "Al"), and debye temperature value in kelvin."""
        assert float( debye_temperature ) > 0.0
        self.__name,self.__dt = str(name),float(debye_temperature)

    @property
    def name(self):
        """Name of atom."""
        return self.__name

    @property
    def debye_temperature(self):
        """Temperature value  (kelvin) associated with the msd value"""
        return self.__dt

def produceNCMAT( cif, dynamics, rawformat, atomdb, cif_uiso_temperature = None, display_labels = None ):

    #FIXME: Get rid of this deprecated function in favour of CIFAnalyser etc.

    cif = LoadedCIF.fromAnySrc(cif)

    if not dynamics:
        dynamics = dict(comments=[],dyninfos=[])

    #In case source had deuterium directly and we have to remap it via the
    #atomdb parameter, we must do a trick and assign the element another
    #name. This is because the @ATOMDB section can not reassign isotope names
    #another meaning according to the NCMAT format.
    remap_D = False
    if atomdb:
        _=atomdb.replace(':',' ').split()
        assert len(_)>=3 and _[1]=='is'
        remap_D = _[0]=='D'
        if remap_D:
            _[0] = 'X99'
        atomdb = ' '.join(_)

    elemNameRemap = lambda en : "X99" if (remap_D and en=="D") else en

    extra_atomdb = []
    if display_labels:
        if remap_D and any(v=='X99' for k,v in display_labels.items()):
            raise ValueError('Can not use label X99 when also trying to remap special element "D" (deuterium)')
        _orig_elemNameRemap = elemNameRemap
        elemNameRemap = lambda en : display_labels[en] if (en in display_labels) else _orig_elemNameRemap(en)
        extra_atomdb = list( '%s is %s'%(v,k) for k,v in display_labels.items())

    #fixme: no need for this extra warning tracking?
    emitted_warnings = []
    def produceWarning(w):
        nc_common.warn(w)
        cif._append_warnlist( [nc_common.WarningSpy._fmtwarning( str(w), NC.NCrystalUserWarning )] )

    def elemName(site):
        species = site.species
        _=site.as_dict()['species']
        if not species.is_element or not len(species.elements)==1 or not len(_)==1:
            raise ValueError('NCrystal CifUtils currently only supports a single well defined element at each site')
        occu=_[0]['occu']
        if occu > 0.99:
            produceWarning( f"Encountered species with almost-but-not-quite unit occupancy ({occu}). Treating as 1.0" )
            occu = 1.0
        if occu!=1.0:
            raise ValueError(f'Error: This script only supports species with unit occupancies (encountered occupancy {occu})')
        return _[0]['element']

    elems = []
    formula_elem2count={}
    def fmt(x):
        s='%.14g'%x
        if s.startswith('0.'):
            s=s[1:]
        return s
    def apfmt(x):
        while x<0.0:
            x+=1
        while x>=1.0:
            x-=1.0
        return (fmt if rawformat else nc_common.prettyFmtValue)(x)
    out_atompos_lines = ''
    for s in cif.structure.sites:
        c=s.frac_coords
        en = elemName(s)
        elems += [en]
        out_atompos_lines += f'  {elemNameRemap(en)} {apfmt(c[0])} {apfmt(c[1])} {apfmt(c[2])}\n'
        formula_elem2count[en] = formula_elem2count.get(en, 0) + 1
    elems_unique = sorted(list(set(elems)))

    a,b,c = cif.structure.lattice.abc
    alpha,beta,gamma = cif.structure.lattice.angles
    sgsymb,sgno = cif.structure.get_space_group_info()
    descr = ''.join('#    %s\n'%e for e in cif.descr)

    if 195 <= sgno <= 230:
        cell=f'@CELL\n  cubic {fmt(a)}'
        assert fmt(a)==fmt(b) and fmt(a)==fmt(c), "cubic cell must have all lattice parameters identical: a==b==c"
        assert ( fmt(alpha)=='90' and fmt(beta)=='90' and fmt(gamma)=='90' ), "cubic cell must have all angles exactly 90 degrees"
    else:
        cell=f'@CELL\n  lengths {fmt(a)} {fmt(b) if fmt(b)!=fmt(a) else "!!"} {fmt(c)}\n  angles {fmt(alpha)} {fmt(beta)} {fmt(gamma)}'

    #
    elem2dyninfo={}
    elems_with_fallback_dyninfo = []
    for en in elems_unique:
        #dynamics = dict(comments=[],dyninfos=[])
        for di in dynamics['dyninfos']:
            if hasattr(di,'atomData'):
                name = di.atomData.displayLabel()
            else:
                name = di.name#AtomMSD or AtomDebyeTemp
            if name == en:
                elem2dyninfo[en] = di
                break
        if not en in elem2dyninfo:
            #Check if perhaps this is just an isotope issue (e.g. "D" instead of
            #"H"), and search for element with same Z:
            en_atomData = NC.atomDB(en)
            for di in dynamics['dyninfos']:
                if hasattr(di,'atomData') and di.atomData.isElement() and di.atomData.Z()==en_atomData.Z():
                    elem2dyninfo[en] = di
                    break

        if not en in elem2dyninfo:
            elems_with_fallback_dyninfo.append( en )

    if dynamics['dyninfos'] and not elem2dyninfo:
        raise SystemExit('Error: was not able to transfer any dynamics from source.')#fixme not systemexit

    formula=cif.structure.composition.reduced_formula
    if remap_D:
        formula+='(WARNING:D was not changed automatically)'
    title=f'{formula} ({_classifySG(sgno)}, SG-{sgno} / {sgsymb})'
    orig_comments=''

    def calc_and_format_debyetemp_from_msd(msd,temperature,mass):
        dtval = NC.debyeTempFromIsotropicMSD( msd=msd,
                                              temperature=temperature,
                                              mass = mass )
        return f'debye_temp {dtval:.14g} # From msd={msd:.14g}Aa^2 @ T={temperature:.14g}K (M={mass:.14g}u) \n'

    elems_with_cif_debyetemp = {}
    if elems_with_fallback_dyninfo:
        #Try to dig out Uiso from cif source
        uiso_info = extract_uiso(cif)
        if cif_uiso_temperature is None:
            #FIXME: Warn (and explain about flag???)
            cif_uiso_temperature = 293.15
        for en in elems_with_fallback_dyninfo[:]:
            if en in uiso_info:
                en_atomData = NC.atomDB(en)
                msd = uiso_info[en]
                en_mass = en_atomData.averageMassAMU()
                elems_with_cif_debyetemp[en] = calc_and_format_debyetemp_from_msd(msd,cif_uiso_temperature,en_mass)
#                en_debyeTemp = NC.debyeTempFromIsotropicMSD( msd=msd,
#                                                             temperature=cif_uiso_temperature,
#                                                             mass=en_mass )
#                elems_with_cif_debyetemp[en] = en_debyeTemp
        elems_with_fallback_dyninfo = [e for e in elems_with_fallback_dyninfo if not e in elems_with_cif_debyetemp]

    if elems_with_fallback_dyninfo:
        orig_comments += '\n# WARNING: Dummy Debye temperature values were\n'
        orig_comments += '#          inserted below for at least 1 element!\n'
        orig_comments += '#'
    if dynamics.get('comments',None):
        orig_comments += '\n# NOTICE: For reference comments in the file used as source\n'
        orig_comments += '#         for the @DYNINFO sections are repeated below.\n'
        orig_comments += '#         (please extract relevant sections and merge with\n'
        orig_comments += '#          comments above in order to finalise file).\n'
        orig_comments += '#\n'
        for e in dynamics['comments']:
            orig_comments += f'#   ---> {e}\n'
        orig_comments+= '#'
    atomdb = '\n@ATOMDB\n  %s'%(' '.join(atomdb.replace(':',' ').split())) if atomdb else ''
    if extra_atomdb and not atomdb:
        atomdb = '\n@ATOMDB'
    for e in extra_atomdb:
        atomdb += f'\n  {e}'

    warnstr=''
    if cif.warnings:
        warnstr+='\n# Notice: The following WARNINGS were emitted when loading the CIF data with pymatgen. They were:\n#'
        nwmax=10
        for wmsg,wcat in cif.warnings:
            wmsg=list(e.strip() for e in wmsg.splitlines() if e.strip())
            for i,m in enumerate(wmsg):
                if i==nwmax:
                    break
                pref='  %s :'%wcat if i==0 else ' '*(len(wcat)+4)
                if i+1==nwmax and len(wmsg)>nwmax:
                    m='<%i lines of output hidden>'%(len(wmsg)-nwmax)
                warnstr+='\n# %s %s'%(pref,m)
        warnstr+='\n#'

    out = f"""NCMAT v5
#
# {title}
#
# Structure converted (with NCrystal.cifutils module) from:
#
{descr}#
# IMPORTANT NOTICE: This is an automatic conversion which has not
# been verified!  In particular the @DYNINFO sections might need
# post-editing. Before distributing this file to other people,
# please review this, amend the comments here to document anything
# done, and remove this notice.
#{warnstr}{orig_comments}{atomdb}
{cell}
@SPACEGROUP
  {sgno}
@ATOMPOSITIONS
"""
    out += out_atompos_lines
    #reduce chemical formula:
    gcd = functools.reduce(math.gcd,[v for k,v in formula_elem2count.items()])
    for en in formula_elem2count.keys():
        formula_elem2count[en] //= gcd
    chemform = ''.join(f'{en}{count if count>1 else ""}' for en,count in sorted(formula_elem2count.items()))

    assert len(elem2dyninfo)+len(elems_with_fallback_dyninfo)+len(elems_with_cif_debyetemp) == len(elems_unique)

    for en in elems_unique:
        di = elem2dyninfo.get(en,None)
        fr_str = f'{formula_elem2count[en]}/{len(elems)//gcd}' if len(elems_unique)>1 else '1'
        out += f'@DYNINFO\n  element  {elemNameRemap(en)}\n  fraction {fr_str}\n'
        if di is None:
            cifdt = elems_with_cif_debyetemp.get(en,None)
            if cifdt is not None:
                out += f'  type     vdosdebye\n  '+cifdt#debye_temp {cifdt:g} # Based on Uiso in CIF input (assuming T={cif_uiso_temperature:g}K)\n'
            else:
                out += f'  type     vdosdebye\n  debye_temp {value_fallback_debye_temp:g} # WARNING: Dummy fall-back value added here!!\n'
        elif isinstance(di,AtomMSD):
            atom = NC.atomDB( di.name )
            out += f'  type     vdosdebye\n  '
            out += calc_and_format_debyetemp_from_msd(di.msd,di.temperature,atom.averageMassAMU())
#            mass=atom.averageMassAMU()
#            dtval = NC.debyeTempFromIsotropicMSD( msd=di.msd,
#                                                  temperature=di.temperature,
#                                                  mass = mass )
#            out += f'  type     vdosdebye\n  debye_temp {dtval:.14g} # From msd={di.msd:.14g}Aa^2 @ T={di.temperature:.14g}K (M={mass:.14g}u) \n'
        elif isinstance(di,AtomDebyeTemp):
            out += f'  type     vdosdebye\n  debye_temp {di.debye_temperature:.14g}\n'
        elif isinstance(di,NC.Info.DI_VDOSDebye):
            out += f'  type     vdosdebye\n  debye_temp {di.debyeTemperature():.14g}\n'
        else:
            assert isinstance(di,NC.Info.DI_VDOS)
            out += '  type     vdos\n'
            out += NC.formatVectorForNCMAT('vdos_egrid',di.vdosOrigEgrid())
            out += NC.formatVectorForNCMAT('vdos_density',di.vdosOrigDensity())

    l = list(e.rstrip() for e in out.splitlines())
    #if emitted_warnings:
    #    l2 = []
    #    l2.append( l[0] )#NCMAT vx line
    #    for w in emitted_warnings:
    #        l2.append(f'# WARNING: {w}')
    #    for e in l[1:]:
    #        l2.append( e )
    #    l = l2
    l.append('')
    return '\n'.join(l)

_keepalive_dynamics=[]
def extractDynamics(cfgstr_or_info):
    #returns { 'comments':strlist, 'dyninfos':dyninfos }
    if hasattr(cfgstr_or_info,'dyninfos'):
        cfgstr, info = None, cfgstr_or_info
    else:
        cfgstr = cfgstr_or_info
        info = NC.createInfo(cfgstr)

    #DynInfos:
    _keepalive_dynamics.append(info)
    dyninfos=[ di for di in info.dyninfos
               if ( isinstance(di,NC.Info.DI_VDOS) or isinstance(di,NC.Info.DI_VDOSDebye) ) ]

    #Comments:
    comments = []
    for l in (NC.createTextData(cfgstr) if cfgstr else []):
        l=l.strip()
        if l.startswith('#'):
            comments.append(l[1:])
        if l.startswith('@'):
            break
    #remove any common leading spaces:
    while ( comments
            and any(e.startswith(' ') for e in comments)
            and all((e.startswith(' ') or not e) for e in comments) ):
        for i in range(len(comments)):
            comments[i] = comments[i][1:]
    return dict( comments = comments, dyninfos=dyninfos )

def lookupAndProduce( cifsrc, *,dynamics=None,rawformat=False,atomdb=None,quiet=False,cif_uiso_temperature=None,display_labels=None):
    #source must be a file-path or a tuple with (onlinedatabasename,databaseid).
    #Supported onlinedatabasename's are so far 'COD' and 'MP'.

    #FIXME: Get rid of this deprecated function in favour of CIFAnalyser etc.


    autofn=[]
    cifsrc = CifSource.fromAnySrc( cifsrc )
    if cifsrc.codid is not None:
        autofn.append(f'cod{cifsrc.codid}')
    elif cifsrc.mpid is not None:
        autofn.append(f'mp{cifsrc.mpid}')
    if atomdb:
        _adb=atomdb.replace(':',' ').split()
        if not len(_adb)>=3 or _adb[1]!='is':
            raise SystemExit(f'Error: Invalid ATOMDB format encountered: "{atomdb}" (supports "X is ..." syntax only)')#fixme not systemexit
        descr[-1] += f' [with {_adb[0]}->{" ".join(_adb[2:])}]'
    cif = LoadedCIF.fromAnySrc(cifsrc,quiet=quiet)

    autofn.insert(0,f'sg{cif.structure.get_space_group_info()[1]}')
    autofn.insert(0,cif.structure.composition.reduced_formula)
    autofn.insert(0,'autogen')
    out = produceNCMAT( cif,
                        dynamics = dynamics,
                        rawformat=rawformat,
                        atomdb=atomdb,
                        cif_uiso_temperature = cif_uiso_temperature,
                        display_labels = display_labels
                       )
    return '_'.join(autofn)+'.ncmat', out

def _actual_cif_analyse( loadedcif, *, quiet ):

    cif = loadedcif

    # Decode sites and sort into species with different atom types and occupancy:
    uiso = extract_uiso(cif)

    species_info = {}
    for i,site in enumerate(cif.structure.sites):
        coords = site.frac_coords
        species = site.species
        site_as_dict_species=site.as_dict()['species']
        if not species.is_element or not len(species.elements)==1 or not len(site_as_dict_species)==1:
            #NB: There is no reason we can't support this!! We just need an
            #example of a CIF file with non-element species to test it on!
            raise ValueError('NCrystal CifUtils currently only supports a single well defined element at each site.'
                             +' Please send your CIF file to the NCrystal developers if you wish us to try and remedy this.')

        d=site_as_dict_species[0]#we checked there is only one entry just above
        elementName = str(d['element'])
        site_species_info = dict( element_name = elementName,#we checked "species.is_element" just above!
                                  occupancy = float(d.get('occu',1.0))
                                  #We could split on oxidation state as well, but we don't for now:
                                  #oxidation_state = float(d.get('oxidation_state',None)),
                                 )
        site_species_info['uiso'] = uiso.get(elementName,None)
        #Combine coords of any entries with same site_species_info
        key = tuple( sorted( site_species_info.items() ) )
        if not key in species_info:
            species_info[key] = dict( info = site_species_info, coords = [] )
        species_info[key]['coords'] += [ ( float(coords[0]), float(coords[1]), float(coords[2]) ) ]

    atoms = []
    for i,(_,e) in enumerate(sorted(species_info.items())):
        atoms.append( dict( idx = i,
                            element = e['info']['element_name'],
                            site_occupancy = e['info']['occupancy'],
                            uiso = e['info']['uiso'],
                            atompositions = e['coords'] ) )

    a,b,c = (float(e) for e in cif.structure.lattice.abc)
    alpha,beta,gamma = (float(e) for e in cif.structure.lattice.angles)

    if (not a>0.0) or not (b>0.0) or not (c>0.0):
        raise RuntimeError(f'Structure loaded from cif has invalid lattice parameters: {a}, {b}, {c}')
    if (not 0.0<alpha<180.0) or not (0.0<beta<180.0) or not (0.0<gamma<180.0):
        raise RuntimeError(f'Structure loaded from cif has invalid lattice angles: {alpha}, {beta}, {gamma}')
    if ( alpha<3.1416 and beta<3.1416 and gamma<3.1416 ):
        raise RuntimeError(f'Structure loaded from cif has invalid lattice angles (might be in radians rather than degrees): {alpha}, {beta}, {gamma}')

    sgsymb,sgno = cif.structure.get_space_group_info()
    if not sgno or not 1<=int(sgno)<=230:
        if sgno:
            nc_common.warn(f'Ignoring invalid space group number: {sgno}')
        sgno = None

    if not sgsymb or not str(sgsymb).strip():
        sgsymb = None
    else:
        sgsymb = str(sgsymb).strip()

    #A simple sanity check
    if sgno and 195 <= sgno <= 230:
        #cubic:
        if not a==b or not a==c:
            raise RuntimeError("cubic cell loaded from CIF file does not have all lattice parameters identical (a==b==c)")
        if not alpha==90.0 or not beta==90.0 or not gamma==90.0:
            raise RuntimeError("cubic cell loaded from CIF file does not have all lattice angles equal to 90 degrees")

    structure = dict( spacegroup = sgno,
                      sgsymb = sgsymb,
                      a = float(a), b = float(b), c = float(c),
                      alpha = float(alpha), beta = float(beta), gamma = float(gamma) )

    return atoms, structure

class CIFAnalyser:

    def __init__( self, cifsrc, quiet = False ):
        self._c = LoadedCIF.fromAnySrc(cifsrc,quiet=quiet)
        with nc_common.WarningSpy() as warnlist:
            self.__atoms,self.__structure = _actual_cif_analyse( self._c, quiet=quiet )
        self.__analysis_warnings = warnlist
        #Make immutable:
        import types
        self.__structure = types.MappingProxyType(self.__structure)#read-only view
        self.__atoms = tuple( self.__atoms )

    def suggest_file_basename():
        pass#fixme

    @property
    def loaded_cif( self ):
        return self._c

    @property
    def cifsrc( self ):
        return self._c.cifsrc

    @property
    def warnings( self ):
        """Warnings emitted during analysis of the CIF data (this is in addition to the warnings emitted when loading the CIF data)"""
        return self.__analysis_warnings

    @property
    def atoms( self ):
        """Returns information about atoms in a format like:
            [
             {'idx':0,'element':'Al','site_occupancy':1.0,'atompositions':[<positions>],'uiso':None},
             {'idx':1,'element':'Al','site_occupancy':0.2,'atompositions':[<positions>],'uiso':None},
             {'idx':2,'element':'O','site_occupancy':1.0,'atompositions':[<positions>],'uiso':0.025},
             ...
            ]
        """
        return self.__atoms

    @property
    def structure( self ):
        """Returns information about symmetries and unit cell in a format like:
           { 'spacegroup' : 225, # 1-230 or None
             'sgsymb' :  "Fm-3m", # None if absent
             'a' : 4.012,
             'b' : 4.012,
             'c' : 4.012,
             'alpha' : 90,
             'beta' : 90,
             'gamma' : 90
           }
        """
        #Fixme: should sgsymb belong in the description instead?
        return self.__structure

    @property
    def description( self ):
        """Description deduced from metadata fields and cif src (suitable for NCMAT comments)"""
        pass#fixme. Remove this?


    def create_ncmat( self, *args, **kwargs ):
        """Return NCMAT data corresponding to the currently held crystal
        parameters. All parameters are forwarded (same as calling
        create_ncmat_composer(*args,**kwargs).create_ncmat())

        """
        c = self.create_ncmat_composer(*args,**kwargs)
        return c.create_ncmat()

    def create_ncmat_composer( self, *, uiso_temperature = None, skip_dyninfo = False ):
        """Setup and return an NCMATComposer object, based on the currently held
        crystal parameters.  If you wish to take advantage of any uiso/biso
        values found in the CIF data, you must provide the uiso_temperature
        parameter value (kelvin), to indicate the temperature at which these are valid.

        If skip_dyninfo=True, no ncmat.set_dyninfo_... calls will be performed,
        presumably because the caller will perform such calls subsequently.
        """

        NCMATComposer = nc_experimental.NCMATComposer
        ncmat = NCMATComposer()
        st = self.__structure
        ncmat.set_cellsg( a=st['a'], b=st['b'], c=st['c'],
                          alpha=st['alpha'], beta=st['beta'], gamma=st['gamma'],
                          spacegroup=st['spacegroup'] )

        element_list = [ a['element'] for a in self.__atoms ]
        assert all(nc_experimental.isKnownElement(e) for e in element_list), "non-pure-elements not supported yet"

        all_atompos = []
        for a in self.__atoms:
            occu = a.get('site_occupancy',1.0)
            for c in a['atompositions']:
                all_atompos.append( ( a['element'], c[0], c[1], c[2], occu ) )

        ncmat.set_atompos( all_atompos )

        for a in ([] if skip_dyninfo else self.__atoms):
            uiso = a.get('uiso',None)
            if uiso is not None:
                assert uiso>0.0
                if uiso_temperature is None:
                    pass#fixme: Add warning somehow about unused uiso info!
                else:
                    ncmat.set_dyninfo_msd(a['element'],msd=uiso, temperature=uiso_temperature)
                    continue
            ncmat.set_dyninfo_vdosdebye(a['element'],300.0,comment='WARNING: Dummy fall-back value added here!!')

        return ncmat
