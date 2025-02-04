
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

Internal implementation details for NCMAT utilities in ncmat.py

"""

__all__ = []
from . import core as _nc_core
from . import _common as _nc_common
import math
import copy

_magic_two_space=b'\xc3\x98@'.decode('utf8')

class NCMATComposerImpl:

    def set_plotlabel( self, lbl ):
        if lbl != self.__params.get('plotlabel'):
            self.__dirty()
            self.__params['plotlabel'] = lbl

    @property
    def plotlabel( self ):
        return self.__params.get('plotlabel')

    def add_secondary_phase(self, fraction, cfgstr, normalise ):
        if normalise:
            from .cfgstr import normaliseCfg
            c = normaliseCfg( cfgstr )
        else:
            c = ' '.join(str(cfgstr).split())
        if not c:
            raise _nc_core.NCBadInput('empty cfg-string used to define secondary phase')
        if '#' in c:
            raise _nc_core.NCBadInput('cfg-strings used to define secondary phases in NCMAT data can not contain "#" characters')
        if len(c.splitlines())>1 or '\n' in c or '\r' in c:
            raise _nc_core.NCBadInput('cfg-strings used to define secondary phases in NCMAT data can not contain newline characters')
        if not c.isascii():
            raise _nc_core.NCBadInput('cfg-strings used to define secondary phases in NCMAT data must be pure ASCII.')
        p = self.__params.get('secondary_phases',[])[:]
        p.append( ( float(fraction), c)  )
        ftot = math.fsum([_f for _f,_ in p])
        if ftot >= 1.0-1e-14:
            raise _nc_core.NCBadInput('Secondary phases have too high a fraction (a total of {ftot*100.0}%%)')
        self.__dirty()
        self.__params['secondary_phases'] = p

    def add_hard_sphere_sans_model( self, sphere_radius ):
        assert sphere_radius > 0.0
        self.__dirty()
        self.__params['custom_hardspheresans'] = float( sphere_radius )

    def add_raw_content( self,  content ):
        if not isinstance(content,str):
            raise _nc_core.NCBadInput('Invalid raw content (must be a string)')
        if not content:
            return
        self.__dirty()
        if 'raw_content' in self.__params:
            self.__params['raw_content'] += str(content)
        else:
            self.__params['raw_content'] = str(content)

    def get_raw_content( self ):
        return self.__params.get('raw_content','')

    def clear_raw_content( self ):
        if 'raw_content' in self.__params:
            self.__dirty()
            del self.__params['raw_content']

    def get_custom_section_data( self, section_name ):
        if section_name is None:
            dd = self.__params.get('custom_sections')
            return copy.deepcopy(dd) if dd else {}
        _check_valid_custom_section_name(section_name)
        if 'custom_sections' in self.__params:
            return self.__params['custom_sections'].get(section_name)

    def set_custom_section_data( self, section_name, content ):
        if section_name == 'HARDSPHERESANS':
            raise _nc_core.NCBadInput('For HARDSPHERESANS use the'
                                      ' .add_hard_sphere_sans_model() method'
                                      ' instead of set_custom_section_data.')
        if section_name == 'UNOFFICIALHACKS':
            raise _nc_core.NCBadInput('Do not set the @CUSTOM_UNOFFICIALHACKS'
                                      ' content directly with '
                                      'set_custom_section_data.')
        _check_valid_custom_section_name(section_name)
        if not isinstance(content,str):
            raise _nc_core.NCBadInput('Invalid custom section data'
                                      ' content (must be a string)')
        if self.__params.get('custom_sections',{}).get(section_name) == content:
            return
        self.__dirty()
        if 'custom_sections' not in self.__params:
            self.__params['custom_sections'] = {}
        self.__params['custom_sections'][section_name] = content

    def clear_custom_section_data( self, section_name ):
        if section_name is not None:
            _check_valid_custom_section_name(section_name)
        if 'custom_sections' not in self.__params:
            return
        if section_name is None:
            self.__dirty()
            del self.__params['custom_sections']
        else:
            if section_name not in self.__params['custom_sections']:
                return
            self.__dirty()
            del self.__params['custom_sections'][section_name]

    def __init__(self, *, data, fmt, quiet ):
        self._ichange = 0#increment on all changes, for downstream caching
        self.__loadcache = None
        self.__params = {}
        assert not (fmt is not None and data is None)
        if data is None:
            return
        if isinstance( data, NCMATComposerImpl ):
            self.__params = data.clone().__params
            return
        if fmt == 'via_ase':
            data= NCMATComposerImpl._ase_io_read( data, quiet = quiet )
            fmt = 'ase'
        if fmt is None and isinstance( data, _nc_core.TextData ) and data.dataType == 'ncmat':
            fmt = 'ncmat'
        if fmt is None:
            if isinstance(data,dict) or ( hasattr(data,'keys') and hasattr(data,'get') ):
                pass#Might be a dictionary from to_dict() but to be safe we do not autodetect this!
            elif isinstance(data,_nc_core.Info):
                fmt = 'info'
            elif hasattr(data,'__fspath__'):
                import pathlib
                _ = pathlib.Path(data).name.lower()
                if _.endswith('.ncmat'):
                    fmt='ncmat'
                if _.endswith('.cif'):
                    fmt='cif'
            elif hasattr(data,'startswith') and hasattr(data,'__contains__'):
                if '\n' in data:
                    if data.startswith('NCMAT'):
                        fmt='ncmat'
                    elif 'loop_' in data and '_atom_site' in data:
                        fmt='cif'
                elif data.startswith('codid::') or data.startswith('mpid::') or data.lower().endswith('.cif'):
                    fmt='cif'
                elif data.lower().endswith('.ncmat'):
                    fmt='cfgstr'
                else:
                    fmt='cfgstr'
            elif NCMATComposerImpl._is_ase_like_object( data ):
                fmt = 'ase'

        if fmt=='__internal_state__':
            self.__params = dict( copy.deepcopy( (k,v) ) for k,v in data.items() )
        elif fmt=='cif':
            o = NCMATComposerImpl.from_cif(data,quiet=quiet)
            self.__params = o.__params
        elif fmt=='ncmat':
            o = NCMATComposerImpl.from_ncmat(data)
            self.__params = o.__params
        elif fmt=='cfgstr':
            o = NCMATComposerImpl.from_cfgstr(data)
            self.__params = o.__params
        elif fmt=='info':
            o = NCMATComposerImpl.from_info(data)
            self.__params = o.__params
        elif fmt=='ase':
            o = NCMATComposerImpl.from_ase(data,quiet=quiet)
            self.__params = o.__params
        else:
            raise _nc_core.NCBadInput('Unknown data type (specify with the fmt keyword if this is just an auto-detection issue))')

    @staticmethod
    def from_cif( cifsrc, quiet, mp_apikey = None, uiso_temperature = None, override_spacegroup = None, **kwargs ):
        from . import cifutils as nc_cifutils
        lc = nc_cifutils.CIFLoader( cifsrc, quiet = quiet, mp_apikey = mp_apikey,
                                    override_spacegroup = override_spacegroup )
        return lc.create_ncmat_composer( quiet = quiet,
                                         uiso_temperature = uiso_temperature,
                                         **kwargs )._get_impl_obj()

    @staticmethod
    def from_cfgstr( cfgstr ):
        info = _nc_core.createInfo( cfgstr )
        return NCMATComposerImpl.from_info( info )

    @staticmethod
    def _is_ase_like_object( o ):
        return all( hasattr( o, e ) for e in ( 'cell', 'get_scaled_positions', 'get_atomic_numbers' ) )

    @staticmethod
    def _ase_io_read( obj, quiet, fmt = None ):
        if not quiet:
            _nc_common.print('Trying to load object via ase.io.read')
        ase, ase_io = _import_ase()
        return ase_io.read( obj, format = fmt )

    @staticmethod
    def from_ase( ase_obj, quiet, ase_format = None, **kwargs ):
        _nc_common.warn('Trying to load ASE object. Note that such objects might lack'
                        ' information including space group and atomic displacements.')
        if not NCMATComposerImpl._is_ase_like_object( ase_obj ):
            arg = _nc_common.extract_path( ase_obj )
            if arg is None:
                if ase_format is None:
                    _nc_common.warn('Anonymous text data might not be recognised by ASE.'
                                    ' Consider specifying the format explicitly in this case (e.g. supply ase_format="cif").')
                from .misc import AnyTextData
                td = AnyTextData(ase_obj)
                import io
                with io.BytesIO(td.content.encode()) as memfile:
                    ase_obj = NCMATComposerImpl._ase_io_read( memfile, quiet = quiet, fmt = ase_format )
            else:
                ase_obj = NCMATComposerImpl._ase_io_read( arg, quiet = quiet, fmt = ase_format )

            assert NCMATComposerImpl._is_ase_like_object( ase_obj )

        # Note1: We support this through CIF, since site occupancies are
        #        currently not available in a standard way through the OO interface
        #        (cf. https://gitlab.com/ase/ase/-/issues/263).
        #
        # Note2: Atomic displacement info seems to be ignored.
        #
        # Note3: D is detected as H (cf. https://matsci.org/t/how-to-load-structures-including-deuterium-with-ase/41950/2).
        cifdata = _cifdata_via_ase( ase_obj, ase_format = 'ase', quiet = quiet )
        if 'no_formula_check' not in kwargs:
            #ASE cif writer puts formula which does not take site occupancies into account:
            kwargs['no_formula_check'] = True
        o = NCMATComposerImpl.from_cif( cifsrc = cifdata, quiet=quiet, **kwargs )
        o.add_comments('NOTICE: Data here is based on an ASE (Atoms)'
                       ' structure through a conversion to CIF format.',
                       add_empty_line_divider=True)
        return o

    @staticmethod
    def from_info( info_obj ):
        return _composerimpl_from_info( info_obj )

    @staticmethod
    def from_ncmat( data, keep_header = False ):
        from .misc import AnyTextData
        td = AnyTextData(data)
        if '\n' not in td.content and not td.content.startswith('NCMAT'):
            return NCMATComposerImpl.from_cfgstr( td.content )
        o = _nc_core.directLoad( td, doScatter = False, doAbsorption = False )
        c = NCMATComposerImpl.from_info( o.info )
        comments = None
        if keep_header:
            comments = _extractInitialHeaderCommentsFromNCMATData(td)
            if comments:
                _hr='---------------------------------------------------------'
                comments = [
                    '>'*70,
                    'Comments found in NCMAT data from which the NCMATComposer',
                    'was initialised. Subsequent updates to the NCMATComposer',
                    'object might have invalidated some or all of these:'
                ]+ list('>>> '+e for e in comments) + ['<'*70]
        if comments:
            c.add_comments( comments, add_empty_line_divider = True )
        return c

    def update_atomdb( self, element_or_isotope, data, *,  mass=None,
                       coh_scat_len=None, incoh_xs=None, abs_xs=None ):
        label,res=_decode_update_atomdb( element_or_isotope,
                                         data=data,mass=mass,
                                         coh_scat_len=coh_scat_len,
                                         incoh_xs=incoh_xs,
                                         abs_xs=abs_xs )
        self.__dirty()
        if 'atomdb' not in self.__params:
            self.__params['atomdb'] = {}
        self.__params['atomdb'][label] = res

    def find_label( self, element, allow_multi ):
        import numbers
        from .atomdata import elementNameToZValue
        def _name2z( name ):
            return elementNameToZValue( name, allow_isotopes = True ) or None
        if isinstance(element,numbers.Integral):
            search_Z = int( element )
        else:
            search_Z = _name2z( element ) or 0
        if not 1<=search_Z<=150:
            return [] if allow_multi else None

        compos = self.__params.get('compositions',{})
        def lbl_zval_iter( lbl ):
            c = compos.get(lbl)
            for _,name in ( c or [(None,lbl)] ):
                z = _name2z(name)
                if z:
                    yield z
        ll = []
        for lbl in self.get_labels():
            if any( search_Z == z for z in lbl_zval_iter(lbl) ):
                ll.append( lbl )

        return ( list(sorted(ll))
                 if allow_multi
                 else ( ll[0] if len(ll)==1 else None ) )

    def get_state_of_matter( self ):
        return self.__params.get('state_of_matter')

    def set_state_of_matter( self, state_of_matter ):
        if not state_of_matter:
            state_of_matter = None
        allowed = ['solid','liquid','gas']
        if state_of_matter not in allowed and state_of_matter is not None:
            s='", "'.join(allowed)
            raise _nc_core.NCBadInput(f'Invalid state of matter value "{state_of_matter}" (allowed: "{s}")')
        if state_of_matter is None:
            if 'state_of_matter' not in self.__params:
                return
            self.__dirty()
            del self.__params['state_of_matter']
            return

        if self.__params.get('state_of_matter') == state_of_matter:
            return
        self.__dirty()
        self.__params['state_of_matter'] = state_of_matter

    def lock_temperature( self, value ):
        if value is None:
            if 'temperature' not in self.__params:
                return
            self.__params.pop('temperature',None)
            self.__dirty()
            return
        v=float(value)
        assert v>0.0
        self.__dirty()
        self.__params['temperature'] = { 'value' : v, 'lock' : True }

    def set_default_temperature( self, value ):
        if value is None:
            self.__params.pop('temperature',None)
            return
        v=float(value)
        assert v>0.0
        self.__dirty()
        self.__params['temperature'] = { 'value' : v, 'lock' : False }

    def get_temperature_setting( self ):
        _ = self.__params.get('temperature')
        return (_['value'],_['lock']) if _ else (None,None)

    def set_composition( self, label, *composition ):
        label, composition = _decode_composition(label,*composition)
        self.__dirty()
        if 'compositions' not in self.__params:
            self.__params['compositions'] = {}
        self.__params['compositions'][ label ] = composition

    def remap_atom( self,  element_or_isotope, *composition ):
        elem = _nc_common.check_elem_or_isotope_marker( element_or_isotope )
        if not elem:
            raise _nc_core.NCBadInput('Invalid element/isotope marker: "%s"'%element_or_isotope)

        elem, composition = _decode_composition(elem,*composition)
        compos = self.__params.get('compositions',{})
        set_compos_args = []
        for lbl in self.get_labels():
            lblcompos = compos.get(lbl,None)
            if lblcompos:
                frac_elem = sum( fr for fr,en in lblcompos if elem==en )
                if frac_elem > 0.0:
                    #update existing composition:
                    _ = [ (fr,en) for fr,en in lblcompos if elem!=en ]
                    _ += [ (frac_elem*fr,en) for fr,en in composition ]
                    set_compos_args.append( (lbl,_) )
            else:
                if lbl == elem:
                    set_compos_args.append( (lbl,composition) )
        for lbl, _compos in set_compos_args:
            self.set_composition( lbl, _compos )

    def clear_comments( self ):
        self.__dirty()#not strictly needed but to be safe
        self.__params['top_comments'] = []

    def add_comments( self, comments, add_empty_line_divider ):
        is_str = isinstance(comments,str) or hasattr(comments,'split')
        ll = comments if not is_str else list(comments.splitlines())

        self.__dirty()#not strictly needed as comments do not affect loading

        if 'top_comments' not in self.__params:
            self.__params['top_comments'] = []
        elif ( add_empty_line_divider
               and self.__params['top_comments']
               and self.__params['top_comments'][-1] ):
            self.__params['top_comments'].append('')

        for e in ll:
            self.__params['top_comments'].append( self.__prunecomment(e) )

    def __prunecomment(self,c):
        return str(c).rstrip().replace('\t','    ')

    def clone(self):
        o = NCMATComposerImpl( data = self.__params, fmt = '__internal_state__', quiet = True )
        o._ichange = self._ichange
        o.__loadcache = self.__loadcache
        return o

    def __dirty(self):
        self._ichange += 1
        self.__loadcache = None

    def get_cache_key( self ):
        return ( int(id(self), self._ichange) )

    def to_dict(self):
        return copy.deepcopy(self.__params)

    def set_cellsg( self, *, a,b,c, alpha,beta,gamma, spacegroup ):
        #TODO: input validation! Also ensure high precision in numbers below
        self.__dirty()
        self.__params['cellsg'] = dict(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,spacegroup=spacegroup)

    def set_cellsg_cubic( self, a, *, spacegroup ):
        if not (spacegroup is None or 195<=int(spacegroup)<=230):
            raise _nc_core.NCBadInput(f'Spacegroup {spacegroup} is not a cubic spacegroup (must be integer in range 195..230)')
        self.set_cellsg( a=a,b=a,c=a,alpha=90.,beta=90.,gamma=90.,
                         spacegroup = ( int(spacegroup) if spacegroup else None ) )

    def set_density( self, value, unit ):
        allowed = ['g/cm3','kg/m3','atoms/Aa3']
        if unit not in allowed:
            s='", "'.join(allowed)
            raise _nc_core.NCBadInput(f'Invalid density unit "{unit}" (allowed: "{s}")')
        if not ( value>0.0 ):
            raise _nc_core.NCBadInput(f'Invalid density value (must be >0): {value} {unit}')
        self.__dirty()
        self.__params['density'] = ( value, unit )

    def set_fraction( self, label, value ):
        _checklabel(label)
        if not 0 < value <= 1.0:
            raise _nc_core.NCBadInput(f'Invalid component fraction is not in (0,1]: {value}')
        self.__dirty()
        if 'fractions' not in self.__params:
            self.__params['fractions'] = {}
        self.__params['fractions'][label] = value

    def __lines_cellsg(self,cellsg):
        if not cellsg:
            return ''
        c = cellsg
        spacegroup = c.get('spacegroup',None)
        _sa = f"{c['a']:g}"
        _sb = f"{c['b']:g}"
        _sc = f"{c['c']:g}"
        if spacegroup and 195<=spacegroup<=230:
            if not ( _sa == _sb and _sb == _sc ):
                raise _nc_core.NCBadInput(f'Invalid lattice parameters for cubic spacegroup ({spacegroup}): a={_sa}, b={_sb}, c={_sc}')
            if not ( c['alpha']==90 and c['beta']==90 and c['gamma']==90 ):
                raise _nc_core.NCBadInput(f"Invalid lattice angles for cubic spacegroup ({spacegroup}): alpha={c['alpha']}, beta={c['beta']}, gamma={c['gamma']}")
            ll = f"""
            @CELL
            cubic {_sa}
            """
        else:
            if _sa == _sb:
                _sb = '!!'
                if _sc == _sa:
                    _sc = '!!'
            elif _sb == _sc:
                _sc = '!!'

            ll = f"""
            @CELL
            lengths {_sa} {_sb} {_sc}
            angles {c['alpha']:g} {c['beta']:g} {c['gamma']:g}
            """
        if spacegroup:
            ll+=f"""
            @SPACEGROUP
            {spacegroup}
            """
        return ll

    def __add_dyninfo( self, label, fraction, dyninfo ):
        _checklabel(label)
        self.__dirty()
        if 'dyninfos' not in self.__params:
            self.__params['dyninfos'] = { label : dyninfo }
        else:
            self.__params['dyninfos'][label] = dyninfo
        if fraction is not None:
            self.set_fraction( label, fraction )

    def allow_fallback_dyninfo( self, debye_temp ):
        dt = _nc_common._decodeflt( debye_temp )
        if not ( dt and dt>0.0 ):
            raise _nc_core.NCBadInput(f'Invalid Debye temperature value: {debye_temp}')
        self.__dirty()
        self.__params['fallback_debye_temp'] = dt

    def set_dyninfo_vdosdebye( self, label, debye_temp, *, comment, fraction ):
        self.__add_dyninfo( label, fraction, dict( ditype='vdosdebye',
                                                       debye_temp = float(debye_temp),
                                                       comment = None if not comment else self.__prunecomment(comment) ) )

    def set_dyninfo_freegas( self, label, *, comment, fraction ):
        self.__add_dyninfo( label, fraction, dict( ditype='freegas',
                                                       comment = None if not comment else self.__prunecomment(comment) ) )

    def set_dyninfo_sterile( self, label, *, comment, fraction ):
        self.__add_dyninfo( label, fraction, dict( ditype='sterile',
                                                       comment = None if not comment else self.__prunecomment(comment) ) )

    def set_dyninfo_vdos( self, label, vdos_egrid, vdos, *, comment, fraction ):
        self.__add_dyninfo( label, fraction, dict( ditype='vdos',
                                                       vdos_egrid = _copyarray_or_None(vdos_egrid),
                                                       vdos = _copyarray_or_None(vdos),
                                                       comment = None if not comment else self.__prunecomment(comment) ) )

    def set_dyninfo_scatknl( self, label, *, alphagrid, betagrid, temperature,
                             sab = None, sab_scaled = None, egrid = None,
                             comment = None, fraction = None ):
        def present( x ):
            return x is not None and _is_nonempty_array(x)
        if not present(sab) and not present(sab_scaled):
            raise _nc_core.NCBadInput('Missing either sab or sab_scaled arguments')
        if present(sab) and present(sab_scaled):
            raise _nc_core.NCBadInput('Do not specify both sab and sab_scaled')
        self.__add_dyninfo( label, fraction, dict( ditype='scatknl',
                                                       temperature = float(temperature),
                                                       egrid = _copyarray_or_None(egrid),
                                                       sab = _copyarray_or_None(sab),
                                                       sab_scaled = _copyarray_or_None(sab_scaled),
                                                       alphagrid = _copyarray_or_None(alphagrid),
                                                       betagrid = _copyarray_or_None(betagrid),
                                                       comment = None if not comment else self.__prunecomment(comment) ) )

    def set_dyninfo_msd( self, label, *, msd, temperature, comment, fraction ):
        self.__add_dyninfo( label, fraction, dict( ditype='msd',
                                                       msd_value = float(msd),
                                                       msd_temperature = float(temperature),
                                                       comment = None if not comment else self.__prunecomment(comment) ) )

    def set_dyninfo_from_object( self, label, source_dyninfo, comment = None, fraction = None ):
        di = source_dyninfo
        if not isinstance(di,_nc_core.Info.DynamicInfo):
            raise _nc_core.NCBadInput('source_dyninfo object must be an object derived from NCrystal.Info.DynamicInfo')
        _infoobj = di._info_wr() if di._info_wr else None
        if not _infoobj:
            raise _nc_core.NCBadInput('source_dyninfo object seems to be associated with expired Info object (or it was improperly initialised!)')

        if comment is None:
            comment = f'Transferred from "{di.atomData.displayLabel()}" in existing NCrystal.DynamicInfo object'


        if isinstance(di,_nc_core.Info.DI_VDOSDebye):
            self.set_dyninfo_vdosdebye( label,
                                        debye_temp=di.debyeTemperature(),
                                        comment=comment,
                                        fraction=fraction )
            return

        if isinstance(di,_nc_core.Info.DI_VDOS):
            self.set_dyninfo_vdos( label,
                                   vdos_egrid = di.vdosOrigEgrid(),
                                   vdos = di.vdosOrigDensity(),
                                   comment = comment,
                                   fraction = fraction)
            return

        if isinstance(di,_nc_core.Info.DI_Sterile):
            self.set_dyninfo_sterile( label,
                                      comment = comment,
                                      fraction = fraction )
            return

        if isinstance(di,_nc_core.Info.DI_FreeGas):
            self.set_dyninfo_freegas( label,
                                      comment = comment,
                                      fraction = fraction )
            return

        if isinstance(di,_nc_core.Info.DI_ScatKnlDirect):
            knl=di.loadKernel()
            #NB: Ignores knl['suggestedEmax'].
            #NB: From DI object we can only extract the sab, even if the input
            #file originally had e.g. sab_scaled.
            self.set_dyninfo_scatknl( label,
                                      temperature = di.temperature,
                                      alphagrid = knl['alpha'],
                                      betagrid = knl['beta'],
                                      sab = knl['sab'],
                                      egrid = knl['egrid'],
                                      comment = comment,
                                      fraction = fraction )
            return

        from .exceptions import nc_assert
        nc_assert(False,'set_dyninfo_from_object not implemented yet'
                  ' for %s type DynamicInfo'%di.__class__.__name__ )

    def transfer_dyninfo_objects( self, source, mapping, allow_none ):
        lbls = set(mapping.keys()) if mapping else self.get_labels()
        from ._ncmatimpl import extract_dyninfo_objects
        _keepalive, lbl2dyninfos = extract_dyninfo_objects( lbls,
                                                            self.__params.get('compositions',{}),
                                                            source,
                                                            mapping )
        if not lbl2dyninfos and not allow_none:
            raise _nc_core.NCBadInput('No dynamic info was transferred from source')
        for lbl, di in sorted( lbl2dyninfos.items() ):
            self.set_dyninfo_from_object( lbl, di['obj'], comment = di['comment'] )

    def as_spglib_cell( self ):
        cellsg = self.__params.get('cellsg',None)
        atomposlist = self.__params.get('atompos',{}).get('pos',[])
        if not cellsg or not atomposlist:
            return None, None
        lattice = _cellparams_to_spglib_lattice(cellsg)
        lbl2atomidx = {}
        atomic_points, atomic_types = [], []
        for lbl,x,y,z in sorted(atomposlist):
            if lbl not in lbl2atomidx:
                lbl2atomidx[lbl] = len(lbl2atomidx)
            #idx = lbl2atomidx[lbl]
            atomic_points.append( (x,y,z) )
            atomic_types.append( lbl2atomidx[lbl] )
        return ( lattice, atomic_points, atomic_types ), dict( (v,k) for k,v in lbl2atomidx.items())

    def refine_crystal_structure( self, symprec, quiet ):
        self._impl_refine( mode_refine = True, symprec = symprec, quiet=quiet )

    def verify_crystal_structure( self, symprec, quiet ):
        self._impl_refine( mode_refine = False, symprec = symprec, quiet=quiet )

    def _impl_refine( self, *, mode_refine, symprec, quiet ):
        cellsg = self.__params.get('cellsg',None)
        atompos = self.__params.get('atompos',None)
        #atomposlist = atompos.get('pos',None) if atompos else None
        if not cellsg and not atompos:
            return#do nothing if not crystal
        if int(bool(cellsg))+int(bool(atompos)) == 1:
            raise _nc_core.NCBadInput('Must set both unit cell parameters and atom'
                                      ' positions before trying to refine or verify a crystal structure')
        assert cellsg and atompos
        sgnumber = cellsg.get('spacegroup',None)
        if sgnumber is None and not mode_refine:
            raise _nc_core.NCBadInput('Must provide a space group number (or invoke'
                                      ' .refine_crystal_structure()) before it is'
                                      ' possible to verify a crystal structure')
        spglib_cell, spglib_idx2lbl = self.as_spglib_cell()
        d =  _spglib_refine_cell( spglib_cell, symprec = symprec ) if spglib_cell else None
        if not d:
            raise _nc_core.NCBadInput(f'Failed to {"refine" if mode_refine else "verify"} crystal structure with spglib.')
        assert len(d)==7
        if mode_refine:
            if not quiet:
                for m in d['msgs']:
                    _nc_common.print(m)
                for w in d['warnings']:
                    _nc_common.warn(w)
        if mode_refine:
            #NB: We do not (yet) have anisotropic atomic properties (like
            #anisotropic displacements), so we can simply update just the
            #positions. This is different from the case in cifutils.py where
            #d['can_keep_anisotropic_properties'] is important!
            new_cellsg = dict(d['cellparams_snapped'].items())
            new_cellsg['spacegroup'] = d['sgno']#NB: ignoring sgsymb_hm here!
                                                #Could we have spacegroup_hm as
                                                #optional field in the cellsg
                                                #array?
            ll = []
            for pos, atomidx in zip(d['refined_cell'][1],d['refined_cell'][2]):
                lbl = spglib_idx2lbl[atomidx]
                ll.append( (lbl,
                           _remap_fract_pos(pos[0]),
                           _remap_fract_pos(pos[1]),
                           _remap_fract_pos(pos[2]) ) )
            ll.sort()
            new_atompos = dict( (k,copy.deepcopy(v))
                                for k,v in atompos.items()
                                if k!='pos' )
            new_atompos['pos'] = ll
            self.__dirty()
            self.__params['cellsg'] = new_cellsg
            self.__params['atompos'] = new_atompos
        else:
            if sgnumber != d['sgno']:
                raise _nc_core.NCBadInput(f'Failed to verify crystal structure with spglib. Expected SG-{sgnumber}, got SG-{d["sgno"]} ({d["sgsymb_hm"]}).')
            rdl = _reldiff_cellparams( cellsg, d['cellparams_snapped'] )
            rda = _reldiff_atompos( spglib_cell, d['refined_cell'] )
            #rd = None if (rdl is None or rda is None) else max(rdl,rda)
            if rdl is None or rda is None or max(rdl,rda) > 0.01:
                raise _nc_core.NCBadInput('Failed to verify crystal structure with spglib.')

    def set_atompos( self, atompos ):
        pos,occumap = [],{}
        for e in atompos:
            assert len(e) in (4,5)
            lbl,x,y,z = str(e[0]),_remap_fract_pos(e[1]),_remap_fract_pos(e[2]),_remap_fract_pos(e[3])
            _checklabel(lbl)
            pos.append( (lbl,x,y,z) )
            occu = float(e[4]) if len(e)==5 else 1.0
            if not ( 0.0 < occu <= 1.0 ):
                raise _nc_core.NCBadInput('site_occupancy values must be in (0.0,1.0]')
            if lbl not in occumap:
                occumap[lbl]=occu
            else:
                if not occumap[lbl] == occu:
                    raise _nc_core.NCBadInput('All site occupancies for label "%s" are not identical'%lbl)

        #remove occu==1.0 entries from map (and sort):
        occumap = dict( sorted( (k,v) for k,v in occumap.items() if v!=1.0 ) )
        if occumap:
            _nc_common.warn('Support for site_occupancy is highly experimental. Do *not* attempt to directly override'
                            +' value with the "density" cfg-parameter for this material (unless it is to simply scale it)')
        self.__dirty()
        self.__params['atompos'] = dict( pos=pos, occumap=occumap )
#
    def get_labels( self ):
        #from positions:
        atompos = self.__params.get('atompos',None)
        s = set( lbl for lbl,x,y,z in atompos['pos'] ) if atompos else set()
        #from dyninfo:
        for lbl in self.__params.get('dyninfos',{}).keys():
            s.add(lbl)
        #from composition:
        for lbl in self.__params.get('compositions',{}).keys():
            s.add(lbl)
        #from fractions:
        for lbl in self.__params.get('fractions',{}).keys():
            s.add(lbl)
        return s

    def __lines_atompos(self, lbl_map, atompos ):
        if not atompos:
            return ''
        ll=[]
        fmt = _nc_common.prettyFmtValue
        #NB: We sort on the val->string->val values, not the original
        #values. Because otherwise we get irreproducibilities when in
        #principle two x values are almost identical and their order might
        #vary slightly.
        def extractval(s):
            if '/' not in s:
                return float(s)
            p=s.split('/')
            assert len(p)==2
            return float(p[0]) / float(p[1])
        for lbl,x,y,z in atompos['pos']:
            fx,fy,fz = fmt(x), fmt(y), fmt(z)
            ll.append( ( lbl_map.get(lbl,lbl),
                         extractval(fx), extractval(fy), extractval(fz),
                         fx, fy, fz ) )
        ll.sort()
        out="@ATOMPOSITIONS\n"
        for lbl,_,_,_,fx,fy,fz in ll:
            out += f'{lbl} {fx} {fy} {fz}\n'
        return out

    def __lines_dyninfo(self,lbl_map,fractions):
        dyninfos = self.__params.get('dyninfos',{})

        lines=''
        natoms_with_fallback_dyninfo = 0

        def formatVector(name,values):
            return formatVectorForNCMAT(name,values,_magic_two_space)

        def transform_msd_to_vdosdebye( lbl, dyninfo ):
            if not dyninfo or not dyninfo.get('ditype','') == 'msd':
                return dyninfo
            d = {}
            d['ditype'] = 'vdosdebye'
            mass = calc_mass(lbl,self.__params.get('compositions',{}).get(lbl,None))
            if not mass:
                raise _nc_core.NCBadInput( f'Failed to calculate the atomic mass associated with the label {lbl}'
                                           ' (needed for dyninfo of type "msd"). Please ensure that the label has a'
                                           ' proper composition set (with .set_composition(..)) and that all component'
                                           ' atoms have known masses (otherwise use .update_atomdb(..) to provide them)')
                assert mass>0.0
            from .vdos import debyeTempFromIsotropicMSD
            msd = dyninfo['msd_value']
            msd_temp = dyninfo['msd_temperature']
            dt = debyeTempFromIsotropicMSD( msd = msd, mass = mass,
                                            temperature = msd_temp )
            assert dt > 0.0
            d['debye_temp'] = dt
            c = f'Debye temperature value derived from msd={msd:g}Aa^2 @ T={msd_temp:g}K (mass={mass:g}u)'
            oldcomment = dyninfo.get('comment',None)
            d['comment'] = c if not oldcomment else f'{oldcomment.strip()}\n{c}'
            return d

        for lbl,fracinfo in sorted(fractions.items()):
            if not isinstance(fracinfo,float) and len(fracinfo)==3:
                _fracval_flt,frac_prettyprint,_ = fracinfo
            else:
                _fracval_flt,frac_prettyprint = fracinfo, _nc_common.prettyFmtValue( fracinfo )
            lines += '@DYNINFO\n'
            dyninfo = dyninfos.get(lbl,None)
            dyninfo = transform_msd_to_vdosdebye( lbl, dyninfo )

            if not dyninfo:
                fbdt = self.__params.get('fallback_debye_temp',None)
                if fbdt is not None:
                    natoms_with_fallback_dyninfo += 1
                    dyninfo = dict( ditype='vdosdebye',
                                    debye_temp = fbdt,
                                    comment = 'WARNING: Using fallback Debye temperature value!' )
            if not dyninfo:
                raise _nc_core.NCMissingInfo(f'Missing dyninfo for atom label "{lbl}" (add via call to set_dyninfo_...)')
            comment = dyninfo['comment']
            if comment:
                for e in comment.splitlines():
                    lines += f'  # {e.strip()}\n'

            lines += f'element {lbl_map.get(lbl,lbl)}\n'
            lines += f'fraction {frac_prettyprint}\n'
            ditype = dyninfo['ditype']
            lines += f'type {ditype}\n'
            if ditype == 'vdosdebye':
                dt=dyninfo['debye_temp']
                lines += f'debye_temp {dt:g}\n'
            elif ditype == 'vdos':
                x = dyninfo['vdos_egrid']
                from ._common import _grid_is_linspace
                if len(x)>2 and _grid_is_linspace(x):
                    lines += f'vdos_egrid {_fmtprecisenum(x[0])} {_fmtprecisenum(x[-1])}\n'
                else:
                    lines += formatVector('vdos_egrid',x)
                lines += formatVector('vdos_density',dyninfo['vdos'])
            elif ditype in ('sterile','freegas'):
                pass
            elif ditype == 'scatknl':
                lines += f'temperature {_fmtprecisenum(dyninfo["temperature"])}\n'
                if dyninfo.get('egrid') is not None:
                    from ._common import _grid_is_linspace
                    x = dyninfo['egrid']
                    if ( len(x) > 3 and _grid_is_linspace(x) ):
                        lines += f'egrid {_fmtprecisenum(x[0])} {_fmtprecisenum(x[-1])} {len(x)}\n'
                    elif len(x) <= 3:
                        if len(x)==3 and x[0]==0.0 and x[2]==0.0:
                            x = [ x[1] ]#single Emax value
                        _=' '.join(_fmtprecisenum(e) for e in x)
                        lines += f'egrid {_}\n'
                    else:
                        lines += formatVector('egrid',x)
                lines += formatVector('alphagrid',dyninfo['alphagrid'])
                lines += formatVector('betagrid',dyninfo['betagrid'])
                if dyninfo['sab'] is not None:
                    assert dyninfo['sab_scaled'] is None
                    lines += formatVector('sab',dyninfo['sab'])
                else:
                    assert dyninfo['sab_scaled'] is not None
                    lines += formatVector('sab_scaled',dyninfo['sab_scaled'])
            else:
                assert False,f'Unexpected dyninfo type : "{ditype}"'
        return lines, natoms_with_fallback_dyninfo

    def register_as( self, virtual_filename, cfg_params ):
        from .datasrc import registerInMemoryFileData
        registerInMemoryFileData( virtual_filename, self.create_ncmat( cfg_params = cfg_params ) )
        return virtual_filename

    def write( self, path, cfg_params ):
        import pathlib
        p = pathlib.Path( path )
        assert p.name.endswith('.ncmat') and len(p.name) > 6
        _nc_common.write_text( p,
                               self.create_ncmat( cfg_params = cfg_params ) )
        return p

    def load(self,cfg_params,force ):
        if force or self.__loadcache is None or self.__loadcache[0] != cfg_params:
            self.__loadcache = ( cfg_params, _nc_core.directLoad( self.create_ncmat(), cfg_params ) )
        return self.__loadcache[1]

    def plot_xsect( self, composer, cfg_params, kwargs_plot_xsect ):
        from .misc import MaterialSource
        from .plot import plot_xsect
        ms = MaterialSource( composer, cfg_params = cfg_params )
        return plot_xsect( ms, **kwargs_plot_xsect )

    def inspect( self, composer, cfg_params, kwargs_plot_xsect ):
        from .misc import MaterialSource
        from .plot import plot_xsect
        ms = MaterialSource( composer, cfg_params = cfg_params )
        loadedmat = ms.load()
        import sys
        def flush():
            ( sys.stdout.flush(),sys.stderr.flush() )
        flush()
        loadedmat.info.dump()
        flush()
        _nc_common.print("Absorption process (objects):")
        flush()
        loadedmat.absorption.dump(prefix='  ')
        flush()
        _nc_common.print("Scattering process (objects):")
        flush()
        loadedmat.scatter.dump(prefix='  ')
        flush()
        return plot_xsect( ms, **kwargs_plot_xsect )

    def __determine_atompos_fractions( self, atompos ):
        natoms = len(atompos['pos']) if atompos else 0
        if not natoms:
            return
        out = {}
        from collections import Counter as _collections_Counter
        counts = _collections_Counter(list(e[0] for e in atompos['pos']))
        for lbl,count in counts.items():
            frac = float(count) / natoms
            gcd = math.gcd(count,natoms)
            out[ lbl ] = ( frac, ( '1' if count==natoms else '%i/%i'%( count // gcd, natoms // gcd ) ), count )
        return out

    def get_chemical_composition( self, as_str ):
        atompos = self.__params.get('atompos')
        _lbl_counts = None
        if atompos:
            _afr = self.__determine_atompos_fractions( atompos )
            _lbl_counts = list( (lbl,count) for lbl,(_,_,count) in _afr.items() )
        elif self.__params.get('cellsg'):
            #We could in principle proceed below, but for crystals, we promise that the
            #return value is always numbers per unit cell!
            raise _nc_core.NCMissingInfo("Can not provide chemical composition of incomplete crystalline materials (missing atomic positions).")
        if _lbl_counts is None:
            fractions = self.__params.get('fractions')
            dyninfos = self.__params.get('dyninfos')
            if not fractions and not atompos and len(dyninfos)==1:
                fractions = {list(dyninfos.keys())[0] : 1.0 }
            if fractions:
                _lbl_counts = list( sorted( fractions.items()) )
        if not _lbl_counts:
            return None

        compos = self.__params.get('compositions')
        occumap = atompos['occumap'] if atompos else None
        ll = []
        for lbl, weight in _lbl_counts:
            c = compos.get(lbl) if compos else None
            _occufactor = occumap.get(lbl,1.0) if occumap else 1.0
            if not c:
                assert _nc_common.check_elem_or_isotope_marker( lbl ) is not None
                ll.append( (lbl, weight * _occufactor ) )
            else:
                for frac,elemisomarker in c:
                    ll.append( (elemisomarker, weight * frac * _occufactor ) )
        d = {}
        for k,v in ll:
            if k in d:
                d[k] += v
            else:
                d[k] = v
        res = list( sorted( (k,v) for k,v in d.items() ) )
        if as_str:
            return _nc_common.format_chemform( res )
        else:
            return res

    def _unofficial_vdos2sab_ignore( self, *, order_low, order_high = None, mode = None ):
        assert mode is None or mode in ('coherent','incoherent')
        import numbers
        assert isinstance(order_low,numbers.Integral)
        assert order_high is None or isinstance(order_low,numbers.Integral)
        order_high = int(order_high) if order_high is not None else None
        order_low = int(order_low)
        assert order_high is None or order_high >= order_low
        line = ['vdos2sab_ignorecontrib',str(order_low)]
        if order_high is not None:
            line.append( str(order_high) )
        if mode is not None:
            line.append(mode)
        old = [e for e in self.__params.get('unofficial_hacks',[]) if ( e and not e[0]=='vdos2sab_ignorecontrib' ) ]
        self.__dirty()
        self.__params['unofficial_hacks'] = old  + [ line ]

    def create_ncmat( self, *, cfg_params='', meta_data = False,
                      verify_crystal_structure = True ):
        #docstr: set meta_data to also return a dictionary with a bit of high
        #level metadata like chemical composition, spacegroup, etc.
        md = {} if meta_data else None

        cellsg = self.__params.get('cellsg')
        atompos = self.__params.get('atompos')

        is_crystal = bool(cellsg or atompos)
        did_verify_xtal_struct = False
        if verify_crystal_structure and cellsg and atompos:
            self.verify_crystal_structure(symprec = 0.01, quiet = False)
            did_verify_xtal_struct = True

        if verify_crystal_structure and cellsg and cellsg.get('spacegroup',None):
            _check_cell_sg_consistency( cellsg['spacegroup'],
                                        cellsg['a'], cellsg['b'], cellsg['c'],
                                        cellsg['alpha'], cellsg['beta'], cellsg['gamma'] )

        dyninfos = self.__params.get('dyninfos',None)

        if md is not None:
            md['cellsg'] = dict( cellsg.items() ) if cellsg else None

        if atompos and not cellsg:
            raise _nc_core.NCMissingInfo("Atom positions in unit cell added but unit cell info missing (add via call to set_cellsg)")
        if cellsg and not atompos:
            raise _nc_core.NCMissingInfo("Cell/spacegroup information added but atom positions in unit cell are missing (add via call to set_atompos)")
        is_crystal = bool(cellsg or atompos)
        if md is not None:
            md['is_crystal'] = is_crystal

        density = self.__params.get('density',None)
        fractions = self.__params.get('fractions',None)
        if is_crystal and density:
            raise _nc_core.NCBadInput('Density should not be added explicitly for crystalline materials (if'
                                +' needed, it can be modified via the "density" cfg-string parameter)')

        if is_crystal and dyninfos:
            for lbl,di in sorted(dyninfos.items()):
                _ditype=di.get('ditype',None)
                if _ditype not in ('vdos','vdosdebye','msd'):
                    raise _nc_core.NCBadInput('Crystalline material are currently only supported with dynamic info type of'
                                        +f' either "vdos", "vdosdebye", or "msd" (offending type "{_ditype}" for "{lbl}")')

        if not is_crystal and len(dyninfos)==1 and not fractions:
            #in case of a single dyninfo we can assume the fraction is 1.0
            fractions = {list(dyninfos.keys())[0] : 1.0 }


        atompos_fractions = None if not is_crystal else self.__determine_atompos_fractions( atompos )
        if atompos_fractions and fractions:
            for lbl,frac in sorted(fractions.items()):
                afrac,_,_ = atompos_fractions.get(lbl,None)
                if afrac is None:
                    raise _nc_core.NCBadInput(f'Label "{lbl}" has fraction explicitly set, but the material'
                                        +' is crystalline and the label has no atomic positions sets')
                if abs(afrac-frac)>1e-5:
                    raise _nc_core.NCBadInput(f'Label "{lbl}" has both explicit and implicit (via atomic'
                                              +f' position count) fractions provided, and they do NOT match up ({afrac} vs. {frac})')
        if not is_crystal:
            if not density:
                raise _nc_core.NCBadInput('Density must be set explicitly for non-crystalline materials (add via call to set_density)')
            if not dyninfos:
                raise _nc_core.NCBadInput('Material incompletely specified.')
            if set(dyninfos.keys()) != set(lbl for lbl,fv in (fractions or {}).items() if fv is not None):
                raise _nc_core.NCBadInput('For non-crystalline materials with more than one component, all components must have fractions specified.')
            _fracsum = math.fsum(val for lbl,val in sorted((fractions or {}).items()))
            if abs( 1.0 - _fracsum ) > 1e-10:
                raise _nc_core.NCBadInput(f'Invalid material - sum of dyninfo fractions is not 1 (it is {_fracsum})')

        lbl_map,atomdb_lines = determine_labels_and_atomdb( self.__params, fractions = fractions )

        ll  = self.__lines_cellsg( cellsg )
        ll += self.__lines_atompos( lbl_map, atompos )
        if density:
            _dval,_dunit = density
            _ncmat_dunit = {'g/cm3': 'g_per_cm3',
                            'kg/m3': 'kg_per_m3',
                            'atoms/Aa3': 'atoms_per_aa3' }[_dunit]
            ll+=f'@DENSITY\n{_dval} {_ncmat_dunit}\n'

        som = self.get_state_of_matter()
        if som is not None:
            assert som in ('solid','liquid','gas')
            ll += f'@STATEOFMATTER\n{som}\n'

        _t = self.__params.get('temperature',None)
        _v = '%.14g'%_t['value'] if _t else None
        if _t and ( _t['lock'] or _v != '293.15' ):
            ll += '@TEMPERATURE\n'
            if not _t['lock']:
                ll += 'default '
            ll += f'{_v}\n'


        if atomdb_lines:
            ll += '@ATOMDB\n'
            ll += '\n'.join(atomdb_lines)
            ll += '\n'

        secondary_phases = self.__params.get('secondary_phases')
        if secondary_phases:
            ll += '@OTHERPHASES\n'
            for frac,cfgstr in secondary_phases:
                ll += f'{frac} {cfgstr}\n'

        custom_hardspheresans = self.__params.get('custom_hardspheresans')
        if custom_hardspheresans:
            if not secondary_phases:
                raise _nc_core.NCBadInput('Material with hard-sphere SANS'
                                          ' enabled must have at least one'
                                          ' secondary phase added.')
            ll += '@CUSTOM_HARDSPHERESANS\n'
            ll += f'{custom_hardspheresans} #sphere radius in angstrom.\n'

        unofficial_hacks = self.__params.get('unofficial_hacks')
        if unofficial_hacks:
            ll += '@CUSTOM_UNOFFICIALHACKS\n'
            for e in unofficial_hacks:
                ll += ' '.join(e)+'\n'

        for sn,cnt in sorted(self.__params.get('custom_sections',{}).items()):
            ll += f'@CUSTOM_{sn}\n'
            ll += cnt
            if not cnt.endswith('\n'):
                ll += '\n'

        #Dyninfo lines last, since they might contain huge arrays of data, and
        #people might only look at the top of the file:
        ld, natoms_with_fallback_dyninfo = (
            self.__lines_dyninfo( lbl_map, atompos_fractions or fractions )
        )
        ll += ld

        out=["NCMAT v7"]
        comments = copy.deepcopy(self.__params.get('top_comments',[]))

        if natoms_with_fallback_dyninfo:
            _ = 'values were' if natoms_with_fallback_dyninfo>1 else 'value was'
            if comments and comments[-1]:
                comments.append('')
            comments += ['WARNING: Fallback (dummy) Debye temperature %s used for '%_
                         +'%i atom%s!'%(
                             natoms_with_fallback_dyninfo,
                             's' if natoms_with_fallback_dyninfo>1 else ''),'']

        if did_verify_xtal_struct:
            comments += ['NOTICE: crystal structure was verified with spglib to be self-consistent.']
        elif is_crystal:
            comments += ['WARNING: crystal structure was not verified automatically with spglib.']


        #determine chemical formula
        if atompos_fractions is not None:
            _lbl_counts = list( (lbl,count) for lbl,(_,_,count) in atompos_fractions.items() )
        else:
            assert fractions
            _lbl_counts = list( sorted( fractions.items()) )

        _chem_compos = self.get_chemical_composition( as_str = False)
        _chemform = _nc_common.format_chemform( _chem_compos )

        _title =  _chemform
        if cellsg:
            _sgno = cellsg['spacegroup']
            if _sgno:
                _title += f' ({_nc_common._classifySG(_sgno)}, SG-{_sgno})'

        if md is not None:
            md['chemform'] = _chemform
            #md['chemform_per_unit_cell'] = _chemform_per_cell
            md['title'] = _title

        disable_autotitle = False
        if comments and comments[0]=='<<disableautotitle>>':
            #hidden feature (might not actually be used yet)
            comments = comments[1:]
            disable_autotitle = True

        if comments and comments[0]=='<<disableautogennotice>>':
            #hidden feature for our command-line scripts.
            comments = comments[1:]
        else:
            _ = self.__class__.__name__
            if _.endswith('Impl') and len(_)>4:
                _ = _[:-4]
            out.append('# Autogenerated by %s'%_)

        if not disable_autotitle:
            out += [ '#',f'# {_title}','#' ]
            if is_crystal:
                _ = '+'.join(f'{count:g}x{lbl}' for lbl,count in _chem_compos)
                out += [ f'# Atoms per unit cell: {_}','#' ]

        for tc in comments:
            out.append(('# %s'%tc).rstrip())

        if cfg_params:
            if out and out[-1] != '#':
                out += [ '#' ]
            out += [ f'# NCRYSTALMATCFG[{str(cfg_params).strip()}]', '#']

        if out and out[-1] != '#':
            out += [ '#' ]

        for e in ll.splitlines():
            e=e.strip()
            if e.startswith('@'):
                out.append(e)
            elif e:
                out.append('  '+e.replace(_magic_two_space,'  '))

        out = '\n'.join(e.rstrip() for e in out)+'\n'

        raw_cnt = self.__params.get('raw_content')
        if raw_cnt:
            out += raw_cnt
        return out if md is None else ( out, md )

def _determine_dyninfo_mapping( labels, composition, dilist ):
    #mapping is composer label -> srclbl

    def extract_z2lbl( datalist, extractfct, allow_multi = False  ):
        z2lbl = {}
        for d in datalist:
            atomdata_or_z,lbl = extractfct(d)
            if atomdata_or_z is None:
                z = None
            elif hasattr(atomdata_or_z,'isElement'):
                z = atomdata_or_z.Z() if atomdata_or_z.isElement() else None
            else:
                z = int(atomdata_or_z) if atomdata_or_z else None
            if z is not None:
                if z in z2lbl:
                    z2lbl[z].append(lbl)
                else:
                    z2lbl[z] = [ lbl ]
        if allow_multi:
            return z2lbl
        return dict( (k,v[0]) for k,v in z2lbl.items() if len(v)==1 )

    def lookup_zval_for_label(lbl):
        c = composition.get(lbl,None)
        if not c:
            return _atomdb(lbl)
        zvals = set()
        for frac,s in c:
            ad = _atomdb(s)
            zvals.add(ad.Z() if (ad and ad.isElement()) else None)
        if len(zvals)==1:
            return zvals.pop()

    z_2_dilbl = extract_z2lbl( dilist, lambda di : (di.atomData, di.atomData.displayLabel()) )
    z_2_atomlbls = extract_z2lbl( labels, lambda lbl : ( lookup_zval_for_label(lbl), lbl ), allow_multi = True )
    d = {}
    for z,lbls in z_2_atomlbls.items():
        dilbl = z_2_dilbl.get(z,None)
        if dilbl:
            for lbl in lbls:
                d[lbl] = dilbl
    return d

def extract_dyninfo_objects( labels, compositions, source, mapping ):

    keepalive = None
    if isinstance( source, _nc_core.Info.DynamicInfo ):
        description, dilist = 'DynamicInfo object', [source]
    elif ( hasattr(source,'__len__')
         and not hasattr(source,'keys') #NB: '__len__' but not 'keys': list-like but not dict-like.
         and ( len(source)==0 or all(isinstance(e, _nc_core.Info.DynamicInfo) for e in source ) ) ):
        #sequence of dyninfo objects
        description, dilist = 'list of DynamicInfo objects', list( e for e in source )
    else:
        from .misc import MaterialSource
        ms = MaterialSource(source)
        description = ms.description
        _loadedmat = ms.load( doScatter = False, doAbsorption = False )
        keepalive = _loadedmat
        if not _loadedmat.info:
            raise _nc_core.NCBadInput('Material source has no Info object'
                                      ' (and therefore no DynamicInfo objects)')
        dilist = _loadedmat.info.dyninfos

    if not mapping:
        #mapping is srclbl -> label
        mapping = _determine_dyninfo_mapping( labels, compositions, dilist )

    out = {}
    dilbl_2_di = dict( (di.atomData.displayLabel(),di) for di in dilist )
    for lbl in labels:
        target_dilbl = mapping.get(lbl,None)
        if not target_dilbl:
            continue
        di = dilbl_2_di.get(target_dilbl,None)
        if not di:
            _='", "'.join(sorted(dilbl_2_di.keys()))
            raise _nc_core.NCBadInput(f'No display label "{target_dilbl}" found in source (source has display labels "{_}").')
        comment = f'Transferred from "{target_dilbl}" in "{description}"'
        out[lbl] = dict( obj=di, comment=comment )

    return keepalive, out

def _atomdb(lbl):
    from .atomdata import atomDB
    return atomDB( lbl, throwOnErrors=False )

def calc_mass( label, composition ):
    if not composition:
        _ = _atomdb( label )
        return _.averageMassAMU() if _ else None
    elif len( composition )==1:
        from .atomdata import atomDB
        return atomDB( composition[0][1] ).averageMassAMU()
    else:
        ll=[]
        fracs = []
        for frac, _lbl in composition:
            _ = _atomdb( _lbl )
            if not _:
                return None
            fracs.append( frac )
            ll.append( frac * _.averageMassAMU() )
        return math.fsum( ll ) / math.fsum( fracs) #NB: sum(fracs) is likely
                                                   #already guaranteed 1 here,
                                                   #playing it safe

def _composerimpl_from_info( infoobj ):
    if isinstance(infoobj,_nc_core.Info):
        i = infoobj
    elif hasattr(infoobj,'info') and isinstance(infoobj.info,_nc_core.Info):
        i = infoobj.info
    else:
        raise _nc_core.NCBadInput('NCMATComposer.from_info got unexpected argument of wrong type')

    if not i.isSinglePhase():
        raise _nc_core.NCBadInput('NCMATComposer can only be initialsed from single-phase materials')
    o = NCMATComposerImpl(data=None,fmt=None,quiet=True)
    top_atomdata = []

    #Transfer temperature. In case of gasses we also lock the value, since
    #density depends very strongly on temperature for such materials (it is of
    #course true for any state of matter, but for gasses it is extreme):
    if i.stateOfMatter() == _nc_core.StateOfMatter.Gas:
        o.lock_temperature( i.getTemperature() )
    else:
        o.set_default_temperature( i.getTemperature() )

    if i.isCrystalline():
        si = i.structure_info
        si_natoms = si['n_atoms']
        o.set_cellsg( a = si['a'], b = si['b'], c = si['c'],
                      alpha = si['alpha'], beta = si['beta'], gamma = si['gamma'],
                      spacegroup = ( si['spacegroup'] or None ) )
        all_atompos = []
        for ai in i.atominfos:
            ad = ai.atomData
            top_atomdata.append(ad)
            lbl = ad.displayLabel()
            di = ai.dyninfo
            o.set_dyninfo_from_object( lbl, di )
            for x,y,z in ai.positions:
                all_atompos.append( (lbl,x,y,z) )#NB: site_occupancy=1 here!
        from .exceptions import nc_assert
        nc_assert(len(all_atompos)==si_natoms)
        o.set_atompos( all_atompos )
    else:
        o.set_density( i.density, unit='g/cm3' )
        if i.stateOfMatter() == _nc_core.StateOfMatter.Solid:
            o.set_state_of_matter( 'solid' )
        elif i.stateOfMatter() == _nc_core.StateOfMatter.Liquid:
            o.set_state_of_matter( 'liquid' )
        elif i.stateOfMatter() == _nc_core.StateOfMatter.Gas:
            o.set_state_of_matter( 'gas' )
        else:
            assert i.stateOfMatter() == _nc_core.StateOfMatter.Unknown
        for di in i.dyninfos:
            ad = di.atomData
            top_atomdata.append(ad)
            lbl = ad.displayLabel()
            o.set_dyninfo_from_object( lbl, di )
            o.set_fraction( lbl, di.fraction )

    #Handle ATOMDB:
    def _extractad( ad ):
        if not ad.isComposite():
            return [ ( 1.0, ad )], [ ad ]
        basicads = []
        compos = []
        for comp in ad.components:
            _frac, _ad = comp.fraction, comp.data
            if not _ad.isComposite():
                compos.append( ( _frac, _ad ) )
                basicads.append( _ad )
            else:
                _subcompos, _subbasicads = _extractad( _ad )
                basicads += _subbasicads
                for _subfrac, _subad in _subcompos:
                    compos.append( (_frac*_subfrac, _subad) )
        return compos, basicads

    basic_atomdata = []#all non-composite atomdatas
    for ad in top_atomdata:
        _compos, _bads  = _extractad( ad )
        basic_atomdata +=  _bads
        o.set_composition( ad.displayLabel(), [(frac,ad.description(False)) for frac,ad in _compos] )

    #any non-builtin atom data:
    _seen = {}
    for ad in basic_atomdata:
        assert ad.isElement() or ad.isIsotope()
        key = (ad.Z(),ad.A())
        _adstr = ad.to_atomdb_str()
        _seenstr = _seen.get(key,None)
        if _seenstr is not None:
            if _adstr != _seenstr:
                raise _nc_core.NCBadInput('Atom with (Z,A)=(%i,%i) appears in multiple roles with different data values in material. Such materials are not supported by the NCMATComposer.'%key)
            continue
        _seen[ key ] = _adstr
        from .atomdata import atomDB
        _builtin = atomDB(Z=key[0],A=key[1],throwOnErrors=False)
        if not _builtin or _adstr != _builtin.to_atomdb_str():
            o.update_atomdb( ad.description(False), _adstr )

    #Custom sections:
    for nm,cnt_lists in (i.customsections or []):
        if nm == 'HARDSPHERESANS':
            _nc_common.warn('Ignoring @CUSTOM_HARDSPHERESANS sections in input.')
            continue
        if nm == 'UNOFFICIALHACKS':
            _nc_common.warn('Ignoring @CUSTOM_UNOFFICIALHACKS sections in input.')
            continue
        if o.get_custom_section_data(nm) is not None:
            raise _nc_core.NCBadInput(f'Multiple @CUSTOM_{nm} sections in input'
                                      'is not supported by NCMATComposer.')
        cnt=''
        for linedata in cnt_lists:
            cnt += ' '.join(linedata)
            cnt += '\n'
        o.set_custom_section_data(nm,cnt)

    return o

def _checklabel(lbl):
    if not lbl:
        raise ValueError(f'invalid label "{lbl}"')

def _decode_composition(label,*composition):
    _checklabel(label)

    label = ' '.join(str(label).split())
    errmsg='Invalid composition syntax'

    if ' is ' in label:
        _=label.split()
        if len(_)>2 and _[1]=='is':
            if composition:
                raise _nc_core.NCBadInput(errmsg)
            return _decode_composition( _[0], _[2:] )

    if not composition:
        raise _nc_core.NCBadInput(errmsg)

    _decodeflt = _nc_common._decodeflt
    #We want composition to end up in the final form: [(f0,name0),(f1,name1),...].
    def is_final_form(c):
        return ( all(hasattr(e,'__len__') for e in c)
                 and all(len(e)==2 for e in c)
                 and all(_decodeflt(e[0]) for e in c)
                 and all(hasattr(e[1],'split') for e in c) )

    single_arg = composition[0] if len(composition)==1 else None
    if single_arg and is_final_form(single_arg):
        single_arg,composition = None,single_arg

    if single_arg:
        if hasattr(single_arg,'split'):
            #assume string form like: 'Al' or '0.2 Al 0.8 Cr'
            p = single_arg.split()
            if len(p)!=1:
                return _decode_composition( label, *p )
            #single element form like: 'Al'
            norm_ident = _nc_common.check_elem_or_isotope_marker( p[0] )
            if not norm_ident:
                raise _nc_core.NCBadInput(errmsg+': invalid element/isotope marker "%s"'%p[0])
            return label, [(1.0,norm_ident)]
        elif ( hasattr(single_arg,'__len__')
               and len(single_arg)>=2 and len(single_arg)%2==0
               and all(_decodeflt(e) for e in single_arg[::2])
               and all(hasattr(e,'split') for e in single_arg[1::2]) ):
            #assume single_arg was a a list form like: [0.2,'Al',0.8,Cr]
            return _decode_composition( label, list(zip(single_arg[::2],single_arg[1::2])) )
        #unknown single arg form
        raise _nc_core.NCBadInput(errmsg)


    assert len(composition)>1 or (len(composition)==1 and is_final_form(composition))

    if ( not is_final_form(composition)
         and len(composition)%2 == 0
         and all(hasattr(e,'split') for e in composition[1::2])
         and all(_decodeflt(e) for e in composition[::2]) ):
        #assume form [f0,name0,f1,name1,...]
        composition = list(zip(composition[::2],composition[1::2]))

    if not is_final_form(composition):
        raise _nc_core.NCBadInput(errmsg)

    #composition is now in the form: [(f0,name0),(f1,name1),...].

    ll=[]

    for frac_orig, ident in composition:
        frac_val = _decodeflt( frac_orig )
        if frac_val is None:
            raise _nc_core.NCBadInput(errmsg+': invalid fraction specification "%s"'%frac_orig)
        if not (0<frac_val<=1.0):
            raise _nc_core.NCBadInput(errmsg+': fraction specification "%s" is not in (0,1]'%frac_orig)
        norm_ident = _nc_common.check_elem_or_isotope_marker( ident )
        if not norm_ident:
            raise _nc_core.NCBadInput(errmsg+': invalid element or isotope identifier "%s"'%ident)
        ll.append( ( frac_val, norm_ident ) )
    fractot = math.fsum(f for f,lbl in ll)
    if abs(fractot-1.0)>1e-5:
        raise _nc_core.NCBadInput(errmsg+': fractions do not sum to 1')

    #Now merge identical entries and snap fraction sum to 1:
    d = {}
    for fr,nme in ll:
        if nme not in d:
            d[nme] = [fr]
        else:
            d[nme] += [fr]
    ll = []
    for mfr, nme in sorted( (-math.fsum(frs)/fractot,nme) for nme,frs in d.items()):
        if mfr:
            ll.append( ( -mfr, nme ) )
    return label, ll

def determine_labels_and_atomdb( _self_params, *, fractions, allow_siteoccu_ncmatv5_hack = True ):
    compositions = _self_params.get('compositions',{})
    #compositions = copy.deepcopy(compositions)#we will modify a bit below and need to retain immutability
    atomdb = _self_params.get('atomdb',{})
    atompos = _self_params.get('atompos',None)
    occumap = atompos['occumap'] if atompos else None
    dyninfos = _self_params.get('dyninfos',{})
    param_fallback_debye_temp = _self_params.get('fallback_debye_temp',None)

    #labels:
    dyninfo_lbls = set( dyninfos.keys() )
    atompos_lbls = set( lbl for lbl,x,y,z in atompos['pos'] ) if atompos else set()
    fractions_lbls = set( (fractions or {}).keys() )
    direct_lbls = dyninfo_lbls.union( atompos_lbls ).union( fractions_lbls )

    #fill expand implicit compositions to be included explicitly:
    def _expand_implicit_compositions( direct_lbls, compositions ):
        for lbl,_ in compositions.items():
            if lbl not in direct_lbls:
                raise _nc_core.NCBadInput('Label "%s" was given a composition but is not actually used in the material'%lbl)

        compositions = copy.deepcopy(compositions)#retain immutability
        for lbl in direct_lbls:
            if lbl not in compositions:
                _ = _nc_common.check_elem_or_isotope_marker( lbl )
                if not _:
                    raise _nc_core.NCBadInput('Not able to determine composition associated with label '
                                     f'"{lbl}" since it was not an element or isotope. Please add a definition with .set_composition(..).')
                compositions[lbl] = [(1.0,_)]
        return compositions
    compositions = _expand_implicit_compositions( direct_lbls, compositions )

    #check completeness of labels in other sections:
    if atompos:
        if (dyninfo_lbls - atompos_lbls):
            raise _nc_core.NCBadInput('Some atoms with dynamic information are missing atomic positions: "%s"'%('", "'.join(dyninfo_lbls - atompos_lbls)))
        if not param_fallback_debye_temp and (atompos_lbls-dyninfo_lbls):
            raise _nc_core.NCBadInput('Missing dynamic information for some atoms: "%s"'%('", "'.join(atompos_lbls-dyninfo_lbls)))
    else:
        if (dyninfo_lbls-fractions_lbls):
            raise _nc_core.NCBadInput('Some atoms with dynamic information are missing fractions: "%s"'%('", "'.join(dyninfo_lbls - fractions_lbls)))

    lbls_with_nonunit_occu = set( lbl for lbl,occu in occumap.items() if occu != 1.0 ) if occumap else set()

    #SPECIAL HACK BEGIN
    if lbls_with_nonunit_occu and not allow_siteoccu_ncmatv5_hack:
        raise _nc_core.NCBadInput('Can not support site occupancies != 1.0 without the special hack for NCMAT v5 (and allow_siteoccu_ncmatv5_hack was set to False)')

    if lbls_with_nonunit_occu:
        #Handle via (non-persistent) modifications of compositions + atomdb:
        atomdb = copy.deepcopy(atomdb)
    atomdb_hack_comments = []
    for i,lbl in enumerate(sorted(lbls_with_nonunit_occu)):
        origelem_mass = calc_mass(lbl,compositions[lbl])
        compos = compositions[lbl][:]
        origelem_lbl = lbl if compos is None else _nc_common.format_chemform(list((v,k) for k,v in compos))
        assert 0.0 < origelem_mass < 2000.0
        #In principle I should check this is not already used, but if anyone
        #actually needs Og999 or thereabouts, they are obviously just trying
        #to break our hack by being smartasses :-)
        hijackedIsotope = 'Og%i'%(299-i)
        assert hijackedIsotope not in atomdb
        atomdb[hijackedIsotope] = _nc_core.AtomData.fmt_atomdb_str( origelem_mass, 0.0, 0.0, 0.0 )
        occu = occumap[lbl]
        compositions[ lbl ] = [ (occu*frac,elem) for frac,elem in compos ] + [ ( 1.0 - occu, hijackedIsotope ) ]
        atomdb_hack_comments.append ( f'#Note: Emulating site_occupancy={occu:g} for {origelem_lbl} by hijacking {hijackedIsotope} to play the role as "sterile {origelem_lbl}".' )
    #SPECIAL HACK END

    atomdb_lines = []
    atomdb_lines += atomdb_hack_comments
    for lbl,paramstr in sorted( atomdb.items() ):
        atomdb_lines.append('%s %s'%(lbl,paramstr))

    #We must remap exactly those labels which are not a single element or
    #isotopes, or those single/element isotopes that have more than one
    #role. So first we fill the composition for ALL:
    lblmap = {}
    Xidx_next = 1
    for lbl in sorted(direct_lbls):
        compos = compositions[lbl]
        if len(compos) == 1 and sum(1 for lbl,e in compositions.items() if (len(e)==1 and e[0][1]==compos[0][1]) ) == 1:
            lblmap[ lbl ] = compos[0][1]
            continue
        #must use special label X1..X99/X:
        if Xidx_next>100:
            raise _nc_core.NCBadInput('Material requires more than 100 special labels which is'
                             +' not supported by NCMAT (allowed are only X1,X2,...,X99 and X')
        speciallbl = 'X' if Xidx_next==100 else 'X%i'%Xidx_next
        Xidx_next += 1
        lblmap[lbl] = speciallbl
        if len(compos)==1:
            assert compos[0][0]==1.0
            atomdb_lines.append( '%s is %s'%(speciallbl,compos[0][1]))
        else:
            s=''
            for frac,name in compos:
                s += f' {frac:.12g} {name}'
            atomdb_lines.append( '%s is%s'%(speciallbl,s))
    return lblmap,atomdb_lines

def _is_nonempty_array( x ):
    return hasattr(x,'__len__') and len(x)>0

def _fmtprecisenum( v, strip_leading_zero = True ):
    _ = f'{v:.14g}' if v else '0'
    return _[1:] if (strip_leading_zero and _.startswith('0.')) else _

def _copyarray_or_None( x ):
    if x is None or (hasattr(x,'__len__') and len(x)==0):
        return None
    from ._numpy import _np
    return copy.deepcopy(_np.asarray(x,dtype=float) if _np else x)

def _decode_update_atomdb(element_or_isotope, data, *,  mass=None, coh_scat_len=None, incoh_xs=None, abs_xs=None ):
    e = _nc_common.check_elem_or_isotope_marker( element_or_isotope )
    if not e:
        raise _nc_core.NCBadInput(f'Invalid element or isotope label (must be of form "Al", "O16", "D", ...): {e}')
    label = e

    nphys = int(bool(mass))+int(bool(coh_scat_len))+int(bool(incoh_xs))+int(bool(abs_xs))
    hasdata = data is not None
    if not hasdata and not nphys:
        raise _nc_core.NCBadInput('Missing data parameter(s).')
    if (hasdata and nphys) or nphys not in (0,4):
        raise _nc_core.NCBadInput('Inconsistent set of parameters provided.')
    if hasdata:
        if hasattr(data,'to_atomdb_str'):
            #Handle AtomData objects (and others with to_atomdb_str method)
            res = data.to_atomdb_str()
        else:
            res = ' '.join(str(data).replace(':',' ').strip().split())
    else:
        assert nphys==4
        res = _nc_core.AtomData.fmt_atomdb_str( mass = mass,
                                                coh_scat_len = coh_scat_len,
                                                incoh_xs = incoh_xs,
                                                abs_xs = abs_xs )
    assert res and 'fm ' in res and 'u ' in res and 'b ' in res
    return label,res

def _lattice_params_to_vectors( a, b, c, alpha, beta, gamma ):
    if max(max(alpha,beta),gamma)<=math.pi:
        raise _nc_core.NCBadInput("angles are expected in degrees not radians")
    def _sincos( x ):
        if x == 0:
            return 0.0, 1.0
        if x == 90.0:
            return 1.0, 0.0
        if x == 120.0:
            return 0.86602540378443864676372317075293618347140262690519, -0.5
        return math.sin( x*(math.pi/180) ), math.cos( x*(math.pi/180) )
    sa,ca = _sincos( alpha )
    sb,cb = _sincos( beta )
    sg,cg = _sincos( gamma )
    assert sg > 0.0
    vect_a = tuple( a*e for e in ( 1, 0.0, 0.0 ) )
    vect_b = tuple( b*e for e in ( cg, sg, 0.0 ) )
    ux = cb
    assert sg > 0.0
    assert sg**2 > 0.0
    uy = (ca-cb*cg)/sg
    uzsq = sb**2 - uy**2
    if uzsq < -1e-99:
        raise _nc_core.NCBadInput('Inconsistent cell parameters detected')
    vect_c = tuple( c*e for e in ( ux, uy, math.sqrt(max(0.0,uzsq)) ) )#sign of third component positive since we want the three vectors to form a right handed basis.
    return ( vect_a, vect_b, vect_c )

def _lattice_vectors_to_params( vect_a, vect_b, vect_c, as_dict = True ):
    def mag( v ):
        return math.sqrt( math.fsum(e**2 for e in v) )
    def dot( v1,v2 ):
        return math.fsum(e1*e2 for e1,e2 in zip(v1,v2))
    def ang( v1,v2 ):
        return math.acos( dot(v1,v2)/(mag(v1)*mag(v2)) )*180/math.pi
    a,b,c = mag(vect_a), mag(vect_b), mag(vect_c)
    alpha,beta,gamma = ang(vect_b,vect_c), ang(vect_a,vect_c), ang(vect_a,vect_b)
    if as_dict:
        return dict( a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma )
    else:
        return a, b, c, alpha, beta, gamma

def _snap_lattice_params( params_dict, relative_tolerance = 1e-6 ):
    log = []
    p = copy.deepcopy(params_dict)
    def rd( x,y ):
        return abs(x-y)/(max(1e-300,abs(x)+abs(y)))
    super_tolerance = 1e-14#almost 64 bit FP prec
    assert relative_tolerance > super_tolerance
    def is_near( x,y ):
        return rd(x,y) < relative_tolerance
    def is_super_near( x,y ):
        return rd(x,y) < super_tolerance

    #"snap" b,c to a,b and angles to 60/90/120:
    for l1,l2 in (('a','b'),('a','c'),('b','c')):
        if is_near(p[l1],p[l2]):
            if not is_super_near(p[l1],p[l2]):
                log.append(f'Snapping lattice parameter {l2} to {l1} since it was very near.')
            p[l2] = p[l1]
    for commonangle in (45,60,90,120):
        for angname in ['alpha','beta','gamma']:
            if is_near(p[angname],commonangle):
                if not is_super_near(p[angname],commonangle):
                    log.append(f'Snapping lattice angle {angname} to {commonangle} since it was very near.')
                p[angname] = commonangle
    return log, p

def _cellparams_to_spglib_lattice( cellsg ):
    return _lattice_params_to_vectors( cellsg['a'], cellsg['b'], cellsg['c'],
                                       cellsg['alpha'], cellsg['beta'], cellsg['gamma'] )

def _reldiff_cellparams( c1, c2 ):
    #compares a,b,c,alpha,beta,gamma in the two cellsg arrays.
    def rd( x,y ):
        return abs(x-y)/(max(1e-300,abs(x)+abs(y)))
    def rdf( s ):
        return rd(c1[s],c2[s])
    return max( rdf(s) for s in ('a','b','c','alpha','beta','gamma') )

def _reldiff_atompos( spglib_cell1, spglib_cell2 ):
    #input in spglib cell format, returns largest dist between positions, or
    #None if not the same type+number of atoms.

    for c in ( spglib_cell1, spglib_cell2 ):
        assert c is not None
        assert len(c) >= 3
        assert c[0] is not None
    _,allpos1,idxlist1 = spglib_cell1
    _,allpos2,idxlist2 = spglib_cell2
    assert len(idxlist1) == len(allpos1)
    assert len(idxlist2) == len(allpos2)
    if list(sorted(idxlist1)) != list(sorted(idxlist2)) or len(idxlist1)!=len(idxlist2):
        return None

    def calc_maxdistsq( l1, l2 ):
        n = len(l1)
        assert n == len(l2)
        distsqmax = 0.0
        for e1 in l1:
            dsqmin_e1 = None
            for e2 in l2:
                _ = _unit_cell_point_distsq(e1,e2)
                if dsqmin_e1 is None or _ < dsqmin_e1:
                    dsqmin_e1 = _
                    if _ == 0.0:
                        break
            assert dsqmin_e1 is not None
            distsqmax = max( dsqmin_e1, distsqmax )
        return distsqmax

    max_distsq_seen = None
    for idx in sorted(set(idxlist1)):
        def fixp( p ):
            return (_remap_fract_pos(p[0]),
                    _remap_fract_pos(p[1]),
                    _remap_fract_pos(p[2]))
        l1 = list( sorted( fixp(p) for i,p in zip(idxlist1,allpos1) if i == idx ) )#nb: fragile sort order!
        l2 = list( sorted( fixp(p) for i,p in zip(idxlist2,allpos2) if i == idx ) )#nb: fragile sort order!
        d = calc_maxdistsq( l1, l2 )
        if d is None:
            return None
        if max_distsq_seen is None or d > max_distsq_seen:
            max_distsq_seen = d
    return math.sqrt(max_distsq_seen) if max_distsq_seen is not None else None

def _import_spglib( *, sysexit = False ):
    try:
        import spglib#both available on pypi and conda-forge
    except ImportError:
        m = ( 'Could not import spglib module needed to standardise and verify crystal structures.'
              +' The spglib package is available on both PyPI ("python3 -mpip install'
              +' spglib") and conda ("conda install -c conda-forge spglib")' )
        if sysexit:
            raise SystemExit(m)
        else:
            raise ImportError(m)
    return spglib.spglib

def _import_ase( *, sysexit = False ):
    try:
        import ase#both available on pypi and conda-forge
        import ase.io
    except ImportError:
        m = ( 'Could not import ase module.'
              +' The ase package is available on both PyPI ("python3 -mpip install'
              +' ase") and conda ("conda install -c conda-forge ase")' )
        if sysexit:
            raise SystemExit(m)
        else:
            raise ImportError(m)
    return ase, ase.io

def _spglib_extractsg( spglib_symdata ):
    sd = spglib_symdata
    sgno = ( getattr(sd,'number')
             if hasattr(sd,'number')
             else sd['number'] )
    sgsymb_hermann_mauguin = ( getattr(sd,'international')
                               if hasattr(sd,'international')
                               else sd['international'] )
    c = ( getattr(sd,'choice')
          if hasattr(sd,'choice')
          else sd.get('choice',None))
    if c:
        sgsymb_hermann_mauguin += f':{c}'
    return sgno, sgsymb_hermann_mauguin

def _check_cell_sg_consistency( sgnumber, a, b, c, alpha, beta, gamma ):
    sgclass = _nc_common._classifySG(int(sgnumber))
    if sgclass in ( 'orthorhombic', 'tetragonal', 'cubic'):
        if sgclass=='cubic':
            if a!=b or a!=c:
                raise _nc_core.NCBadInput(f'Cubic space group {sgnumber} requires cell parameters a==b==c (found a={a:g}, b={b:g}, c={c:g})')
        if not all( (_ == 90) for _ in (alpha, beta, gamma) ):
            raise _nc_core.NCBadInput(f'Space group {sgnumber} requires cell parameters alpha=beta=gamma=90 '
                                +f'(found alpha={alpha}, beta={beta}, gamma={gamma})')
    elif sgclass in ('trigonal','hexagonal'):
        if not ( alpha == 90 and beta == 90 and gamma == 120 ):
            raise _nc_core.NCBadInput(f'Space group {sgnumber} requires cell parameters alpha=beta=90 and gamma=120 '
                                +f'(found alpha={alpha}, beta={beta}, gamma={gamma})')
    else:
        if not all( (0<a<180) for a in ( alpha, beta, gamma ) ):
            raise _nc_core.NCBadInput('Unit cell must have all angles between 0 and 180'
                                +f' (found alpha={alpha}, beta={beta}, gamma={gamma})')
    if int(sgnumber) >= 75:
        if a!=b:
            raise _nc_core.NCBadInput(f'Space group {sgnumber} requires cell parameters a==b (found a={a:g}, b={b:g}, c={c:g})')

def _spglib_refine_cell( spglib_cell, symprec = 0.01 ):
    assert spglib_cell
    p_orig = _lattice_vectors_to_params( *spglib_cell[0], as_dict = True )

    warnings, msgs = [], []
    refined_cell, symdata = _impl_spglib_refine_cell( spglib_cell, symprec=symprec, warnings=warnings, msgs=msgs )

    #check difference from first to final, to warn/msg about corrections:
    p_final = _lattice_vectors_to_params( *refined_cell[0], as_dict = True )
    _,p_final = _snap_lattice_params( p_final )
    sgno, sgsymb_hm = _spglib_extractsg( symdata )
    rdl = _reldiff_cellparams( p_orig, p_final )
    rda = _reldiff_atompos( spglib_cell, refined_cell )
    rd = None if (rdl is None or rda is None) else max(rdl,rda)
    discard_aniso = False
    if rd is None or rd > 1e-2:
        discard_aniso = True
        _ = 'possibly due to conversion from primitive cell' if rd is None else f'at the {100.0*rd:g}% level'
        warnings.append(f'Structure received major corrections by spglib ({_})')
    elif rd > 1e-6:
        warnings.append(f'Structure received minor corrections by spglib (at the {100.0*rd:g}% level)')
    else:
        msgs.append('Self-consistency of structure was verified by spglib')

    #remove all but the last of any kind of error message (otherwise we might get many repetitions of 'snapping b=a' kind of messages, due to the recursion.

    warnrev = []
    for w in warnings[::-1]:
        if w not in warnrev:
            warnrev.append(w)
    warnings = warnrev[::-1]

    return dict( refined_cell = refined_cell,
                 cellparams_snapped = p_final,
                 sgno = sgno,
                 sgsymb_hm = sgsymb_hm,
                 warnings = warnings,
                 msgs = msgs,
                 can_keep_anisotropic_properties = not discard_aniso )

def _impl_spglib_refine_cell( spglib_cell, *, symprec, warnings, msgs, nrepeat = 0 ):
    spglib = _import_spglib()

    orig_cell = spglib_cell

    orig_cell_copy = copy.deepcopy( orig_cell )#copy is just a safeguard
    refined_cell = spglib.standardize_cell( orig_cell_copy, symprec = symprec )
    symdata = spglib.get_symmetry_dataset( orig_cell_copy )
    if not refined_cell or not symdata or refined_cell[0] is None:
        raise _nc_core.NCBadInput('Could not standardize provided crystal structure with spglib')

    p = _lattice_vectors_to_params( *refined_cell[0], as_dict = True )

    snaplog, p = _snap_lattice_params( p )
    if snaplog or nrepeat == 0:
        #This is either the first try, of we snapped some parameter
        #(e.g. alpha=89.99999 -> alpha=90.0). After this a second refinement
        #might find a different (higher) symmetry:
        warnings += snaplog[:]
        refined_cell = tuple( [ _cellparams_to_spglib_lattice( p ) ] + list( e for e in refined_cell[1:] ) )
        if nrepeat <= 3:
            return _impl_spglib_refine_cell( refined_cell, symprec = symprec, warnings=warnings, msgs=msgs, nrepeat = nrepeat + 1 )

    return refined_cell, symdata


def _remap_fract_pos_pt(xyz):
    return ( _remap_fract_pos(xyz[0]),
             _remap_fract_pos(xyz[1]),
             _remap_fract_pos(xyz[2]) )

def _remap_fract_pos(x):
    x = float(x)
    if -1.0<=x<0.0:
        x += 1.0
    if 1.0<=x<2.0:
        x -= 1.0
    return x

def _unit_cell_point_dist( p1, p2 ):
    return math.sqrt( _unit_cell_point_distsq( p1, p2 ) )

def _unit_cell_point_distsq( p1, p2 ):
    #take into account unit cell periodicity, so e.g. (a,b,epsilon) and
    #(a,b,1.0-epsilon) are actually only 2*epsilon from each other, not
    #1.0-epsilon.
    def remap( x ):
        return x-math.floor(x)
    d = [ abs(remap(e1)-remap(e2)) for e1,e2 in zip(p1,p2)]
    return math.fsum(  min(e,1.0-e)**2 for e in d )

def _cifdata_via_ase( data_or_file, ase_format = None, quiet = False ):
    import io
    ase,ase_io = _import_ase()
    ase_obj = data_or_file if (ase_format=='ase' or NCMATComposerImpl._is_ase_like_object( data_or_file )) else None
    prfct = ( lambda *a, **kwa : None ) if quiet else _nc_common.print

    if not ase_obj and hasattr( data_or_file, 'load_data' ) and hasattr( data_or_file, 'codid' ) and hasattr( data_or_file, 'mpid' ):
        #assume object is CIFSource object:
        prfct('Attempting to load data via ASE')
        ase_obj = ase_io.read( io.StringIO( data_or_file.load_data( quiet = quiet) ),
                               format = ase_format or 'cif' )

    isbytes = hasattr( data_or_file, 'decode' )
    str_to_sb = ( lambda s : s ) if not isbytes else ( lambda s : s.encode() )
    if not ase_obj and ( hasattr(data_or_file,'__fspath__')
                         or str_to_sb('\n') not in data_or_file ):
        #on-disk data:
        _ = 'data'
        if ( hasattr(data_or_file,'__fspath__') or isinstance(data_or_file,str) ):
            import pathlib
            _ = pathlib.Path(data_or_file).name
        prfct(f'Attempting to load {_} via ASE')
        ase_obj = ase_io.read( data_or_file, format = ase_format )
    if not ase_obj:
        #fall back to in-mem data
        if isbytes:
            dof_bytes = data_or_file
            try:
                dof_str = data_or_file.decode('utf8')
            except UnicodeDecodeError:
                dof_str = None
        else:
            dof_str = data_or_file
            try:
                dof_bytes = data_or_file.encode('utf8')
            except UnicodeEncodeError:
                dof_bytes = None

        if not ase_format:
            import ase.io.formats
            if dof_bytes is not None:
                try:
                    ase_format = ase.io.formats.match_magic( dof_bytes ).name
                except ase.io.formats.UnknownFileTypeError:
                    ase_format = None
            if ase_format is None and dof_str is not None and '\n' not in dof_str:
                try:
                    ase_format = ase.io.formats.filetype( dof_str ).name
                except ase.io.formats.UnknownFileTypeError:
                    ase_format = None
            if not ase_format:
                _nc_common.warn('Could not determine ASE format of'
                                ' input data. Assuming CIF.')
                ase_format = 'cif'
        assert isinstance(ase_format,str)
        prfct('Attempting to load data via ASE')
        if dof_str:
            assert isinstance(dof_str,str)
            with io.StringIO(dof_str) as memfile:
                ase_obj = ase_io.read( memfile, format = ase_format )
        else:
            assert isinstance(dof_bytes,bytes)
            with io.BytesIO(dof_bytes) as memfile:
                ase_obj = ase_io.read( memfile, format = ase_format )

    assert ase_obj is not None

    with io.BytesIO() as s:
        ase_io.write(s, ase_obj, format='cif')
        return s.getvalue().decode()

def formatVectorForNCMAT(name,values,indent):
    def provideFormattedEntries():
        def _fmtnum(num):
            _ = '%g'%num if num else '0'#avoid 0.0, -0, etc.
            if _.startswith('0.'):
                _=_[1:]
            return _
        from ._numpy import _np
        v = _np.asarray(values,dtype=float).flatten() if _np else values
        i, nv = 0, len(v)
        while i < nv:
            fmt_vi=_fmtnum(v[i])
            #check if is repeated:
            irepeat=i
            while irepeat+1<nv:
                if _fmtnum(v[irepeat+1])==fmt_vi:
                    irepeat+=1
                else:
                    break
            yield '%sr%i'%(fmt_vi,1+irepeat-i) if irepeat>i else '%s'%fmt_vi
            i=irepeat+1#advance
    out=''
    line='  %s'%name
    collim=80
    for e in provideFormattedEntries():
        snext=' %s'%e
        line_next=line+snext
        if len(line_next)>collim:
            out += line
            out += '\n'
            line = indent+snext
        else:
            line = line_next
    if line:
        out += line
        out += '\n'
    return out

def _stripCommentsAndPackNCMATData( ncmat_data ):
    if hasattr( ncmat_data, 'rawData' ):
        return _stripCommentsAndPackNCMATData( ncmat_data.rawData )
    if not ncmat_data.startswith('NCMAT'):
        return ncmat_data#pass through broken data unchanged!
    out = ['']
    def addline( s ):
        if s:
            out[0] += ( s + '\n' )
    for e in ncmat_data.splitlines():
        p = e.split('#',1)
        if not p:
            continue
        p0 = ' '.join(p[0].split())
        if len(p)==2 and 'NCRYSTALMATCFG[' in p[1]:
            c = p[1].split('NCRYSTALMATCFG[',1)
            if len(c)==2:
                i = c[1].rfind(']')
                if i == -1:
                    #seems to be a syntax error, just keep everything:
                    addline(p0 + '#' + p[1])
                else:
                    addline(p0 + '#NCRYSTALMATCFG[' + c[1][0:i].strip() + ']')
        else:
            addline(p0)
    return out[0]

def _extractInitialHeaderCommentsFromNCMATData( ncmat_data, dedent = True ):
    #Extracts initial comments, dedented appropriately. Returns None if data
    #does not look like NCMAT data (an empty list returned is not a failure, it
    #just means no comments).
    from .misc import AnyTextData
    td = AnyTextData(ncmat_data)
    if '\n' not in td.content and not td.content.startswith('NCMAT'):
        return None
    comments = []
    first = True
    for line in td.content.splitlines():
        if first:
            first = False
            continue
        if line.strip().startswith('@'):
            break
        line=line.lstrip()
        if line.startswith('#'):
            comments.append(line[1:].rstrip())
    #dedent and return:
    if dedent:
        while all( ( e.startswith(' ') or not e) for e in comments ):
            comments = [ e[1:] for e in comments ]
    return comments

def _check_valid_custom_section_name( name ):
    if ( not name
         or not isinstance(name,str)
         or not name.isalpha()
         or not name.isupper() ):
        raise _nc_core.NCBadInput(f'Invalid custom section name: "{name}" '
                                  '(must be non-empty and contain only'
                                  ' capitalised letters A-Z)')
