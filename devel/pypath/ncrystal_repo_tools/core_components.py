
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2026 NCrystal developers                                   ##
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

from .dirs import coreroot as _coreroot
_srcroot = _coreroot / 'src'
_incroot = _coreroot / 'include'

def is_valid_component_name( s ):
    return ( 3 <= len(s) <= 14
             and s==s.lower()
             and s.isidentifier()
             and s[0].isalpha()
             and s[-1].isalpha() )

class Component:
    def __lt__( self, o ):
        return self.name < o.name

    def __init__(self, name ):
        if not is_valid_component_name(name):
            raise SystemExit(f'Not a valid component name: "{name}"')
        self.name = name
        self.srcdir = _srcroot / name
        assert self.srcdir.is_dir(), f"Not a directory: {self.srcdir}"
        self.depfile = self.srcdir/ 'dep.txt'
        assert self.depfile.exists(), f"not found: {self.depfile}"
        def globsafe( path, pattern ):
            for e in path.glob(pattern):
                if '#' not in e.name and '~' not in e.name:
                    yield e
        self.srcfiles = tuple(sorted(globsafe(self.srcdir,'*.cc')))
        self.local_hdrs = tuple(sorted(globsafe(self.srcdir,'*.hh')))
        pub_hdrdir = _incroot.joinpath('NCrystal',name)
        internal_hdrdir = _incroot.joinpath('NCrystal','internal',name)
        if pub_hdrdir.exists() and internal_hdrdir.exists():
            raise RuntimeError('Can not have both '
                               f'{pub_hdrdir} and {internal_hdrdir}')
        if pub_hdrdir.exists():
            assert pub_hdrdir.is_dir()
            hdrdir = pub_hdrdir
            self.is_internal = False
        elif internal_hdrdir.exists():
            assert internal_hdrdir.is_dir()
            hdrdir = internal_hdrdir
            self.is_internal = True
        else:
            hdrdir = None
            self.is_internal = True
        self.hdrdir = hdrdir
        if hdrdir is not None:
            self.hdrfiles = tuple(sorted(list(globsafe(hdrdir,'*.h'))
                                         +list(globsafe(hdrdir,'*.hh'))))
            if not self.hdrfiles:
                raise SystemExit(f'ERROR empty dir: {hdrdir}')
        else:
            self.hdrfiles = tuple([])

        self.hdrfiles_icc = tuple([])#TODO: support?

        self.direct_depnames = sorted(set(self.depfile.read_text().split()))
        #To be filled later:
        self.deps = None
        self.depnames = None
        self.direct_deps = None
        _ispref = ( f'NCrystal/internal/{self.name}/'
                    if self.is_internal
                    else f'NCrystal/{self.name}/' )
        self.__incstatements_to_pkg = set( _ispref + h.name
                                           for h in
                                           sorted( self.hdrfiles ) )

    def _init_deps( self, name2comp, block = None ):
        if self.deps is not None:
            return#done
        if block is None:
            block = set()
        else:
            if self.name in block:
                raise RuntimeError('Circular dependency involving'
                                   f' component "{self.name}"')
        block.add(self.name)
        depnames = set()
        deps = []
        direct_deps = []
        for dn in self.direct_depnames:
            dc = name2comp.get(dn)
            if dc is None:
                raise RuntimeError(f'Unknown dependency "{dn}" in '
                                   f'component "{self.name}"')
            direct_deps.append( dc )
            dc._init_deps( name2comp, block )
            if dn not in depnames:
                depnames.add( dn )
                deps.append( dc )
            for dc2 in dc.deps:
                if dc2.name not in depnames:
                    depnames.add( dc2.name )
                    deps.append( dc2 )
        self.depnames = depnames
        self.deps = deps
        self.direct_deps = set(direct_deps)
        self._incstm_in_pkgincdir = None
        self._incstm_in_pkgsrcdir = None

    def calc_minimal_deps( self ):
        #direct_deps with those deps provided indirectly removed
        indirect_deps = set()
        for d in self.direct_deps:
            indirect_deps.update( d.deps )
        return self.direct_deps - indirect_deps

    def allowed_incs_in_pkgincdir( self ):
        if self._incstm_in_pkgincdir is not None:
            return self._incstm_in_pkgincdir
        res = set()
        #Always our own of course:
        res.update( self.__incstatements_to_pkg )
        for d in self.deps:
            if not self.is_internal and d.is_internal:
                #A non-internal pkg header can not include internal headers,
                #because that would leak internal headers to public api!
                continue
            res.update( d.__incstatements_to_pkg )
        self._incstm_in_pkgincdir = res
        return res

    def allowed_incs_in_pkgsrcdir( self ):
        if self._incstm_in_pkgsrcdir is not None:
            return self._incstm_in_pkgsrcdir
        res = set()
        #Always our own of course:
        res.update( self.__incstatements_to_pkg )
        #From all dep pkgs as well:
        for d in self.deps:
            res.update( d.__incstatements_to_pkg )
        #And finally the local headers:
        for h in ( self.local_hdrs or [] ):
            res.add( h.name )
        self._incstm_in_pkgsrcdir = res
        return res

    def all_file_iter( self ):
        for filelist in [self.hdrfiles,self.srcfiles,self.local_hdrs]:
            for f in filelist or []:
                yield f

    def sloc_count( self, headers_only = False ):
        files = self.hdrfiles if headers_only else self.all_file_iter()
        return sum( file_calc_sloc_count(f) for f in files )

def file_calc_sloc_count( f ):
    #TODO: Improve
    return len(f.read_text().splitlines())

def load_components( *, init_deps = True ):
    name2comp = dict( (d.name, Component( d.name ))
                      for d in _srcroot.iterdir()
                      if d.is_dir() and '.' not in d.name )
    if init_deps:
        for n,c in name2comp.items():
            c._init_deps( name2comp )
    return name2comp

def fix_deps( dryrun = False ):
    from .extract_includes import get_include_staments_from_file
    for c in load_components( init_deps = False ).values():
        deplist = set()
        for f in c.all_file_iter():
            for i in get_include_staments_from_file(f):
                if not i.startswith('NCrystal/'):
                    continue
                p = i.split('/')
                if len(p) == 4 and i.startswith('NCrystal/internal/'):
                    deplist.add( p[2] )
                elif len(p) == 3 and i.startswith('NCrystal/'):
                    deplist.add( p[1] )
        deptxt = '\n'.join(sorted(e for e in deplist if e != c.name)) + '\n'
        if c.depfile.read_text() != deptxt:
            if dryrun:
                print("Would update",c.depfile)
            else:
                print("Updating",c.depfile)
                c.depfile.write_text( deptxt )
