
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2024 NCrystal developers                                   ##
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

import common.dirs
srcroot = common.dirs.coreroot / 'src'
incroot = common.dirs.coreroot / 'include/NCrystal'

class Component:
    def __lt__( self, o ):
        return self.name < o.name

    def __init__(self, name ):
        self.name = name
        self.srcdir = srcroot / name
        assert self.srcdir.is_dir(), f"Not a directory: {self.srcdir}"
        depfile = self.srcdir/ 'dep.txt'
        assert depfile.exists(), f"not found: {depfile}"
        self.srcfiles = sorted(self.srcdir.glob('*.cc'))#ok with just a dep.txt and
                                                        #no actual src files
        self.srcdir_hdrs = sorted(self.srcdir.glob('*.hh'))
        pub_hdrdir = incroot / name
        internal_hdrdir = incroot / 'internal' / name
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
        if hdrdir is not None:
            self.hdrfiles = sorted(list((hdrdir).glob('*.h'))
                                   +list((hdrdir).glob('*.hh')))
            assert self.hdrfiles, f"empty dir: {hdrdir}"
        else:
            self.hdrfiles = None

        self.direct_depnames = sorted(set(depfile.read_text().split()))
        #To be filled later:
        self.deps = None
        self.depnames = None

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
        for dn in self.direct_depnames:
            dc = name2comp.get(dn)
            if dc is None:
                raise RuntimeError(f'Unknown dependency "{dn}" in '
                                   f'component "{self.name}"')
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

def load_components():
    name2comp = dict( (d.name ,Component( d.name ))
                      for d in srcroot.iterdir()
                      if d.is_dir() )
    for n,c in name2comp.items():
        c._init_deps( name2comp )
    return name2comp
