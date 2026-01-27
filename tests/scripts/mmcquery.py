#!/usr/bin/env python3

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

# NEEDS: numpy

from NCrystalDev._common import json_query_cpplayer as ncquery
import numpy # noqa F401 (before fpe)
import NCTestUtils.enable_fpe # noqa F401
from NCrystalDev.core import NCBadInput
from NCTestUtils.common import ensure_error

def test(*query):
    print()
    print(">>> Sending query: %s"%repr(list(e for e in query)))
    r = ncquery(query)
    import pprint
    pprint.pprint(r)

def main():
    with ensure_error(NCBadInput,
                      'Missing or empty query'):
        test()
    with ensure_error(NCBadInput,
                      'Invalid JSON query key: "foo"'):
        test('foo')

    with ensure_error(NCBadInput,
                      "Invalid query (not all entries are strings):"
                      " ('foo', 117)"):
        test('foo',117)

    for a in ('bar'+'\x07','\x07','\x07\x07','\x07'+'bar'):
        with ensure_error(NCBadInput,
                          'Can not use character 0x07 (ASCII BEL) in JSON'
                          ' query strings'):
            test('foo',a)
    with ensure_error(NCBadInput,
                      'Invalid MiniMC JSON query: ["mmc","bla"]'):
        test('mmc','bla')
    with ensure_error(NCBadInput,
                      'Invalid MiniMC JSON query: ["mmc","scenario"] (correct'
                      ' usage: ["mmc","scenario",CFGSTR,SCENARIOSTR])'):
        test('mmc','scenario')
    with ensure_error(NCBadInput,
                      'Invalid MiniMC JSON query: ["mmc","scenario","1.8Aa"]'
                      ' (correct'
                      ' usage: ["mmc","scenario",CFGSTR,SCENARIOSTR])'):
        test('mmc','scenario','1.8Aa')


    test('mmc','scenario','Al_sg225.ncmat','1.8Aa')
    test('mmc','scenario','Al_sg225.ncmat','')

    with ensure_error(NCBadInput,
                      'Invalid MiniMC JSON query: ["mmc","inspectcfg"] (correct'
                      ' usage: ["mmc","inspectcfg","src|geom|engine",STRCFG])'):
        test('mmc','inspectcfg')
    with ensure_error(NCBadInput,
                      'Invalid MiniMC JSON query: ["mmc","inspectcfg","bla",""]'
                      ' (correct'
                      ' usage: ["mmc","inspectcfg","src|geom|engine",STRCFG])'):
        test('mmc','inspectcfg','bla','')



    with ensure_error(NCBadInput,
                      'Invalid src cfg: ""'):
        test('mmc','inspectcfg','src','')

    test('mmc','inspectcfg','engine','tally=cosmu,cosmu')

    test('mmc','inspectcfg','src','constant; z=-10;wl=2.0')
    test('mmc','inspectcfg','src','isotropic; z=-10;ekin=0.025')
    test('mmc','inspectcfg','src','circular; r=0.1;z=-1;wl=1.8;ux=1')

    #Test length formatting by testing various scale of thicknesses:
    test('mmc','inspectcfg','geom','slab;dz=1e-10')
    test('mmc','inspectcfg','geom','slab;dz=2e-6')
    test('mmc','inspectcfg','geom','slab;dz=2e-3')
    test('mmc','inspectcfg','geom','slab;dz=2e-2')
    test('mmc','inspectcfg','geom','slab;dz=0.2')
    test('mmc','inspectcfg','geom','slab;dz=10.2')
    test('mmc','inspectcfg','geom','slab;dz=2000')

    #Other shapes:
    test('mmc','inspectcfg','geom','box;dx=1;dy=0.1;dz=0.001')
    test('mmc','inspectcfg','geom','sphere;r=17.234')
    test('mmc','inspectcfg','geom','cyl;r=0.5')
    test('mmc','inspectcfg','geom','cyl;r=0.5;dy=12.0')

    with ensure_error(NCBadInput,
                      'Invalid MiniMC JSON query: ["mmc","tallylist","bla"]'
                      ' (no arguments should come after: ["mmc","tallylist"])'):
        test('mmc','tallylist','bla')
    test('mmc','tallylist')

if __name__ == '__main__':
    main()
