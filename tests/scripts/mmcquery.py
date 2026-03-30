#!/usr/bin/env python3

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

# NEEDS: numpy

from NCrystalDev.misc import evaluate_query as ncquery
import numpy # noqa F401 (before fpe)
import NCTestUtils.enable_fpe # noqa F401
from NCrystalDev.core import NCBadInput
from NCTestUtils.common import ensure_error
import pprint

def test(*query):
    print()
    print(">>> Sending query: %s"%repr(list(e for e in query)))
    r = ncquery(query)
    pprint.pp(r)

def test_examples(subject):
    assert subject in ('geom','src','srcenergy')
    is_cpe = False
    if subject=='srcenergy':
        is_cpe = True
        subject = 'src'


    exlist = ncquery(['mmc','cfgdoc',subject])['commonpars_energy_examples'
                                               if is_cpe else 'examples']
    for exval, descr in exlist:
        if is_cpe:
            exval=f'constant;{exval}'
        test('mmc','inspectcfg',subject,exval)

def main():

    with ensure_error(NCBadInput,
                      'Missing or empty query'):
        test()
    with ensure_error(NCBadInput,
                      'Invalid JSON query key: "foo"'):
        test('foo')

    with ensure_error(NCBadInput,
                      'Invalid JSON query key: ""'):
        test('')
    with ensure_error(NCBadInput,
                      'Invalid MiniMC JSON query: ["mmc",""]'):
        test('mmc','')

    with ensure_error(NCBadInput,
                      'Query items can not start with a "-" character.'
                      ' Problem found in: "-bla"'):
        test('mmc','-bla')
    with ensure_error(NCBadInput,
                      'Query items can not start with a "-" character.'
                      ' Problem found in: "-bla"'):
        test('mmc','-bla')
    with ensure_error(NCBadInput,
                      'Query items can not start with a "-" character.'
                      ' Problem found in: " -bla"'):
        test('mmc',' -bla')
    with ensure_error(NCBadInput,
                      'Invalid MiniMC JSON query: ["mmc","bla-bla"]'):
        test('mmc','bla-bla')

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

    test('mmc','inspectcfg','src','constant; z=-10;wl=2.0')
    test('mmc','inspectcfg','src','isotropic; z=-10;ekin=0.025')
    test('mmc','inspectcfg','src','circular; r=0.1;z=-1;wl=1.8;ux=1')

    test('mmc','inspectcfg','src','constant; z=-10;wl=2.0-3.0')
    test('mmc','inspectcfg','src','constant; z=-10;ekin=0.025-1.0')
    test('mmc','inspectcfg','src','constant; z=-10;wl=2.0+-0.1')
    test('mmc','inspectcfg','src','constant; z=-10;ekin=0.025+-0.001')
    test('mmc','inspectcfg','src','constant; z=-10;ekin=thermal:50K')
    energy_helpstr="""
    Examples for how to set: "ekin=0.025" (fixed in eV), "wl=1.8" (fixed in Aa),
    "ekin=0.025-0.05" (uniform range in eV), "wl=1.8-5" (uniform range in Aa),
    "ekin=0.025+-0.001" (log-normal of given mean and rms in eV), "wl=1.8+-0.01"
    (log-normal of given mean and rms in Aa), "ekin=thermal:77" (thermal Maxwell
    of given temperature in K).
    """
    energy_helpstr=' '.join(energy_helpstr.split())
    with ensure_error(NCBadInput,
                      f'Invalid value "1eV" for parameter "ekin". {energy_helpstr}'):
        test('mmc','inspectcfg','src','constant; z=-10;ekin=1eV')
    with ensure_error(NCBadInput,
                      f'Invalid value "1Aa" for parameter "wl". {energy_helpstr}'):
        test('mmc','inspectcfg','src','constant; z=-10;wl=1Aa')
    with ensure_error(NCBadInput,energy_helpstr):
        test('mmc','inspectcfg','src','constant; z=-10;wl=help')
    with ensure_error(NCBadInput,energy_helpstr):
        test('mmc','inspectcfg','src','constant; z=-10;ekin=help')
    with ensure_error(NCBadInput,
                      'For a Maxwell distribution, do not use the "wl" para'
                      'meter (use "ekin" instead, like "ekin=thermal:300K").'):
        test('mmc','inspectcfg','src','constant; z=-10;wl=thermal:50K')

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
                      'Invalid parameter for chosen engine: "whatever"'):
        test('mmc','inspectcfg','engine','whatever = 0')
    with ensure_error(NCBadInput,
                      'Invalid tally "nobreakdown"'
                      ' (use "tallybreakdown=0" to disable breakdowns)'):
        test('mmc','inspectcfg','engine','tally=nobreakdown,theta,theta')

    test('mmc','inspectcfg','engine','tally=theta,theta;tallybreakdown=0')

    test('mmc','inspectcfg','engine','')
    test('mmc','inspectcfg','engine','nthreads=4')
    test('mmc','inspectcfg','engine','roulette=0.8,0.001,1')
    test('mmc','inspectcfg','engine','nthreads=0 ;   ignoremiss =1')
    test('mmc','inspectcfg','engine','std;nthreads=0;ignoremiss=1')
    test('mmc','inspectcfg','engine',
         '\t  std ; nthreads=  auto  \t\n;ignoremiss= 0')

    test('mmc','inspectcfg','engine','tally=q,theta;tallybins=+,;tallybreakdown=0')

    test('mmc','cfgdoc','geom')
    test('mmc','cfgdoc','src')
    test('mmc','cfgdoc','engine')
    test('mmc','cfgdoc','tally')

    test_examples('geom')
    test_examples('src')
    test_examples('srcenergy')


if __name__ == '__main__':
    main()
