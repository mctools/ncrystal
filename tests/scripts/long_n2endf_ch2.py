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

# NEEDS: numpy endf-parserpy

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.reprint_escaped_warnings # noqa F401
from NCTestUtils.ncmat2endf_utils import test_cfg
import NCrystalDev._ncmat2endf_impl as ncmat2endf_impl
from NCrystalDev.ncmat2endf import EndfMetaData
from NCrystalDev.core import createTextData as nc_createTextData
import pathlib

ncmat2endf_impl.unit_test_chop_vals[0] = True

d = {'matnum':{"C":37, "H": 38},
     'edate':'JUL01',
     'ddate':'JUL01',
     'rdate':'JUL01',
     'alab':'TestLab',
     'libname':'TestLib',
     'auth': 'Jane Doe',
     'reference': 'Aaaaaa, et al.',
     'nlib': 0,
     'nver': 4
}
m = EndfMetaData()
m.update_from_dict(d)

metadata = EndfMetaData()
metadata.set_value( 'MATNUM', {"C":37, "H": 38} )

longfn = ('Polyethylene_as_an_amorphous_polymeric_material_consisting_of'
          '_a_long_chain_of_carbon_atoms_with_two_hydrogen_atoms_for_each'
          '_of_those_carbon_atoms_thus_yielding_the_formula_CH2.ncmat')

longstuff = ('# This is a very long comment line talking about all the'
             ' wonders of polymers and the authors childhood experiences'
             ' gazing at the sky and pondering polymers and all the'
             ' other wonders of material science, hoping to one day grow'
             ' up and work on Monte Carlo code evolving around such'
             ' marvelous materials. Also discussed is why underscores'
             ' are superior to spaces.')
testlines = """#
# some hdr
# ------------------------------------------------------------------------------
# bla
# ===== foo <many=>:
#
#  bla bla bla.
#
# bla
#
# <many-> foo <many=>:
#
#  bla bla bla.
#
# bla
#
# == foo <many=>:
#  bla bla bla.
#
#<many~>   :  \t
# yiha
# <many*>
#
# a lot of
# lines of text
# here
#
"""

testlines_expanded = []
for line in testlines.splitlines():
    for c in '+=-~*':
        line = line.replace(f'<many{c}>',c*80)
    testlines_expanded.append(line)

data = []
for i,line in enumerate(nc_createTextData('stdlib::Polyethylene_CH2.ncmat')):
    data.append( line )
    if i == 5:
        data.append( longstuff )
    if i == 8:
        data.append( longstuff.replace(' ','_') )
    if i == 10:
        data += testlines_expanded

#print( '\n'.join(data[:100]) )

pathlib.Path(longfn).write_text( '\n'.join(data) )

test_cfg( f'{longfn};vdoslux=1', material_name='CH2',
          check_teff=True, include_gif=True,
          ref_parsed={'tsl_H_in_CH2.endf':'0 0 1 451 7 2 7 4 7 451',
                      'tsl_C_in_CH2.endf':'0 0 1 451 7 2 7 4 7 451'},
          endf_metadata=metadata, dump_file=True)
