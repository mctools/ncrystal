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

#test file with same element in different roles (here differing by MSD value)
import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
import sys

o = NC.directMultiCreate( """NCMAT v5
#Non-crystalline solid, one element as free gas
@DENSITY
  1 g_per_cm3
@DYNINFO
  element H
  fraction 2/3
  type vdosdebye
  debye_temp 300
@DYNINFO
  element C
  fraction 1/3
  type freegas
""")
o.info.dump()

o = NC.directMultiCreate( """NCMAT v4
#Non-crystalline solid, one element as free gas
@DENSITY
  1 g_per_cm3
@DYNINFO
  element H
  fraction 2/3
  type freegas
@DYNINFO
  element C
  fraction 1/3
  type freegas
""")
o.info.dump()


try:
    o = NC.directMultiCreate( """NCMAT v5
#crystal with missing MSD for some elements -> should fail.
@CELL
  cubic 4.0
@ATOMPOSITIONS
  Al 0 0 0
  C .5 .5 .5
@DYNINFO
  element Al
  fraction 1/2
  type vdosdebye
  debye_temp 300
@DYNINFO
  element C
  fraction 1/2
  type freegas
""")
except NC.NCBadInput:
    pass
else:
    raise SystemExit('ERROR: did not fail as required')


#help(o.scatter)
##o.absorption
##o.info
##
##NC.createInfo('Al_sg225_multi.ncmat;dcutoff=0.6').dump()


def teststrip(s):
    from NCrystalDev._ncmatimpl import _stripCommentsAndPackNCMATData
    sys.stdout.write( '-'*40+' orig ' + '-'*40 +'\n' )
    sys.stdout.write( s.rawData if hasattr(s,'rawData') else s )
    sys.stdout.write( '-'*40+' packed ' + '-'*40 +'\n' )
    sys.stdout.write( _stripCommentsAndPackNCMATData( s ) )

teststrip(NC.createTextData('stdlib::Al_sg225.ncmat'))
teststrip(NC.createTextData('C_sg194_pyrolytic_graphite.ncmat'))
teststrip("""NCMAT v5
# bla bla NCRYSTALMATCFG[temp=300
# bla bla
@CELL
 cubic 200
""")
teststrip("""N C M A T v5
# bla bla NCRYSTALMATCFG[temp=300
# bla bla
@CELL
 cubic 200
""")


def testcomments(s):
    import NCrystalDev._ncmatimpl as n
    f = n._extractInitialHeaderCommentsFromNCMATData
    print("Comments raw:")
    for c in f(s,dedent=False):
        print('  >>%s'%repr(c))
    print("Comments dedent:")
    for c in f(s,dedent=True):
        print('  >>%s'%repr(c))

testcomments("""NCMAT v5
#   bla bla NCRYSTALMATCFG[temp=300
#       foo
#   bla bla
@CELL
 cubic 200
""")
testcomments("""NCMAT v5
@CELL
 cubic 200
""")
