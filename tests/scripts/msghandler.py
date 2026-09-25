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

# Test the (internal) _msg._setMsgHandler, including resetting to the C++
# default handler by passing None.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
from NCrystalDev._msg import _setMsgHandler, _default_pymsghandler
import sys

#Debye temperature + VDOS triggers a deterministic NCRYSTAL_WARN on load:
data = """NCMAT v7
@CELL
 cubic 4.04958
@SPACEGROUP
 225
@ATOMPOSITIONS
 Al 0 1/2 1/2
 Al 0 0 0
 Al 1/2 1/2 0
 Al 1/2 0 1/2
@DEBYETEMPERATURE
 Al 400
@DYNINFO
 element Al
 fraction 1
 type vdos
 vdos_egrid 0.01 0.04
 vdos_density 1 4 9 16 25 20 5
"""

nload = [0]
def load():
    #Unique comment to avoid any caching:
    nload[0] += 1
    NC.directLoad( data + f'#{nload[0]}\n',
                   doScatter = False, doAbsorption = False )

unraisable = []
sys.unraisablehook = lambda u : unraisable.append(u)

def flush():
    sys.stdout.flush()
    sys.stderr.flush()

print('--- Custom handler:')
got = []
_setMsgHandler( lambda msg, msgtype : got.append( (msgtype,msg) ) )
load()
print('Custom handler got:',got)
assert len(got)==1

print('--- C++ default handler (set via None):')
flush()
_setMsgHandler( None )
load()
flush()
assert not unraisable, unraisable

print('--- Python default handler:')
_setMsgHandler( _default_pymsghandler )
load()
assert len(got)==1
assert not unraisable
