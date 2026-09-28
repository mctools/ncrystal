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

# Test that invalid values of environment variables result in exceptions
# which can be caught. This includes variables which used to be processed
# already during initialisation of global variables, where invalid values
# made it impossible to even load the NCrystal library.

import subprocess
import sys

import NCTestUtils.enable_fpe  # noqa F401
from NCrystalDev._common import expand_envname
from NCTestUtils.env import ncsetenv


def test( envvar, code ):
    #Must run in fresh processes, since variables are only read once:
    full = ( 'import NCrystalDev as NC\n'
             'print("NCrystal loaded OK")\n'
             'try:\n'
             + ''.join( f'    {line}\n' for line in code.splitlines() ) +
             '    print("Did not fail")\n'
             'except NC.NCBadInput as e:\n'
             '    print("NCBadInput: %s"%e)\n' )
    ncsetenv(envvar,'abc')
    try:
        rv = subprocess.run( [ sys.executable, '-c', full ],
                             capture_output = True, text = True,
                             check = False )
    finally:
        ncsetenv(envvar,None)
    out = rv.stdout.strip().splitlines()
    print(f'-- Invalid value of {envvar}:')
    for line in out:
        #Actual name of variable depends on the NCrystal namespace:
        print('   ' + line.replace( expand_envname(envvar), f'<{envvar}>' ) )
    assert rv.returncode == 0
    assert len(out) == 2
    assert out[0] == 'NCrystal loaded OK'
    assert out[1].startswith('NCBadInput: Invalid value of environment')
    assert envvar in out[1]

def main():
    loadal = 'NC.createScatter("stdlib::Al_sg225.ncmat;comp=inelas")'
    test( 'DEBUG_PHONON', loadal )
    test( 'DEBUG_FACTORY', loadal )
    test( 'DEBUGINFO', loadal )
    test( 'NCMAT_NOWARNFORCUSTOM',
          'NC.registerInMemoryFileData("custom.ncmat","NCMAT v3\\n'
          '@CUSTOM_BLA\\n  hello\\n@DENSITY\\n  1 g_per_cm3\\n'
          '@DYNINFO\\n  element H\\n  fraction 1\\n'
          '  type freegas\\n")\n'
          'NC.createInfo("custom.ncmat")' )

if __name__ == '__main__':
    main()
