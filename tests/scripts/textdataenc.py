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

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC

# Test for line encodings (even just printing
# getTextData('stdlib::somefile.ncmat').rawData gave problems on Windows before,
# due to encodings or line endings:

def query_print( key, load = True ):
    lines=[]
    print()
    print( "===============================================================" )
    print( "===============================================================" )
    print( "===============================================================" )
    print( f"=== NC.createTextData( {repr(key)} ===" )
    td = NC.createTextData(key)
    print( "=== Printing line by line ===" )
    for line in td:
        lines.append(line)
        print(line)
    print( "=============================" )
    assert '\r' not in td.rawData
    print( "=== Printing .rawData ===" )
    print( td.rawData, end='' )
    print( "=============================" )
    print( "=== Loading ===" )
    ok = False
    try:
        NC.load(f'{key};dcutoff=1.0;vdoslux=0')
        ok = True
    except NC.NCBadInput as e:
        if not load:
            print(f" ... got expected error: {e}")
        else:
            raise
    if ok and not load:
        raise SystemExit('Load did not fail as expected')
    if ok:
        print(" ... success.")
    print( "=============================" )
    return lines,td.rawData

query_print('stdlib::Ti_sg194.ncmat')

dummy_data = [
    b'NCMAT v7',
    b'@DENSITY',
    b'  1 g_per_cm3',
    b'@DYNINFO',
    b'  element  C',
    b'  fraction 1',
    b'  type     vdosdebye',
    b'  debye_temp 300',
]

dummy_data = b'<ENDOFLINE>'.join( dummy_data+[b''] )

dummy_data_decoded = dummy_data.decode('ascii').replace('<ENDOFLINE>','\n')
NC.registerInMemoryFileData( 'unix.ncmat',
                             dummy_data.replace(b'<ENDOFLINE>', b'\x0a' ) )
NC.registerInMemoryFileData( 'dos.ncmat',
                             dummy_data.replace(b'<ENDOFLINE>', b'\x0d\x0a' ) )
NC.registerInMemoryFileData( 'mixed.ncmat',
                             dummy_data
                             .replace( b'@DENSITY<ENDOFLINE>',
                                       b'@DENSITY\x0d\x0a')
                             .replace( b'<ENDOFLINE>',
                                       b'\x0a' )
                            )
#weird AKA "Ancient pre-OSX mac line endings are single \r". These can NOT be
#loaded, as they are disallowed by the NCMAT spec.:
#
# "Line endings can be Unix (LF) or DOS (CRLF), but not the ancient Macintosh
#  endings (CR)."

NC.registerInMemoryFileData( 'weird.ncmat',
                             dummy_data
                             .replace( b'@DENSITY<ENDOFLINE>',
                                       b'@DENSITY\x0d')
                             .replace( b'<ENDOFLINE>',
                                       b'\x0a' )
                            )
for enc in ('unix','dos','mixed','weird'):
    do_load = ( enc != 'weird' )
    lines,rawdata = query_print(f'virtual::{enc}.ncmat',
                                load = do_load )
    assert lines == dummy_data_decoded.splitlines(), f'Unexpected {enc}.ncmat lines'
    assert rawdata == dummy_data_decoded, f'Unexpected {enc}.ncmat .rawData'
