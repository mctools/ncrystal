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
from NCTestUtils.loadlib import Lib
import NCrystalDev.exceptions as nc_exceptions
import pathlib
#from NCTestUtils.common import ( work_in_tmpdir,
#                                 #explicit_unicode_str,
#                                )

lib = Lib('ncmat')
lib.dump()
assert hasattr(lib,'nctest_tryParseNCMAT')

def tryParseNCMATData( data ):
    res = lib.nctest_tryParseNCMAT( data )
    if not res:
        return
    parts = res.split('@@@',1)
    if len(parts)==1:
        ename, emsg = parts[0], '<unknown reason>'
    else:
        ename, emsg = parts
    if not ename.startswith('NC'):
        raise RuntimeError(emsg)
    assert hasattr( nc_exceptions, ename )
    raise getattr( nc_exceptions, ename )( emsg )

def tryParseNCMATFromPath( path ):
    return tryParseNCMATData( pathlib.Path(path).read_text() )

def main( verbose ):
    testall(verbose = verbose)

def testall( verbose ):
    def print_sep(*args):
        print((('==== '+' '.join(args)+' ') if args else '').ljust(120,'='))

    from NCTestUtils.ncmatparser_testdata import testdata

    for i,data in enumerate(testdata):
        expect_badinput, data = ( (False,data[0])
                                  if isinstance(data,tuple)
                                  else (True,data) )
        print(f'\n\n==> Test {i}')
        #print_sep()
        def show_data():
            print(">>>TEST DATA Begin:")
            print(data)
            print(">>>TEST DATA End:")
        if verbose:
            show_data()

        if expect_badinput:
            try:
                tryParseNCMATData( data )
            except ( nc_exceptions.NCBadInput,
                     nc_exceptions.NCDataLoadError ) as e:
                print(f'==> Correct! {e.__class__.__name__}: {e}')
            except:
                print ('==> Ended in wrong type of exception!')
                show_data()
                raise
            else:
                print ('==> Parsed when it should not have!')
                show_data()
                raise SystemExit('Did not end in exception as it should!')
        else:
            try:
                tryParseNCMATData( data )
            except:
                print ('==> Ended with exception but it should not!')
                show_data()
                raise
            else:
                print("==> Ok.")

        #print_sep('Test %i done'%i)
    print()
    #doall = '--all' in sys.argv[1:]
    files = []
    #FIXMEfiles = (sorted(pathlib.Path(os.environ['SBLD_DATA_DIR']).glob('NCData/*.ncmat'))
    #FIXME         +sorted(pathlib.Path(os.environ['SBLD_DATA_DIR']).glob('NCDataRefNC1/*.ncmat'))
    #FIXME         +sorted(pathlib.Path(os.environ['SBLD_DATA_DIR']).glob('NCDataRefNC2d5/*.ncmat'))
    #FIXME         )

    #bad = False
    for p in files:
        tryParseNCMATFromPath(p)
#
#        #if p.name in {'CF2_sg14.ncmat','Al.ncmat'}:
#        #    #some actual bad file
#        #    continue
#        #if '@VDOS' in p.read_text() or '@SCATTERINGKERNEL' in p.read_text():
#        #    continue
#        sf = "%s/%s"%(p.parent.name,p.name) if not doall else str(p)
#        print("Now loading %s"%sf,end='')
#        ec=subprocess.run(['sb_nctests_parsencmat',str(p)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT).returncode
#        if ec!=0:
#            bad=True
#            print("  <========== FAILED!!!!!!")
#        else:
#            print()
#
#    if not doall:
#        print()
#        print ('Checked all NCData/*.ncmat files (supply --all to test ALL .ncmat files instead)!!\n')
#
    print()
    print('\nAll tests completed succesfully!')

if __name__ == '__main__':
    import sys
    main(verbose = ('-v' in sys.argv[1:]))
