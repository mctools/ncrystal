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

# Test inference of namespace from NCRYSTAL_LIB (only file name should
# matter, not dots in parent dirs).

#Staged progress prints below (flushed): this test was observed to hang
#with zero output on the GitHub windows-2025 runners, so make the next
#such hang reveal how far it got (the print-interleaved imports are
#deliberate, hence the noqa markers):
print('locatelib: begin imports',flush=True)
import NCTestUtils.enable_fpe # noqa F401
print('locatelib: enable_fpe imported',flush=True)
from NCrystalDev._locatelib import _search_env_overrides # noqa: E402, I001
print('locatelib: _locatelib imported',flush=True)
import os # noqa: E402, I001
import pathlib # noqa: E402
import tempfile # noqa: E402

def main():
    orig = { k: os.environ.get(k) for k in
             ('NCRYSTAL_LIB','NCRYSTAL_LIB_NAMESPACE_PROTECTION') }
    os.environ.pop('NCRYSTAL_LIB_NAMESPACE_PROTECTION',None)
    try:
        with tempfile.TemporaryDirectory() as td:
            print('locatelib: tempdir created',flush=True)
            for sub in ['plain','.venv','a.b/c.d']:
                d = pathlib.Path(td) / sub
                d.mkdir(parents=True)
                for fn in ['libNCrystal.so','libNCrystal-foo.so',
                           'libNCrystal-foo.so.4.4.7','libNCrystal-foo.dylib',
                           'NCrystal-foo.dll','NCrystal.dll']:
                    f = d / fn
                    f.touch()
                    os.environ['NCRYSTAL_LIB'] = str(f)
                    lib, ns, _version = _search_env_overrides()
                    assert pathlib.Path(lib) == f
                    print(f'{sub+"/"+fn:>35} -> namespace {ns!r}')
    finally:
        for k,v in orig.items():
            if v is None:
                os.environ.pop(k,None)
            else:
                os.environ[k] = v

def test_failed_load():
    #A failed library load must fail again on retry (not return False):
    import subprocess
    import sys
    code = '''
from NCrystalDev import _chooks
for i in range(2):
    try:
        _chooks._get_raw_cfcts()
        print(f'call {i}: no error')
    except OSError:
        print(f'call {i}: OSError')
'''
    with tempfile.TemporaryDirectory() as td:
        f = pathlib.Path(td) / 'libNCrystal.so'
        f.write_text('not a shared library')
        env = os.environ.copy()
        env['NCRYSTAL_LIB'] = str(f)
        env['NCRYSTAL_SLIMPYINIT'] = '1'
        rv = subprocess.run( [sys.executable,'-c',code], env = env,
                             capture_output = True, text = True,
                             check = True )
    print('Loading bogus library twice:')
    print(rv.stdout,end='')

if __name__ == '__main__':
    main()
    test_failed_load()
