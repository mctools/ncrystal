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

import NCTestUtils.enable_fpe # noqa F401
from NCrystalDev._locatelib import _search_env_overrides
import os
import pathlib
import tempfile

def main():
    orig = dict( (k,os.environ.get(k)) for k in
                 ('NCRYSTAL_LIB','NCRYSTAL_LIB_NAMESPACE_PROTECTION') )
    os.environ.pop('NCRYSTAL_LIB_NAMESPACE_PROTECTION',None)
    try:
        with tempfile.TemporaryDirectory() as td:
            for sub in ['plain','.venv','a.b/c.d']:
                d = pathlib.Path(td) / sub
                d.mkdir(parents=True)
                for fn in ['libNCrystal.so','libNCrystal-foo.so',
                           'libNCrystal-foo.so.4.4.7','libNCrystal-foo.dylib',
                           'NCrystal-foo.dll','NCrystal.dll']:
                    f = d / fn
                    f.touch()
                    os.environ['NCRYSTAL_LIB'] = str(f)
                    lib, ns, version = _search_env_overrides()
                    assert pathlib.Path(lib) == f
                    print(f'{sub+"/"+fn:>35} -> namespace {ns!r}')
    finally:
        for k,v in orig.items():
            if v is None:
                os.environ.pop(k,None)
            else:
                os.environ[k] = v

if __name__ == '__main__':
    main()
