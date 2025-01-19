#!/usr/bin/env python3

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2024 NCrystal developers                                   ##
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

def platform_so_ending():
    import platform
    return '.dylib' if platform.system() == 'Darwin' else '.so'

def pkg_libname( pkgname ):
    libname = 'libPKG__%s'%pkgname
    return libname + platform_so_ending()

def main():
    import os
    pl = os.environ['SBLD_LIB_DIR']+'/'+pkg_libname('NCTestPlugin')
    if not os.path.exists(pl):
        raise SystemExit(f'Not found: {pl}')
    os.environ['NCRYSTALDEV_DEBUG_PLUGIN']='1'
    os.environ['NCRYSTALDEV_PLUGIN_LIST']=pl
    os.environ['NCRYSTALDEV_PLUGIN_RUNTESTS']='1'
    os.environ['NCRYSTALDEV_REQUIRED_PLUGINS']='DummyPlugin'
    for k,v in os.environ.items():
        if k.startswith('NCRYSTAL'):
            print(k,v)
    import NCrystalDev as NC
    pls = [e for e in NC.browsePlugins() if e[0]=='DummyPlugin']
    if not pls:
        raise SystemExit('Could not load DummyPlugin as expected')
    assert len(pls) == 1
    print(f'Verified loading of DummyPlugin from {pl}.')
    print('All ok')

if __name__ == '__main__':
    main()
