
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

def main():
    modnameprefix = 'ncrystal_plugin_'
    import pkgutil
    names = [ name
              for finder, name, ispkg
              in pkgutil.iter_modules()
              if name.startswith(modnameprefix) ]
    if not names:
        return

    import importlib
    from pathlib import Path

    plugins = set()
    datadirs = set()
    for pymodname in names:
        mod = importlib.import_module(pymodname)
        moddir = Path(mod.__file__).parent
        for f in moddir.joinpath('plugins').glob('*NCPlugin*.*'):
            plugins.add( str(f.resolve().absolute()) )
        for datadir in moddir.glob('data*'):
            if datadir.is_dir() and any( True for p in datadir.iterdir() ):
                datadirs.add( (pymodname, str(datadir.resolve().absolute())) )

    entries = []
    for p in sorted(plugins):
        entries.append( p )
    nmnp = len(modnameprefix)
    for n,d in sorted(datadirs):
        entries.append( ':DATA:%s:%s'%(n[nmnp:],d) )
    print( ';\n'.join(entries) )

if __name__ == '__main__':
    main()
