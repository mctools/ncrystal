
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
    import pkgutil
    names = [ name
              for finder, name, ispkg
              in pkgutil.iter_modules()
              if name.startswith('ncplugin_') ]
    if not names:
        return

    import importlib
    from pathlib import Path

    plugins = set()
    for pymodname in names:
        mod = importlib.import_module(pymodname)
        for f in Path(mod.__file__).parent.joinpath('plugins').glob('*NCPlugin*.*'):
            plugins.add( str(f.resolve().absolute()) )
    print(';'.join(sorted(plugins)))

if __name__ == '__main__':
    main()

#
#discovered_plugins = {
#    name: importlib.import_module(name)
# . ...
