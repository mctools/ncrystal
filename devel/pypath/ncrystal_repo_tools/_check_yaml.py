
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

def main():
    import yaml
    from .srciter import all_files_iter
    from .dirs import reporoot

    for f in all_files_iter('yaml'):
        print("  Trying to load %s"%f.relative_to(reporoot))
        with f.open() as fh:
            yaml.safe_load(fh)

    #Extra check of plugin database file:
    from . import _climode_plugindb
    print(f"  Extra verification of {_climode_plugindb.dbfilerelpath}")
    pdb_data = _climode_plugindb.load_and_check_data()
    s = _climode_plugindb._produce_wiki( pdb_data )
    assert len(s.splitlines())>20
    print('all ok')

if __name__=='__main__':
    main()
