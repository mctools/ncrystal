
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

"""

   Module providing directories needed during testing.

"""

from pathlib import Path as _Path

_pymoddir = _Path(__file__).resolve().absolute().parent

def _find_data_dir():
    from .modeinfo import is_simplebuild_mode, is_ncrystalverify_mode
    if is_simplebuild_mode():
        import os
        ddir = _Path(os.environ['SBLD_DATA_DIR'])/'NCTestUtils'
    elif is_ncrystalverify_mode():
        ddir = _pymoddir.parent.parent.joinpath('data')
    else:
        ddir = _pymoddir.parent.parent.joinpath('data')
    assert ddir.is_dir()
    return ddir

test_data_dir = _find_data_dir()

def get_named_test_data_dir(name, for_updates = False ):
    #name of subdir of /tests/data, None means test_data_dir itself.  Set for
    #for_updates to True if running a test script in "--update" mode, where the
    #script will afterwards (over)write one or more files in it.
    tdir = test_data_dir
    if for_updates:
        tdir = _test_data_dir_for_updates()
        if not tdir:
            raise SystemExit('Can not update test data unless'
                             ' in NCrystal git repo')
    if for_updates:
        tdir = tdir.joinpath(name) if name else tdir
        tdir.mkdir( exist_ok = True )
        return tdir
    else:
        if not name:
            return test_data_dir
        from .modeinfo import is_simplebuild_mode, is_ncrystalverify_mode
        if is_simplebuild_mode():
            import os
            ddir = _Path(os.environ['SBLD_DATA_DIR'])/f'NCTestData_{name}'
        elif is_ncrystalverify_mode():
            ddir = _pymoddir.parent.parent.joinpath('data',name)
        else:
            ddir = _pymoddir.parent.parent.joinpath('data',name)
    if not ddir.is_dir():
        #We did not run for_updates, but the error message will depend on
        #whether or not we could have done so:
        tdir = _test_data_dir_for_updates()
        if tdir:
            raise SystemExit(f"ERROR: test directory ({tdir.joinpath(name)})"
                             " does not  exist.\nIf this is the first time"
                             " running the test, you can rerun with --update"
                             " to create the directory.")
        else:
            raise SystemExit(f"Error: test directory does not exist: {ddir}")
    return ddir

def _test_data_dir_for_updates():
    #If in NCrystal dev repo, return real tests/data dir in which updated files
    #can be placed, otherwise return None
    reporoot = _pymoddir.parent.parent.parent
    if not reporoot.joinpath('ncrystal_metapkg','pyproject.toml').is_file():
        return False
    if not reporoot.joinpath('ncrystal_core','pyproject.toml').is_file():
        return False
    if not reporoot.joinpath('.git','config').is_file():
        return False
    tdir = reporoot.joinpath('tests','data')
    if not tdir.is_dir():
        return False
    return tdir
