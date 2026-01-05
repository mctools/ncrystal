
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


def load():
    from .srciter import all_files_iter
    from .dirs import testroot
    py = set( all_files_iter( '*.py', root = testroot.joinpath('scripts') ) )
    scripts = {}
    for f in py:
        bn = f.stem
        log = f.parent.joinpath(f'{f.stem}.log')
        scripts[bn] = dict( pyfile = f,
                            logfile = log if log.exists() else None )

    pypath = testroot.joinpath('pypath')
    pymods = set( all_files_iter( '*.py', root = pypath ) )

    testmods = {}
    for f in pymods:
        if f.parent.parent.samefile( pypath ):
            modname = f.parent.name
            if modname not in testmods:
                testmods[modname]=set()
            testmods[modname].add( f )

    ddir=testroot.joinpath('data')
    datafiles = set()
    for f in all_files_iter( root = ddir ):
        if '~' in f.name or '#' in f.name:
            continue
        if f.parent == ddir:
            subdir = ''
        else:
            assert f.parent.parent == ddir
            subdir = f.parent.name
        datafiles.add( ( subdir, f ) )

    return dict( scripts = scripts,
                 testmods = testmods,
                 datafiles = datafiles )
