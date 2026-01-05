
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

from contextlib import contextmanager as _ctxmgr

def get_nprocs( nice_factor = 0.9 ):
    import os
    if hasattr(os,'sched_getaffinity'):
        n = len(os.sched_getaffinity(0))
    else:
        import multiprocessing
        n = multiprocessing.cpu_count()
    n = min(1024,max(1,n))
    if n >= 4:
        #Be nice, leave a tiny bit for other tasks on the machine:
        n = round( n * nice_factor )
    return n

@_ctxmgr
def change_dir( path ):
    """Context manager for working in a directory (automatically
    created if doesn't exist) and then switching back"""
    import pathlib
    import os

    the_cwd = os.getcwd()
    p = pathlib.Path(path)
    p.mkdir( parents = True, exist_ok = True )
    try:
        os.chdir( p )
        yield
    finally:
        os.chdir( the_cwd )
    return

@_ctxmgr
def work_in_tmpdir():
    """Context manager for working in a temporary directory (automatically
    created+cleaned) and then switching back"""
    import os
    import tempfile
    the_cwd = os.getcwd()

    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            os.chdir(tmpdir)
            yield
        finally:
            os.chdir(the_cwd)#Important to leave tmpdir *before* deletion, to
                             #avoid PermissionError on Windows.

def path_is_relative_to( p, pother ):
    #Path.is_relative_to(..) was introduced in Python 3.9, this function lets us
    #support python 3.8.
    import pathlib
    assert isinstance( p, pathlib.Path )
    if hasattr(p,'is_relative_to'):
        #Python 3.9+:
        return p.is_relative_to(pother)
    else:
        #Python 3.8:
        try:
            p.relative_to(pother)
            return True
        except ValueError:
            return False
