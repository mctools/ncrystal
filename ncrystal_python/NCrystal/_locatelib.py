
"""Internal module providing ctypes-based hooks into the compiled NCrystal shared library"""

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

import pathlib

_cache = [None]
def get_libpath_and_namespace():
    if _cache[0] is None:
        _cache[0] = _search()
    return _cache[0]

def _detect_monolithic_installation():
    #Determine if these python modules were installed via "pip install
    #<ncrystalrepo>" (monolithic package with bore core and python parts), or if
    #using a standard ncrystal-python package. Also detect conflicts in case a
    #user did both, and thus messed up their environment.

    d = pathlib.Path(__file__).parent
    is_std = d.joinpath('_is_std.py').is_file()
    is_mono = d.joinpath('_is_monolithic.py').is_file()
    assert is_mono or is_std
    if not ( is_mono and is_std ):
        return is_mono

    import textwrap
    conflict_msg = textwrap.dedent("""

    ERROR: Broken environment detected.

    The current environment has traces of both a normal NCrystal installation
    (where ncrystal-core and ncrystal-python modules are installed separately),
    and a monolithic installation (most likely from running "pip install"
    directly on the root of the ncrystal source tree). This latter method is
    primarily intended for developers, and should never be performed in an
    environment where NCrystal is already installed through other conventional
    methods.

    Most likely your environment is now broken. You can attempt to fix it by
    uninstalling all NCrystal packages (with the same package manager that did
    the original installation of that package), and then installing NCrystal
    again from one source only. Or you might simply wish to remove and recreate
    the environment from scratch.
    """)
    #We can not simply emit ImportError, since it might be taken by clients as
    #"ncrystal not installed". We can choose RuntimeError or SystemExit. To be
    #100% this triggers a warning, we pick SystemExit.
    raise SystemExit(conflict_msg)

def _search():
    #Try to determine the location of the NCrystal shared library. In order of
    #preference, we try:
    #
    #  1) The NCrystal_LIB environment variable (intended for specialised usage,
    #     like CTests).
    #     Note that an associated NCRYSTAL_LIB_NAMESPACE_PROTECTION variable must
    #     be used for any namespace (fixme: decode the namespace from the
    #     NCrystal_LIB file name?).
    #
    #  2) Attempt to get the location via a python module providing it (in case
    #     ncrystal-core was installed via PyPI for instance).
    #     FIXME: still local _nclibpath for now!
    #
    #  3) Invoke "ncrystal-config --show libpath namespace" for the information.
    #

    #always invoke _detect_monolithic_installation() since it also detects
    #broken installations.
    verbose = False#fixme
    if verbose:
        from ._common import print
        print('NCrystal._locatelib: Starting search for NCrystal shared library')
    is_monolithic = _detect_monolithic_installation()
    if verbose:
        print('NCrystal._locatelib: monolithic installation'
              f' = {"yes" if is_monolithic else "no"}')

    v = _search_env_overrides()
    if verbose and v:
        print('NCrystal._locatelib: Succesfully searched via method: env vars')
    if not v:
        v = _search_core_info_mod( is_monolithic )
        if verbose and v:
            print('NCrystal._locatelib: Succesfully searched via method: pymod')

    if not v:
        v = _search_nccfgapp()
        if verbose and v:
            print('NCrystal._locatelib: Succesfully searched'
                  ' via method: ncrystal-config')

    if not v:
        from .exceptions import NCFileNotFound
        raise NCFileNotFound('Could not locate the NCrystal shared library'
                             ' file. Have you installed the ncrystal-core'
                             ' package?')
    lib, namespace, version = v
    if verbose:
        print(f'NCrystal._locatelib: namespace = "{namespace}"')
        print(f'NCrystal._locatelib: lib = "{lib}"')

    lib = pathlib.Path(lib)
    if not lib.is_file():
        from .exceptions import NCFileNotFound
        raise NCFileNotFound('Problems locating the NCrystal shared library'
                             f' file. Tried "{lib}" unsuccesfully.')
    version_ok = True
    if version is not None:
        from . import __version__ as _nc_version
        version_ok = _nc_version.strip() == version.strip()
    if version_ok:
        return ( lib.absolute().resolve(), ( namespace or '' ) )

    #Found it, but there was a version mismatch!
    import textwrap
    msg = textwrap.dedent(f"""
    ERROR: Inconsistent environment detected.

    The version of the installed core ncrystal modules ({version}) and the
    ncrystal Python API ({_nc_version}) are not identical. This is not
    supported.
    """)
    #The environment is not necessarily completely broken, but the NCrystal
    #packages can not be used. Hence, we emit an ImportError here. However, we
    #also emit a warning, in case the client code is hiding ImportError's, on
    #the assumption it simply means NCrystal is not available.
    from ._common import warn
    warn(msg)
    raise ImportError(msg)

def _search_env_overrides():
    import os
    lib = os.environ.get('NCRYSTAL_LIB')
    ns = os.environ.get('NCRYSTAL_LIB_NAMESPACE_PROTECTION','')
    if lib:
        if not ns and '.' in lib and 'NCrystal-' in lib:
            #Try to infer the namespace from the library name (so it is enough
            #to set NCRYSTAL_LIB).
            ll = lib.split('.')[0]
            if 'NCrystal-' in ll:
                ll = lib.split('NCrystal-')[-1]
                if ll and 'NCrystal-' not in ll and '.' not in ll:
                    ns = ll
        lib = pathlib.Path(lib)
        if not lib.exists() or lib.is_dir():
            from .exceptions import NCFileNotFound
            raise NCFileNotFound('NCRYSTAL_LIB environment variable is set'
                                 f' ("{lib}") but does not point'
                                 ' to an actual file.')
        return lib, ns, None

def _search_core_info_mod( is_monolithic ):
    #Look for the standard module, installed by the ncrystal-core package.
    try:
        if is_monolithic:
            import _ncrystal_core_monolithic.info as mod
        else:
            import _ncrystal_core.info as mod
    except ModuleNotFoundError:
        mod = None
    if mod:
        return mod.libpath(), mod.namespace(), mod.version()

def _search_nccfgapp():
    #Try to query ncrystal-config script:
    import subprocess
    cmdname='ncrystal-config'
    res = subprocess.run([cmdname,'--show','shlibpath','namespace','version'],
                         capture_output=True)
    if res.returncode == 0:
        lines = res.stdout.decode('utf8').splitlines()
        if len(lines) == 3:
            p = pathlib.Path(lines[0])
            if p.is_file():
                return p, ( lines[1].strip() or None ), lines[2].strip()
