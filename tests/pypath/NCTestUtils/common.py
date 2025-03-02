
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

"""
   Various utilities to be used when implementing tests.
"""

import contextlib as _contextlib
@_contextlib.contextmanager
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

_is_windows = [None]
def is_windows():
    if _is_windows[0] is None:
        import platform
        _is_windows[0] = platform.system().lower()=='windows'
    return _is_windows[0]

class ensure_error:

    """
    For usage in testing code where exceptions are expected, like:

    with ensure_error(NC.NCBadInput,'Some error message'):
        do_something_here()
    """
    def __init__(self, exc_type, exc_value = None, printfct = 'ncrystal' ):
        assert exc_type is not None
        self.__et = exc_type
        self.__ev = exc_value
        self.__printfct = printfct

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is None:
            raise SystemExit('Did not emit exception as required!')
        if exc_type != self.__et:
            raise SystemExit(f'Emitted {exc_type.__name__}({exc_value}) instead of the'
                             f' required {self.__et.__name__}!')
        val = exc_value.args
        if isinstance(val,tuple) and len(val)==1:
            val=val[0]
        elif len(val)>1:
            val = str(exc_value) #seems to work for argparse.ArgumentError

        if self.__ev is not None and val != self.__ev:
            raise SystemExit(f'Expected {exc_type.__name__} did not have'
                             ' the correct value (expected'
                             f' {repr(self.__ev)} got {repr(val)}).!')

        pf = _resolve_print_fct(self.__printfct)
        pf(f"Caught expected {exc_type.__name__}({repr(val)})")
        return True

def fmt_args_as_str( *args, **kwargs ):
    return ', '.join( [ repr(a) for a in args ]
                      + [ f'{k}={repr(v)}' for k,v in sorted(kwargs.items()) ] )

def print_text_file_with_snipping(content,
                                  nstart=30,
                                  nend=20,
                                  printfct = 'ncrystal',
                                  prefix=''):
    """Prints text files, but snips out the middle part of larger
    files. Printout includes a checksum of the snipped part."""
    nstart = max(3,nstart)
    nend = max(3,nend)
    print = _resolve_print_fct(printfct)
    lines=content.splitlines()
    if len(lines) < int((nstart+nend)*1.5+1):
        for line in lines:
            print(f'{prefix}{line}')
    else:
        for i in range(nstart):
            print(f'{prefix}{lines[i]}')
        from NCrystalDev._common import _calc_md5hexdigest
        md5 = _calc_md5hexdigest( '\n'.join(lines[nstart:-nend]) )
        def nleading_spaces( s ):
            return len(s)-len(s.lstrip(' '))
        nspaces = min(nleading_spaces(lines[nstart-1]),
                      nleading_spaces(lines[-nend]))
        spaces = ' '*nspaces
        print(f"{prefix}{spaces}<<<SNIPPED {len(lines)-nstart-nend} LINES,"
              f" MD5={md5}>>>")
        for i in range(nend):
            print(f'{prefix}{lines[-nend+i]}')

def _resolve_print_fct(printfct):
    if printfct == 'ncrystal':
        from NCrystalDev._common import print as ncprint
        return ncprint
    elif printfct is None:
        from builtins import print as biprint
        return biprint
    else:
        return printfct

def explicit_unicode_char(c):
    #32 is space, <32 are control chars, 127 is DEL.
    return c if 32<=ord(c)<=126 else r'\u{%s}'%(hex(ord(c))[2:])
def explicit_unicode_str(s):
    return ''.join( explicit_unicode_char(c) for c in s)

def fix_ncrystal_version_printouts( filtermap = None ):
    import NCrystalDev._common as nc_common
    import NCrystalDev as NC
    orig = nc_common.get_ncrystal_print_fct()
    if filtermap is None:
        filtermap = ( 'NCrystal v%s'%NC.__version__,
                      'NCrystal v<current>' )
    def version_filter( s ):
        return s.replace(*filtermap) if isinstance(s,str) else s
    def newprint( *a, **kwargs ):
        orig( *( version_filter(e) for e in a ), **kwargs)
    nc_common.set_ncrystal_print_fct(newprint)
