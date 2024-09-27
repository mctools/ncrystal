
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

"""
   Various utilities to be used when implementing tests.
"""

from NCrystal._testimpl import _work_in_tmpdir as work_in_tmpdir # noqa F401

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
            raise SystemExit(f'Emitted {exc_type.__name__} instead of the'
                             f' required {self.__et.__name__}!')
        val = exc_value.args
        if isinstance(val,tuple) and len(val)==1:
            val=val[0]

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
        from NCrystal._common import _calc_md5hexdigest
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
        from NCrystal._common import print as ncprint
        return ncprint
    elif printfct is None:
        from builtins import print as biprint
        return biprint
    else:
        return printfct
