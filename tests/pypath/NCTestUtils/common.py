
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
        val = exc_value.args
        if exc_type != self.__et:
            print('Message was: "%s"'%val)
            raise SystemExit(f'Emitted {exc_type.__name__}({exc_value}) instead of the'
                             f' required {self.__et.__name__}!')
        if isinstance(val,tuple) and len(val)==1:
            val=val[0]
        elif len(val)>1:
            val = str(exc_value) #seems to work for argparse.ArgumentError

        if self.__ev is not None and val != self.__ev:
            raise SystemExit(f'Expected {exc_type.__name__} did not have'
                             ' the correct value (expected'
                             f' {self.__ev!r} got {val!r}).!')

        pf = _resolve_print_fct(self.__printfct)
        pf(f"Caught expected {exc_type.__name__}({val!r})")
        return True

@_contextlib.contextmanager
def change_random_seed( seed ):
    """Context manager for temporarily changing random seed."""
    import random
    state = random.getstate()
    random.seed(seed)
    try:
        yield
    finally:
        random.setstate(state)

def fmt_args_as_str( *args, **kwargs ):
    return ', '.join( [ repr(a) for a in args ]
                      + [ f'{k}={v!r}' for k,v in sorted(kwargs.items()) ] )

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
    import NCrystalDev as NC
    import NCrystalDev._common as nc_common
    orig = nc_common.get_ncrystal_print_fct()
    if filtermap is None:
        filtermap = ( 'NCrystal v%s'%NC.__version__,
                      'NCrystal v<current>' )
    def version_filter( s ):
        return s.replace(*filtermap) if isinstance(s,str) else s
    def newprint( *a, **kwargs ):
        orig( *( version_filter(e) for e in a ), **kwargs)
    nc_common.set_ncrystal_print_fct(newprint)

def reldiff( x, y ):
    import math
    if math.isinf(x):
        return ( float('inf') if ( not math.isinf(y) or
               ( ( x>0 ) != ( y>0 ) ) ) else 0.0 )
    return abs(x-y)/(max(1e-300,abs(x)+abs(y)))

def require_flteq( x, y, tol = 1e-13 ):
    def okfct( a, b ):
        return bool( reldiff( a, b ) < tol )
    if hasattr( x, '__len__' ):
        if ( not len(x) == len(y) or
             any( ( not okfct(a,b) ) for a,b in zip(x,y) )):
            raise RuntimeError('numpy flteq failed for arrays '
                              f'x={x} and y={y}!')
    elif not okfct(x,y):
        raise RuntimeError(f'require_flteq( x={x}, y={y} ) failed!')

def interp1d(x, y):
    """Returns a function which interpolates linearly between {xi,yi}
    points. Extrapolation outside the range of x-values simply yields zero.
    """
    import numpy as np
    x, y = np.asarray(x, float), np.asarray(y, float)
    assert x.ndim == y.ndim == 1
    assert len(x) == len(y) > 0
    assert np.all(np.diff(x) > 0)
    def f(xq):
        xq = np.asarray(xq)
        return np.where( (xq < x[0]) | (xq > x[-1]), 0.0, np.interp(xq, x, y) )
    return f

def interp1d_loglin(x, y):
    """interpolate linearly in log(y), falling back to linear in y when at least
    one y-value is zero."""
    import numpy as np
    x, y = np.asarray(x, float), np.asarray(y, float)
    assert x.ndim == y.ndim == 1
    assert len(x) == len(y) > 1
    assert np.all(np.diff(x) > 0) and np.all(y >= 0)

    def f(q):
        q = np.asarray(q, float)
        inside = (q >= x[0]) & (q <= x[-1])
        q = np.clip(q, x[0], x[-1])

        i = np.clip(np.searchsorted(x, q, side="right") - 1, 0, len(x) - 2)
        t = (q - x[i]) / (x[i + 1] - x[i])
        a, b = y[i], y[i + 1]

        lin = a + t * (b - a)
        pos = (a > 0) & (b > 0)
        loglin = np.exp((1 - t) * np.log(np.where(pos, a, 1)) +
                        t * np.log(np.where(pos, b, 1)))

        return np.where(inside, np.where(pos, loglin, lin), 0.0)

    return f

def powspace(start, stop, num, p):
    import numpy as np
    assert num>=2 and num==int(num)
    if num == 2:
        return np.array([start,stop], dtype=float)
    res = start + (stop - start) * np.linspace(0.0, 1.0, num) ** p
    res[0] = start
    res[-1] = stop
    return res

def calc_reldiff( a,b ):
    import numpy as np
    if a.size != b.size:
        return np.inf
    return np.abs(a - b) / (0.5*(np.abs(b)+np.abs(a))+np.finfo(float).eps)

def thicken_grid(x, n):
    """Thickens grid of x values by inserting n extra points between each pair
    of existing points. New points are spaced out linearly.
    """
    import numpy as np
    x = np.asarray(x, float)
    t = np.arange(n + 1) / (n + 1)
    return np.r_[((x[1:] - x[:-1])[:, None] * t + x[:-1, None]).ravel(), x[-1]]
