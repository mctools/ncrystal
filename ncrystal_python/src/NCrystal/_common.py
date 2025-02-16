
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

"""Module with internal utilities used by several NCrystal modules"""

#Common print function for NCrystal modules (allowing one to capture and
#redirect the output of NCrystal):

def _builtin_print():
    import builtins
    return builtins.print

_stdprint = [_builtin_print()]

def print(*args,**kwargs):
    _stdprint[0](*args,**kwargs)

#TODO: Overlaps somewhat with the (more complete) _setMsgHandler from _msg.py:
def set_ncrystal_print_fct( fct ):
    oldfct = _stdprint[0]
    _stdprint[0] = fct or _builtin_print()
    return oldfct

def get_ncrystal_print_fct():
    return _stdprint[0]

_cache_ns = [None]
def get_namespace():
    if _cache_ns[0] is None:
        from ._locatelib import get_libpath_and_namespace as f
        _p,_ns = f()
        _cache_ns[0] = _ns
    return _cache_ns[0]

_cache_nsev = [None]
def has_namespaced_envvars():
    if _cache_nsev[0] is None:
        import pathlib
        _cache_nsev[0] = ( pathlib.Path(__file__).parent
                           .joinpath('_do_namespace_envvars.py').is_file() )
    return _cache_nsev[0]

def expand_envname( name ):
    assert not name.startswith('NCRYSTAL_')#common mistake
    if has_namespaced_envvars():
        return 'NCRYSTAL' + get_namespace().upper()+'_' + name
    else:
        return 'NCRYSTAL_' + name

def ncgetenv( name, defval = None ):
    import os
    return os.environ.get( expand_envname(name), defval )

def ncgetenv_bool( name ):
    v = ncgetenv( name )
    if v is None or v=='0':
        return False
    return True

def ncgetenv_int( name, defval ):
    assert isinstance( defval, int)
    v = ncgetenv( name )
    if v is None:
        return defval
    try:
        v=int(v)
    except ValueError:
        n = expand_envname( name )
        raise ValueError(f"ERROR: {n} environment variable must"
                         " be integral value if set")
    return v

def ncgetenv_int_nonneg( name, defval ):
    assert isinstance( defval, int) and defval >= 0
    v = ncgetenv_int( name, defval )
    if v < 0:
        n = expand_envname( name )
        raise ValueError(f"ERROR: {n} environment variable must"
                         " be non-negative value if set")
    return v

def ncsetenv( name, val ):
    import os
    assert not name.startswith('NCRYSTAL_')#common mistake
    varname = 'NCRYSTAL' + get_namespace().upper()+'_' + name
    if val is None:
        if varname in os.environ:
            del os.environ[varname]
    else:
        os.environ[varname] = str(val)

class modify_ncrystal_print_fct_ctxmgr:
    """Context manager for modifying the print function (for instance to block
    output temporarily)"""
    def __init__(self, fct ):
        if fct == 'block':
            def fct( *a, **kw ):
                pass
        self.__fct = fct

    def __enter__(self):
        self.__orig = get_ncrystal_print_fct()
        set_ncrystal_print_fct( self.__fct )

    def __exit__(self,*args,**kwargs):
        set_ncrystal_print_fct( self.__orig )
        self.__orig = None
        self.__fct = None


def warn(msg):
    """Emit NCrystalUserWarning via standard warnings.warn function"""
    from .exceptions import NCrystalUserWarning
    import warnings
    warnings.warn( NCrystalUserWarning(str(msg)), stacklevel = 2 )

class WarningSpy:
    """Context manager which spies on any warnings emitted via warnings
    module and returns list of messages and categories.

    It does so by temporarily intercepting the warnings.showwarning
    function.

    Usage example:
        with WarningSpy() as warnlist:
            fct()
        print( "Warnings emitted:", warnlist )
    """

    def __init__(self,block=False,blockfct = None):
        assert not ( bool(block) and blockfct )
        self.__block = bool(block)
        self.__blockfct = blockfct
        self.__state = None

    def __enter__(self):
        #Make sure we create our own warnings.catch_warnings context
        #manager. This not only has the advantage of automatically restoring the
        #global warnings states (including the showwarning and simplefilter
        #function) upon __exit__, but also means that global warning
        #suppressions etc. will be temporarily disabled.
        self.__state = {}
        s = self.__state
        import warnings
        s['warncntxmgr'] = warnings.catch_warnings()
        s['warncntxmgr'].__enter__()
        s['orig'] = warnings.showwarning
        warnings.showwarning = self.__spy
        warnings.simplefilter("always")
        s['l'] = []
        return s['l']

    def __exit__(self,*args,**kwargs):
        assert self.__state
        import warnings
        warnings.showwarning = self.__state['orig']#not strictly needed since next line will also restore it
        self.__state['warncntxmgr'].__exit__()
        self.__state = None

    @staticmethod
    def _fmtwarning( message, category ):
        catname = None
        if category is not None:
            catname = str(getattr(category,'__name__',category)).strip()
        elif isinstance(message,Warning):
            catname = str(message.__class__.__name__).strip()
        if not catname:
            catname='Warning'
        return ( ( str(message).strip() or '<no message>'), catname )

    def __spy( self, message, category, *args, **kwargs ):
        #only store string objects, to decouple lifetimes of more complex objects:
        msg_str, cat_str = WarningSpy._fmtwarning( message, category )
        self.__state['l'].append( ( msg_str, cat_str ) )
        do_block = self.__block or ( self.__blockfct and self.__blockfct(msg_str, cat_str) )
        if not do_block:
            return self.__state['orig']( message, category, *args, **kwargs )

def find_fraction( x, tol = 1e-15, max_denom = 1000000 ):
    import fractions
    f = fractions.Fraction.from_float( x )
    f = f.limit_denominator(max_denom)
    if hasattr(f,'as_integer_ratio'):
        #py 3.8+
        a,b = f.as_integer_ratio()
    else:
        a,b = f.numerator,f.denominator
    if abs( a / b - x ) < tol:
        return a,b
    else:
        return None

def prettyFmtValue(x):
    """Recognises common fractions in truncated values like '0.33333333' or
    '0.4166667' and formats them as '1/3' and '5/12'. This not only saves space,
    but can also reinject accuracy in atom positions (ncrystal_verifyatompos can
    be used to check this, but it seems to be a clear improvement in all cases
    tested so far). Additionally also performs optimisations for storage
    efficiency, in particular omitting leading 0s with no information value
    (e.g. "0.2341" becomes ".2341").

    """
    assert 0.0 <= x <= 1.0
    if x==0.0:
        return '0'
    if x==1.0:
        return '1'
    if x==0.5:
        return '1/2'
    def stripleading0(s):
        return ( s[1:] if (len(s) > 2 and s.startswith('0.')) else s )
    xfmt = stripleading0('%.13g'%x)#.14g leads to irreproducibility issues in
    #our tests but there are sooo many numbers ending with 3333.... or
    #666667.... that we can safely "snap" these to their correct values:
    if xfmt[0]=='.' and xfmt.endswith('33333') and len(xfmt)<18:
        xfmt += '3'*(18-len(xfmt))
    if xfmt[0]=='.' and xfmt.endswith('66667') and len(xfmt)<18:
        xfmt = xfmt[:-1]+'6'*(18-len(xfmt))+'7'
    ff = find_fraction( x, tol = 1e-14, max_denom = 40 )
    if ff is not None and ff[1]!=1:
        assert ff[1] > 0
        v = ff[0] / ff[1]
        v_dbfmt = f'{ff[0]}/{ff[1]}'
        #check if credible that the numbers are identical (must be close AND the str
        #representations must be consistent with truncation):
        expected_at_prec = stripleading0((f'%.{len(xfmt)-1}g')%v)
        if abs(v-x)<1e-7 and expected_at_prec == xfmt:
            return v_dbfmt
        if abs(v-x)<1e-14:
            return v_dbfmt
    #No fraction found:
    if xfmt.isdigit() and float( xfmt ) != x:
        #abort, we don't want 0.99999999999 to print as '1' (which would then be
        #mapped to 0 in a unit cell):
        return stripleading0('%.19g'%x)
    return xfmt

def _split_trailing_digit( s ):
    if not s or not s[-1].isdigit():
        return s,None
    n=0
    while (n+1)<len(s) and s[-(n+1)].isdigit():
        n += 1
    if n:
        return s[0:-n], int( s[-n:] )
    else:
        return s, None

def check_elem_or_isotope_marker( s ):
    """If input is of form "Al", "O16", ..., return it, otherwise return
       None. Ignores excess whitespace in input. Will also map "H2"->"D" and
       "H3"->"T".
    """
    s=' '.join(str(s).split())
    if not s:
        return None
    elem_name, isotope_val = _split_trailing_digit( s )

    #special support for D,T:
    if not isotope_val and elem_name in ('D','T'):
        return elem_name
    if elem_name == 'H' and isotope_val in (2,3):
        return 'D' if isotope_val == 2 else 'T'

    from .atomdata import isElementName

    if isElementName( elem_name) and ( isotope_val is None or ( 0 <= isotope_val < 999 ) ):
        return elem_name if not isotope_val else '%s%i'%(elem_name,isotope_val)


def _hill_sort( chemform ):
    #takes chemform like [('Al',2'),('O',3)]  and sorts in order of Hill system notation

    #Remap H2/H3 and remove duplicates:
    remap = {'H2':'D','H3':'T'}
    if any(k in remap for k,v in chemform) or len(set(k for k,v in chemform))!=len(chemform):
        d={}
        for k,v in chemform:
            k = remap.get(k,k)
            d[k] = d.setdefault(k,0) + v
        return _hill_sort( list(d.items()) )

    has_carbon = any( en=='C' for en,c in chemform )
    if not has_carbon:
        return list( sorted( chemform ) )
    def hillsortkey( e ):
        #if not has carbon, then all in alphabetical order
        #first carbon, then H/D/T, then in alphabetical order
        if e[0]=='C':
            return ( -999999, e )
        is_hydrogen = ( e[0] in ('H','D','T') or ( e[0].startswith('H') and e[0][1:].isdigit() ) )
        return ( -1, e ) if is_hydrogen else (0, e )
    return list(sorted(chemform,key=hillsortkey))

def _gcd( *vals ):
    import math
    if len(vals)<2:
        return vals[0] if vals else None
    gcd = math.gcd(vals[0],vals[1])
    for e in vals[2:]:
        gcd = math.gcd( gcd, e )
    return gcd

def format_chemform( chemform, *, allow_rescaling = True ):
    #takes chemform like [('Al',2),('O',3)] and returns nicely formatted string
    #"Al2O3", with no duplicated element names and sorted according to the Hill
    #system of notation. Additionally, integral counts are divided by greatest
    #common divisor, and an attempt is made to scale up non-integral counts
    #where it makes sense (e.g. Al0.5O1.5 becomes AlO3).
    if len(chemform)==1:
        return str(chemform[0][0])
    cf = _hill_sort(chemform)
    def is_near_int(v):
        return abs(int(v)-v)<1e-15
    def get_non_ints(_l):
        return [ v for k,v in _l if not is_near_int(v) ]

    non_ints = get_non_ints( cf )
    def _find_frac(x):
        return find_fraction( x, max_denom = 100 )
    non_int_ffractions = [ _find_frac(x) for x in non_ints ]
    denoms = [ ff[1] for ff in non_int_ffractions if ff is not None ]
    cf_alt = None
    if allow_rescaling and denoms:
        import math
        factor = 1
        while denoms:
            _d = denoms.pop(0)
            _g = math.gcd( factor, _d )
            factor *= ( _d // _g )
        _cf = []
        for k,v in cf:
            if is_near_int(v):
                _cf.append( (k,v*factor ) )
            else:
                ff = _find_frac( v )
                if ff is None:
                    _cf.append( (k,v*factor ) )
                else:
                    assert factor % ff[1] == 0
                    _cf.append( (k,ff[0] * ( factor // ff[1] ) ) )
        if max( v for k,v in _cf ) < 100:
            cf_alt = _cf

    #find and apply gcd of all integers:
    def final_format( the_cf ):
        ll=[ int(v) for k,v in the_cf if is_near_int(v) ]
        gcd = _gcd( *ll ) if ( ll and allow_rescaling ) else 1
        def wrapiso( x ):
            #nb: these curly braces are not great for filenames...:
            return x if not x[-1].isdigit() else '{%s}'%x
        the_cf = [ (wrapiso(en),(count//gcd if is_near_int(count) else count/gcd)) for en,count in the_cf ]
        return ''.join( (en if count==1 else '%s%g'%(en,int(count) if count==int(count) else count)) for en,count in the_cf )

    f1 = final_format(cf)
    f2 = final_format(cf_alt) if cf_alt else f1
    return f1 if len(f1) < len(f2) else f2


def _classifySG(sgno):
    assert 1<=sgno<=230
    ll=[(195,'cubic'),(168,'hexagonal'),(143,'trigonal'),
        (75,'tetragonal'),(16,'orthorombic'),(3,'monoclinic'),(1,'triclinic')]
    for thr,nme in ll:
        if sgno>=thr:
            return nme
    assert False

#colors inspired by http://www.mulinblog.com/a-color-palette-optimized-for-data-visualization/
#supposedly from Stephen Few's book, "Show Me the Numbers":
_palette_Few = dict(red = "#F15854",
                    blue="#5DA5DA",
                    orange="#FAA43A",
                    green="#60BD68",
                    brown="#B2912F",
                    purple="#B276B2",
                    yellow="#DECF3F",
                    pink="#F17CB0",
                    gray="#4D4D4D")

def _grid_is_linspace( grid, tol = 1e-6 ):
    if len(grid)<=2:
        return len(grid)==2
    from ._numpy import _ensure_numpy, _np
    _ensure_numpy()
    g = _np.asarray(grid,dtype=float)
    bws = g[1:] - g[:-1]
    bmin, bmax = bws.min(), bws.max()
    if ( not ( bmax >= bmin > 0.0 ) ) or _np.isinf(bmin) or _np.isinf(bmax):
        return False
    return ( bmax - bmin ) < tol * ( bmax + bmin )

def extract_path( s ):
    """Returns pathlib.Path path object from argument if it is a path, otherwise
    None. Strings with no newlines will be assumed to be paths."""
    import pathlib
    if hasattr( s, '__fspath__' ):
        return pathlib.Path(s)
    if isinstance( s, str ) and '\n' not in s:
        res_try = _lookup_existing_file(s)
        if res_try:
            return res_try
        return pathlib.Path(s)

def download_url( url,
                  decode_as_utf8_str = True,
                  wrap_exception = True,
                  timeout = None,
                  quiet_network_fail = False ):
    """Download the provided url. The timeout value is in seconds if given.
    Returns None on network failure if quiet_network_fail is True.
    """
    import urllib.request
    import urllib.error
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout = timeout ) as response:
            data = response.read()
    except urllib.error.URLError as e:
        if quiet_network_fail:
            return None
        if wrap_exception:
            from .exceptions import NCException
            raise NCException(f'Error downloading url "{url}": {e}')
        else:
            raise e
    if decode_as_utf8_str:
        try:
            data = data.decode('utf8')
        except UnicodeDecodeError as e:
            if wrap_exception:
                from .exceptions import NCException
                raise NCException(f'Error decoding url to utf8 data "{url}": {e}')
            else:
                raise e
    return data

def _decodeflt(s):
    import numbers
    if s is None:
        return None
    if isinstance(s,numbers.Real):
        return float(s)
    try:
        x = float( s )
    except (TypeError,ValueError):
        return None
    return x

def _calc_md5hexdigest( str_or_bytes, / ):
    import hashlib
    if hasattr(str_or_bytes,'encode'):
        data = str_or_bytes.encode('utf8',errors='backslashreplace')
    else:
        data = str_or_bytes
    return hashlib.md5( data ).hexdigest()

def write_text( path, content ):
    """Like path.write_text(content) but forcing some global settings
    to ensure consistent NCrystal behaviour on any platform. Specifically
    encoding='utf8' and newline='\n'"""
    #The newline parameter was only added for pathlib.Path.write_text in python
    #3.10, so we do instead:
    import pathlib
    with pathlib.Path(path).open( 'wt',
                                  encoding = 'utf8',
                                  newline='\n' ) as fh:
        fh.write(content)

def _lookup_existing_file( path ):
    """Will wrap in pathlib.Path, except in simplebuild development mode where
    it will try $SBLD_DATA_DIR/+path when the path is a string with a single
    path separator. Returns None if unable to find an existing file at the
    location, otherwise an absolute path to the file.

    """
    import pathlib
    p = pathlib.Path( path ).expanduser()
    if p.exists():
        return p.absolute()
    if ( isinstance( path, str )
         and path.count('/')==1
         and ( pathlib.Path(__file__).parent
             .joinpath('_is_sblddevel.py').is_file() ) ):
        import os
        dd = os.environ.get('SBLD_DATA_DIR')
        if dd:
            p = pathlib.Path(dd).joinpath(path)
            if p.exists():
                return p.absolute()
