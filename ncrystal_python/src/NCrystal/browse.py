
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

"""Utilities for browsing and searching the data files available to NCrystal
(e.g. files in the standard data library, in-memory files, or files in the
current directory), optionally based on the physics content of the materials
they describe.

Example:

  import NCrystal.browse as nb
  for e in nb.find( where = "'B' in elements and absxs > 100" ):
      print( e.fullkey, e.description )

The data is provided by the C++ layer (via the JSON queries
["util","browsedb",...]), which also loads materials in parallel.
"""

__all__ = [ 'DataEntry', 'PhysicsProps', 'browse', 'filter_entries', 'find',
            'list_all_factories', 'list_factories', 'physics_props_doc' ]

def list_factories():
    """Returns a dictionary with the names of all TextData factories and the
    number of browsable entries they provide."""
    return _q('browsedb')

def list_all_factories():
    """Returns a dictionary with lists of the names of all factories, by
    type ("textdata", "info", "scatter", and "absorption")."""
    return _q('browsefactories')

def physics_props_doc():
    """Returns list of (name,description) of the properties available on
    PhysicsProps objects (and in "where" expressions)."""
    return list( _propdocs )

class PhysicsProps:
    """Physics properties of a loaded material. They are available as
    attributes (e.g. props.absxs), or as a dictionary via .as_dict(). See
    physics_props_doc() for a list with descriptions."""

    def __init__( self, info ):
        """For internal usage only (info is the "info" dictionary from the
        C++ browsedb query)."""
        from ._common import format_chemform
        from .atomdata import elementZToName
        compos = tuple( ( Z, tuple( (A,fr) for A,fr in isotopes ) )
                        for Z, isotopes in info['composition'] )
        atomfracs = {}
        for Z, isotopes in compos:
            en = elementZToName(Z)
            for A, fr in isotopes:
                lbl = ( en if A == 0 else
                        { (1,2):'D', (1,3):'T' }.get( (Z,A), f'{en}{A}' ) )
                atomfracs[lbl] = atomfracs.get(lbl,0.0) + fr
        elements = frozenset( elementZToName(Z) for Z,_ in compos )
        self.__d = dict(
            elements = elements,
            atoms = frozenset( atomfracs ),
            nelements = len( elements ),
            formula = format_chemform( sorted( atomfracs.items() ) ),
            absxs = info['xsect_absorption'],
            scatxs = info['xsect_free'],
            density = info['density'],
            numdens = info['numberdensity'],
            temp = info['temperature'],
            state = info['stateofmatter'].lower(),
            crystal = info['crystalline'],
            sg = info['spacegroup'],
            natoms = info['natoms_unitcell'],
            dyninfo = frozenset( info['dyninfo_types'] ),
            nphases = info['nphases'],
        )
        self.__compos = compos

    def as_dict( self ):
        """The properties as a (new) dictionary."""
        return dict( self.__d )

    @property
    def composition( self ):
        """The flattened composition as ((Z,((A,fraction),...)),...), where A=0
        indicates natural elements."""
        return self.__compos

    def __getattr__( self, name ):
        d = self.__dict__.get('_PhysicsProps__d')
        if d is not None and name in d:
            return d[name]
        raise AttributeError(f'PhysicsProps has no attribute "{name}"')

    def __str__( self ):
        return 'PhysicsProps(%s)'%( ', '.join( f'{k}={_fmt_prop(v)}'
                                               for k,v in self.__d.items() ) )

    def __repr__( self ):
        return str(self)

class DataEntry:
    """An entry (e.g. a file) available to NCrystal. Objects are created by
    the browse() or find() functions."""

    def __init__( self, data ):
        """For internal usage only (data is a dictionary from the C++
        browsedb query)."""
        self.__d = data
        c = data.get('comments')
        self.__comments = tuple(c) if c is not None else None
        self.__props = ( PhysicsProps( data['info'] ) if 'info' in data
                         else None )

    @property
    def name( self ):
        """Name (e.g. file name) used to request the entry."""
        return self.__d['name']

    @property
    def fullkey( self ):
        """The string "<factory>::<name>", which can be used to explicitly
        request this entry, even if it is hidden."""
        return self.__d['fullkey']

    @property
    def factory( self ):
        """Name of the factory delivering the entry."""
        return self.__d['factory']

    @property
    def source( self ):
        """Description of source (e.g. a directory path)."""
        return self.__d['source']

    @property
    def priority( self ):
        """Priority of the entry (an integer, or the string
        "OnlyOnExplicitRequest")."""
        return self.__d['priority']

    @property
    def hidden( self ):
        """True if the name is shadowed by an entry with the same name from a
        higher priority source (use .fullkey to select it explicitly)."""
        return self.__d['hidden']

    @property
    def display_name( self ):
        """The name, or the full key if the entry can not be selected by its
        name alone (because it is hidden or needs an explicit request)."""
        needs_key = ( self.hidden
                      or self.priority == 'OnlyOnExplicitRequest' )
        return self.fullkey if needs_key else self.name

    @property
    def datatype( self ):
        """Data type (e.g. "ncmat"), or None if data could not be read."""
        return self.__d.get('datatype')

    @property
    def comments( self ):
        """Initial comment lines of NCMAT data (tuple of str, without the
        leading '#' and dedented), or None if not NCMAT data."""
        return self.__comments

    @property
    def description( self ):
        """Short description, based on the first paragraph of the NCMAT
        header comments (empty string if not available)."""
        return _short_descr( self.__comments )

    @property
    def props( self ):
        """PhysicsProps of the loaded material, or None if the material was not
        loaded (or could not be loaded, see .error)."""
        return self.__props

    @property
    def error( self ):
        """Error message if loading failed, otherwise None."""
        return self.__d.get('error')

    def matching_lines( self, search, *, regex = False ):
        """Lines of the header comments matching any of the search words
        (case-insensitive, or regular expressions if regex=True)."""
        res = _compile_words( search, regex )
        return [ ll.strip() for ll in ( self.__comments or [] )
                 if any( r.search(ll) for r in res ) ]

    def textdata( self ):
        """Returns the NCrystal.TextData object of the entry."""
        from .core import createTextData
        return createTextData( self.fullkey )

    def load( self, cfg_params = '' ):
        """Load the material (returns NCrystal.LoadedMaterial). Additional
        cfg parameters can be provided (e.g. ";temp=20K")."""
        from .core import load
        return load( self.fullkey + cfg_params )

    def __str__( self ):
        return f'DataEntry({self.fullkey})'

    def __repr__( self ):
        return str(self)

def browse( factory = None, *, load = False, nthreads = 'auto',
            quiet = True, progress = None ):
    """Returns a list of DataEntry objects for all available entries, or those
    of a given factory. Entries are sorted by priority (highest first),
    factory, source, and name.

    If load=True, all materials are loaded, and their physics properties are
    available via the .props attribute of the entries. This is done in
    parallel with nthreads threads (an integer, or "auto"), which are only
    used temporarily, and only if the user did not already configure
    NCrystal's factory threads (via enableFactoryThreads or the
    NCRYSTAL_FACTORY_THREADS env var, which are then respected). If
    quiet=True, messages from NCrystal during loading are suppressed. If
    progress is a function, it will be called as progress(ndone,ntotal)
    during loading.
    """
    counts = list_factories()
    if factory is not None and factory not in counts:
        from .exceptions import NCBadInput
        raise NCBadInput(f'Unknown TextData factory: "{factory}"')
    facts = [ f for f,n in counts.items()
              if n and ( factory is None or f == factory ) ]
    ntotal = sum( counts[f] for f in facts )
    extra = ( ( f'nthreads={nthreads}', ) if load else ('cheap',) )
    entries, ndone = [], 0
    from ._msg import _suppress_msgs_ctx
    msgctx = _suppress_msgs_ctx() if ( load and quiet ) else _nullctx()
    with msgctx:
        for f in facts:
            #Smaller chunks if progress reporting is needed:
            nchunks = ( 1 if progress is None
                        else max( 1, -(-counts[f]//_progress_chunk) ) )
            for i in range( nchunks ):
                if progress is not None:
                    progress( ndone, ntotal )
                data = _q( 'browsedb', f, str(i), str(nchunks), *extra )
                entries += [ DataEntry(d) for d in data ]
                ndone += len(data)
    if progress is not None:
        progress( ndone, ntotal )
    entries.sort( key = _sortkey )
    return entries

def filter_entries( entries, *, patterns = (), search = (), regex = False,
                    where = () ):
    """Select entries from a list. Entries must match at least one of the
    name patterns (if any), contain all the search words in their name or
    header comments, and fulfil all where conditions.

    Patterns are case-insensitive substrings of names, or glob patterns if
    they contain "*", "?", or "[". Patterns containing "::" are matched
    against the full key (e.g. "stdlib::Al*"). If regex=True, patterns and
    search words are instead (case-insensitive) Python regular expressions,
    where ^ and $ match at line boundaries.

    Where conditions are either functions taking a PhysicsProps object, or
    Python expressions (strings) using the property names (e.g.
    "'B' in elements and absxs > 100"). They require loaded entries, and
    entries which could not be loaded are never selected.
    """
    patterns, search, where = map( _aslist, ( patterns, search, where ) )
    pattern_res = _compile_words( patterns, regex ) if regex else None
    search_res = _compile_words( search, regex )
    where_fcts = [ ( w if callable(w) else _WhereExpr(w) ) for w in where ]
    res = []
    for e in entries:
        if patterns:
            prs = pattern_res or [ None ]*len(patterns)
            if not any( _name_matches( e, p, pr )
                        for p, pr in zip( patterns, prs ) ):
                continue
        if search_res:
            text = '\n'.join( [ e.name ] + list( e.comments or [] ) )
            if not all( r.search( text ) for r in search_res ):
                continue
        if where_fcts:
            if e.props is None:
                if e.error is None:
                    from .exceptions import NCBadInput
                    raise NCBadInput('Where conditions require entries'
                                     ' loaded with load=True')
                continue
            if not all( _eval_where( w, e.props ) for w in where_fcts ):
                continue
        res.append( e )
    return res

def find( *patterns, search = (), regex = False, where = (),
          factory = None, load = None, nthreads = 'auto', quiet = True,
          progress = None ):
    """Browse and filter in one go (see the browse and filter_entries
    functions for details). Materials are loaded if load=True, or if
    load=None (the default) and there are where conditions."""
    if load is None:
        load = bool( _aslist( where ) )
    entries = browse( factory, load = load, nthreads = nthreads,
                      quiet = quiet, progress = progress )
    return filter_entries( entries, patterns = patterns, search = search,
                           regex = regex, where = where )

###############################################################################
# Implementation details:

_progress_chunk = 16

_propdocs = [
    ('elements', 'set of element names, e.g. {"Al","O"}'),
    ('atoms', 'set of atom labels, including isotopes, e.g. {"D","O"}'),
    ('nelements', 'number of different elements'),
    ('formula', 'chemical formula, e.g. "Al2O3"'),
    ('absxs', 'absorption cross section per atom at 2200m/s [barn]'),
    ('scatxs', 'free scattering cross section per atom [barn]'),
    ('density', 'density [g/cm3]'),
    ('numdens', 'number density [atoms/Aa3]'),
    ('temp', 'temperature [K]'),
    ('state', 'state of matter: "solid", "liquid", "gas", or "unknown"'),
    ('crystal', 'True if crystalline (for multiphase: any phase)'),
    ('sg', 'space group number (None if not available)'),
    ('natoms', 'number of atoms in unit cell (None if not available)'),
    ('dyninfo', ('set of dynamic info types present, among "vdos",'
                 ' "vdosdebye", "scatknl" (a full scattering kernel),'
                 ' "freegas", and "sterile"')),
    ('nphases', 'number of phases (1 for single-phase materials)'),
]

_where_funcs = dict( len = len, min = min, max = max, any = any, all = all,
                     abs = abs, round = round, set = set, sorted = sorted )

def _q( *args ):
    from .misc import evaluate_query
    return evaluate_query( ['util'] + list(args) )

def _aslist( x ):
    return [ x ] if isinstance( x, str ) or callable( x ) else list( x )

class _nullctx:
    def __enter__( self ):
        pass
    def __exit__( self, *a ):
        pass

def _sortkey( e ):
    p = e.priority
    rank = p if isinstance( p, int ) else ( -1 if p == 'OnlyOnExplicitRequest'
                                            else -2 )
    return ( -rank, e.factory, e.source, e.name )

def _is_glob( pattern ):
    return any( c in pattern for c in '*?[' )

def _compile_words( words, regex ):
    #Case-insensitive regexes (MULTILINE: ^ and $ match at line boundaries):
    import re
    res = []
    for w in _aslist( words ):
        try:
            res.append( re.compile( w if regex else re.escape(w),
                                    re.IGNORECASE | re.MULTILINE ) )
        except re.error as e:
            from .exceptions import NCBadInput
            raise NCBadInput(f'Invalid regular expression "{w}": {e}') from e
    return res

def _name_matches( entry, pattern, pattern_re = None ):
    s = entry.fullkey if '::' in pattern else entry.name
    if pattern_re is not None:
        return bool( pattern_re.search( s ) )
    import fnmatch
    p, s = pattern.lower(), s.lower()
    return fnmatch.fnmatchcase( s, p ) if _is_glob( p ) else ( p in s )

def _short_descr( comments ):
    #First paragraph of the header comments, as a single line.
    words = []
    for line in ( comments or [] ):
        if not line.strip():
            if words:
                break
            continue
        #Skip ascii-art rulers like "----" or "=====":
        if not any( c.isalnum() for c in line ):
            if words:
                break
            continue
        words += line.split()
    return ' '.join(words)

def _fmt_prop( v ):
    if v is None:
        return 'None'
    if isinstance( v, frozenset ):
        return '{' + ','.join( repr(e) for e in sorted(v) ) + '}'
    if isinstance( v, float ):
        return '%g'%v
    return repr(v)

class _WhereExpr:
    #Where expression, compiled after checking that only known names and no
    #private attributes are used (basic safety, and to catch typos early).
    def __init__( self, expr ):
        import ast

        from .exceptions import NCBadInput
        try:
            tree = ast.parse( expr, mode = 'eval' )
        except SyntaxError as e:
            raise NCBadInput(f'Invalid where expression "{expr}":'
                             f' {e.msg}') from e
        allowed = set( n for n,d in _propdocs ) | set( _where_funcs )
        for node in ast.walk( tree ):
            if isinstance( node, ast.Name ) and node.id not in allowed:
                raise NCBadInput(f'Unknown name "{node.id}" in where'
                                 f' expression "{expr}"')
            if isinstance( node, ast.Attribute ) and node.attr.startswith('_'):
                raise NCBadInput(f'Invalid where expression "{expr}" (private'
                                 ' attributes are not allowed)')
        self.expr = expr
        self.__code = compile( tree, '<where>', 'eval' )
    def __call__( self, props ):
        ns = dict( _where_funcs )
        ns['__builtins__'] = {}
        return eval( self.__code, ns, props.as_dict() )

def _eval_where( fct, props ):
    #Comparisons with unavailable values (None) are considered false:
    try:
        return bool( fct( props ) )
    except TypeError:
        return False
    except Exception as e:
        if not isinstance( fct, _WhereExpr ):
            raise
        from .exceptions import NCBadInput
        raise NCBadInput(f'Error evaluating where expression "{fct.expr}":'
                         f' {e}') from e
