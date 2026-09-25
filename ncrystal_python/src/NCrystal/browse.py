
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

"""Utilities for browsing and searching the data available to NCrystal (e.g.
files in the standard data library or in the current directory, in-memory
data, or data created on-demand like "solid::B4C/2.52gcm3"), optionally based
on the physics content of the materials they describe.

The main tool is the DataBrowser class, providing the same features as the
"ncrystal browse" command-line tool. For example:

  import NCrystal.browse as nb
  b = nb.DataBrowser( factory = 'stdlib' )
  b.match('Al*').dump()                   #listing with short descriptions
  sel = b.search('togo').where("'O' in elements and absxs < 0.1")
  print( sel.table('formula,sg,density') )
  for e in sel.sorted('density'):
      print( e.display_name, e.props.density )
  html = sel.to_html('formula,crystalsystem,description')

Physics properties (see physics_props_doc()) are loaded when first needed.
The underlying data is provided by the C++ layer (via JSON queries), which
also loads materials in parallel. It is available directly via the
query_data() function.
"""

__all__ = [ 'DataBrowser', 'DataEntry', 'PhysicsProps', 'list_all_factories',
            'list_factories', 'physics_props_doc', 'query_data' ]

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
        cell = info['cell'] or {}
        ai = info['atominfo']
        def per_atom( key ):
            if ai is None:
                return None
            return dict( (e['label'],e[key]) for e in ai
                         if e[key] is not None )
        self.__d = dict(
            elements = elements,
            atoms = frozenset( atomfracs ),
            nelements = len( elements ),
            formula = format_chemform( sorted( atomfracs.items() ) ),
            mass = info['mass'],
            absxs = info['xsect_absorption'],
            scatxs = info['xsect_free'],
            cohxs = info['xsect_coh'],
            incohxs = info['xsect_incoh'],
            density = info['density'],
            numdens = info['numberdensity'],
            temp = info['temperature'],
            state = info['stateofmatter'].lower(),
            crystal = info['crystalline'],
            crystalsystem = _crystal_system( info['spacegroup'] ),
            sg = info['spacegroup'],
            natoms = info['natoms_unitcell'],
            a = cell.get('a'),
            b = cell.get('b'),
            c = cell.get('c'),
            alpha = cell.get('alpha'),
            beta = cell.get('beta'),
            gamma = cell.get('gamma'),
            volume = cell.get('volume'),
            braggthreshold = info['braggthreshold'],
            debyetemps = per_atom('debyetemp'),
            msds = per_atom('msd'),
            dyninfo = frozenset( info['dyninfo_types'] ),
            customsections = frozenset( info['customsections'] ),
            nphases = info['nphases'],
        )
        self.__compos = compos

    def as_dict( self, json_compatible = False ):
        """The properties as a (new) dictionary. If json_compatible=True,
        sets are replaced with sorted lists."""
        if not json_compatible:
            return dict( self.__d )
        return dict( ( k, sorted(v) if isinstance(v,frozenset) else
                       ( dict(v) if isinstance(v,dict) else v ) )
                     for k,v in self.__d.items() )

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
    """A data entry available to NCrystal (e.g. a file, in-memory data, or
    data created on-demand). Objects are provided by DataBrowser objects."""

    def __init__( self, data, _source = None ):
        """For internal usage only (data is a dictionary from the C++
        browsedb query, and _source is used to load physics on demand)."""
        self.__d = data
        c = data.get('comments')
        self.__comments = tuple(c) if c is not None else None
        self.__props = ( PhysicsProps( data['info'] ) if 'info' in data
                         else None )
        loaded = 'info' in data or 'error' in data
        self.__source = None if loaded else _source

    def _ensure_loaded( self ):
        src, self.__source = self.__source, None
        if src is not None and src.physics is not False:
            e = src.with_physics( [ self ] )[0]
            self.__d, self.__props = e.__d, e.__props

    @property
    def name( self ):
        """Name (e.g. a file name) used to request the entry."""
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
    def path( self ):
        """Absolute path of the on-disk file with the data, or None if not
        on disk (e.g. in-memory data or data embedded in the library)."""
        return self.__d.get('path')

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
        """PhysicsProps of the loaded material, or None if the material could
        not be loaded (see .error). The material is loaded on first access
        if needed (unless the DataBrowser was created with physics=False, in
        which case this is always None)."""
        self._ensure_loaded()
        return self.__props

    @property
    def error( self ):
        """Error message if loading the material failed, otherwise None
        (loads the material if needed, as for .props)."""
        self._ensure_loaded()
        return self.__d.get('error')

    def matching_lines( self, search, *, regex = False ):
        """Lines of the header comments matching any of the search words
        (case-insensitive, or regular expressions if regex=True)."""
        res = _compile_words( search, regex )
        return [ ll.strip() for ll in ( self.__comments or [] )
                 if any( r.search(ll) for r in res ) ]

    def as_dict( self ):
        """JSON-compatible dictionary with all information about the
        entry (props is None if not loaded)."""
        p = self.__props
        return dict( name = self.name, fullkey = self.fullkey,
                     factory = self.factory, source = self.source,
                     priority = self.priority, hidden = self.hidden,
                     datatype = self.datatype, path = self.path,
                     description = self.description,
                     comments = ( list(self.__comments)
                                  if self.__comments is not None else None ),
                     props = ( p.as_dict( json_compatible = True )
                               if p is not None else None ),
                     error = self.error )

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

def query_data( factory = None, *, physics = True, nthreads = 'auto',
                quiet = True, as_json = False, progress = None ):
    """Returns the raw data about all available entries (or those from a given
    factory) as a list of dictionaries, sorted by priority (highest first),
    factory, source, and name. If as_json=True, a JSON string is returned
    instead.

    If physics=True, all materials are loaded and their physics properties
    included (in "info" dicts, or an "error" message if loading failed). This
    is done in parallel with nthreads threads (an integer, or "auto"), which
    are only used temporarily, and only if the user did not already configure
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
    data = _query( facts, counts, physics = physics, nthreads = nthreads,
                   quiet = quiet, progress = progress )
    if as_json:
        import json
        return json.dumps( data )
    return data

class DataBrowser:
    """Browse and search the data available to NCrystal, in the same way as
    with the "ncrystal browse" command-line tool.

    Selection methods (match, search, where, filter, from_factory, sorted,
    and slicing) always return new DataBrowser objects, so they can be
    chained. The selected entries (DataEntry objects) are available via
    iteration, indexing, or .entries, and can be shown with methods like
    dump(), table(), info(), to_csv(), to_json(), and to_html().

    Physics properties (see physics_props_doc()) of the selected materials
    are loaded when first needed (e.g. by where(..) or tables of physics
    properties), unless physics=True which loads all materials immediately,
    or physics=False which disables loading. The nthreads, quiet and
    progress parameters control the loading (see query_data()).
    """

    def __init__( self, *patterns, factory = None, search = (),
                  regex = False, where = (), physics = None,
                  nthreads = 'auto', quiet = True, progress = None ):
        """Browse all available data, or that of a given factory. Any
        patterns, search words, and where conditions are applied as with the
        corresponding methods."""
        src = _Source( factory, physics = physics, nthreads = nthreads,
                       quiet = quiet, progress = progress )
        self._init( src, src.entries(), (), False )
        b = self
        if physics:
            b = b.with_physics()
        if patterns:
            b = b.match( *patterns, regex = regex )
        if _aslist( search ):
            b = b.search( *_aslist(search), regex = regex )
        if _aslist( where ):
            b = b.where( *_aslist(where) )
        self._init( b._src, b._entries, b._search, b._searchregex )

    def _init( self, src, entries, search, searchregex ):
        self._src = src
        self._entries = tuple( entries )
        self._search = tuple( search )
        self._searchregex = searchregex

    def _derived( self, entries, search = None, searchregex = None ):
        b = object.__new__( DataBrowser )
        b._init( self._src, entries,
                 self._search if search is None else search,
                 self._searchregex if searchregex is None else searchregex )
        return b

    #Selection methods:

    def match( self, *patterns, regex = False ):
        """Select entries whose names match at least one of the patterns
        (case-insensitive substrings, or glob patterns if they contain "*",
        "?", or "["). Patterns containing "::" are matched against the full
        key (e.g. "stdlib::Al*"). If regex=True, patterns are instead
        (case-insensitive) Python regular expressions."""
        if not patterns:
            return self
        prs = ( _compile_words( patterns, True ) if regex
                else [ None ]*len(patterns) )
        return self._derived( [ e for e in self._entries
                                if any( _name_matches( e, p, pr )
                                        for p, pr in zip( patterns, prs ) ) ] )

    def search( self, *words, regex = False ):
        """Select entries containing all the words in their name or NCMAT
        header comments (case-insensitive, or Python regular expressions if
        regex=True, where ^ and $ match at line boundaries). The words are
        remembered, so matching comment lines are shown by dump()."""
        res = _compile_words( words, regex )
        sel = [ e for e in self._entries
                if all( r.search( '\n'.join( [ e.name ]
                                             + list( e.comments or [] ) ) )
                        for r in res ) ]
        if self._search and self._searchregex != regex:
            #Mixed modes: store all as regular expressions:
            import re
            prev = ( self._search if self._searchregex else
                     tuple( re.escape(w) for w in self._search ) )
            new = words if regex else tuple( re.escape(w) for w in words )
            return self._derived( sel, prev + tuple(new), True )
        return self._derived( sel, self._search + tuple(words),
                              regex or self._searchregex )

    def where( self, *conditions ):
        """Select entries whose loaded materials fulfil all the conditions,
        which are either functions taking a PhysicsProps object, or Python
        expressions (strings) using the property names (e.g.
        "'B' in elements and absxs > 100"). Expressions failing due to
        unavailable (None) values are considered false, and entries which
        could not be loaded are never selected."""
        fcts = [ ( w if callable(w) else _WhereExpr(w) ) for w in conditions ]
        if not fcts:
            return self
        b = self.with_physics()
        return b._derived( [ e for e in b._entries if e.props is not None
                             and all( _eval_where( f, e.props )
                                      for f in fcts ) ] )

    def filter( self, fct ):
        """Select entries for which fct(entry) is true, where entry is a
        DataEntry object."""
        return self._derived( [ e for e in self._entries if fct(e) ] )

    def from_factory( self, name ):
        """Select entries delivered by the named factory."""
        return self._derived( [ e for e in self._entries
                                if e.factory == name ] )

    def sorted( self, key, *, reverse = False ):
        """Sort entries by the given key, which is either "name", the name of
        a physics property (see physics_props_doc()), or a function taking a
        DataEntry object. Entries without a value (None, e.g. materials which
        could not be loaded) are always placed last, and ties keep their
        order. Sets are compared by their sorted contents."""
        if callable( key ):
            value = key
        else:
            names = [ n for n,d in _propdocs if n not in _dict_props ]
            if key != 'name' and key not in names:
                from .exceptions import NCBadInput
                raise NCBadInput(f'Invalid sort key "{key}" (must be "name" or'
                                 f' one of: {", ".join(names)})')
            def value( e ):
                if key == 'name':
                    return e.display_name
                if e.props is None:
                    return None
                v = getattr( e.props, key )
                return tuple(sorted(v)) if isinstance( v, frozenset ) else v
        b = self if ( callable(key) or key == 'name' ) else self.with_physics()
        have = [ e for e in b._entries if value(e) is not None ]
        missing = [ e for e in b._entries if value(e) is None ]
        return b._derived( sorted( have, key = value, reverse = reverse )
                           + missing )

    def with_physics( self ):
        """Returns browser with physics properties of all entries loaded
        (entries which could not be loaded have an .error instead)."""
        return self._derived( self._src.with_physics( self._entries ) )

    #Access to entries:

    @property
    def entries( self ):
        """Tuple of the selected DataEntry objects."""
        return self._entries

    def names( self ):
        """List of names (the full keys of hidden entries)."""
        return [ e.display_name for e in self._entries ]

    def factories( self ):
        """List of factories of the selected entries (in order)."""
        res = []
        for e in self._entries:
            if e.factory not in res:
                res.append( e.factory )
        return res

    def __len__( self ):
        return len( self._entries )

    def __bool__( self ):
        return bool( self._entries )

    def __iter__( self ):
        return iter( self._entries )

    def __getitem__( self, idx ):
        if isinstance( idx, slice ):
            return self._derived( self._entries[idx] )
        return self._entries[idx]

    def __repr__( self ):
        n, nf = len(self), len( self.factories() )
        return ( f'DataBrowser({n} entr{"y" if n==1 else "ies"} from {nf}'
                 f' factor{"y" if nf==1 else "ies"})' )

    def suggestions( self, pattern ):
        """Names of entries similar to pattern (e.g. to help with typos)."""
        return _suggestions( self._entries, pattern )

    #Output (NB: dump methods print, other methods return strings):

    def format_listing( self, *, comments = False, props = False,
                        truncate = True, linewidth = 80, highlight = None ):
        """Returns listing of the entries grouped by source, with short
        descriptions, and matching comment lines of any search words. If
        comments=True, full NCMAT header comments are shown instead of
        descriptions, and if props=True the physics properties are shown.
        Text is shortened to the linewidth unless truncate=False. If
        highlight is provided, it is a function applied to names and texts
        (e.g. adding color codes around search hits)."""
        b = self.with_physics() if props else self
        return _format_listing( b, comments = comments, props = props,
                                truncate = truncate, linewidth = linewidth,
                                hl = highlight )

    def dump( self, **kwargs ):
        """Print the result of format_listing(**kwargs)."""
        from ._common import print
        print( self.format_listing( **kwargs ), end = '' )

    def table( self, columns = 'description', *, fmt = 'text',
               truncate = True, linewidth = 80 ):
        """Returns table with the names and the given columns (a list, or a
        comma-separated string, of physics properties and "description").
        The format (fmt) is "text", "csv" (with full numerical precision),
        "json", or "html". Text tables are never truncated, except for the
        description column (unless truncate=False)."""
        cols = self._columns( columns )
        b = ( self.with_physics() if any( c != 'description' for c in cols )
              else self )
        if fmt == 'text':
            return _format_table( b, cols, truncate = truncate,
                                  linewidth = linewidth )
        if fmt == 'csv':
            return _format_csv( b, cols )
        if fmt == 'json':
            import json
            return json.dumps( [ _table_json( e, cols ) for e in b ],
                               indent = 1 ) + '\n'
        if fmt == 'html':
            return _format_html( b, cols )
        from .exceptions import NCBadInput
        raise NCBadInput(f'Invalid table format: "{fmt}" (must be "text",'
                         ' "csv", "json", or "html")')

    def to_csv( self, columns = 'description' ):
        """Same as table(columns,fmt="csv")."""
        return self.table( columns, fmt = 'csv' )

    def to_html( self, columns = 'description' ):
        """Same as table(columns,fmt="html")."""
        return self.table( columns, fmt = 'html' )

    def to_dicts( self ):
        """List of dictionaries (JSON compatible) with all information about
        the entries, including physics properties."""
        return [ e.as_dict() for e in self.with_physics() ]

    def to_json( self ):
        """JSON string with all information about the entries (as
        to_dicts())."""
        import json
        return json.dumps( self.to_dicts(), indent = 1 ) + '\n'

    def info( self, *, linewidth = 80 ):
        """Returns all available information about the entries, including
        physics properties, header comments, and usage examples."""
        return _format_info( self.with_physics(), linewidth = linewidth )

    def dump_info( self, **kwargs ):
        """Print the result of info(**kwargs)."""
        from ._common import print
        print( self.info( **kwargs ), end = '' )

    def _columns( self, columns ):
        cols = ( [ c.strip() for c in columns.split(',') if c.strip() ]
                 if isinstance( columns, str ) else list( columns ) )
        allowed = [ n for n,d in _propdocs ] + ['description']
        for c in cols:
            if c not in allowed:
                from .exceptions import NCBadInput
                raise NCBadInput(f'Invalid column "{c}" (must be "description"'
                                 f' or one of: {", ".join(allowed[:-1])})')
        return cols

###############################################################################
# Implementation details:

_progress_chunk = 16

_propdocs = [
    ('elements', 'set of element names, e.g. {"Al","O"}'),
    ('atoms', 'set of atom labels, including isotopes, e.g. {"D","O"}'),
    ('nelements', 'number of different elements'),
    ('formula', 'chemical formula, e.g. "Al2O3"'),
    ('mass', 'average atomic mass [amu]'),
    ('absxs', 'absorption cross section per atom at 2200m/s [barn]'),
    ('scatxs', 'free scattering cross section per atom [barn]'),
    ('cohxs', 'bound coherent scattering cross section per atom [barn]'),
    ('incohxs', 'bound incoherent scattering cross section per atom [barn]'),
    ('density', 'density [g/cm3]'),
    ('numdens', 'number density [atoms/Aa3]'),
    ('temp', 'temperature [K]'),
    ('state', 'state of matter: "solid", "liquid", "gas", or "unknown"'),
    ('crystal', 'True if crystalline (for multiphase: any phase)'),
    ('crystalsystem', ('crystal system (e.g. "cubic"), based on the space'
                       ' group (None if not available)')),
    ('sg', 'space group number (None if not available)'),
    ('natoms', 'number of atoms in unit cell (None if not available)'),
    ('a', 'unit cell length a [Aa] (None if not available)'),
    ('b', 'unit cell length b [Aa] (None if not available)'),
    ('c', 'unit cell length c [Aa] (None if not available)'),
    ('alpha', 'unit cell angle alpha [degree] (None if not available)'),
    ('beta', 'unit cell angle beta [degree] (None if not available)'),
    ('gamma', 'unit cell angle gamma [degree] (None if not available)'),
    ('volume', 'unit cell volume [Aa^3] (None if not available)'),
    ('braggthreshold', ('Bragg threshold [Aa], i.e. the longest wavelength'
                        ' with Bragg diffraction (None if not available)')),
    ('debyetemps', ('dict of per-atom Debye temperatures [K], e.g.'
                    ' {"Al":412.2} (None if not a crystal)')),
    ('msds', ('dict of per-atom mean-squared-displacements [Aa^2] (None if'
              ' not a crystal)')),
    ('dyninfo', ('set of dynamic info types present, among "vdos",'
                 ' "vdosdebye", "scatknl" (a full scattering kernel),'
                 ' "freegas", and "sterile"')),
    ('customsections', 'set of names of @CUSTOM_ sections in NCMAT data'),
    ('nphases', 'number of phases (1 for single-phase materials)'),
]

def _crystal_system( sg ):
    if not sg:
        return None
    for sgmin, name in ( (195,'cubic'), (168,'hexagonal'), (143,'trigonal'),
                         (75,'tetragonal'), (16,'orthorhombic'),
                         (3,'monoclinic'), (1,'triclinic') ):
        if sg >= sgmin:
            return name

_dict_props = ('debyetemps','msds')#not sortable

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

def _sortkey( d ):
    #Sort key for dicts from browsedb query:
    p = d['priority']
    rank = p if isinstance( p, int ) else ( -1 if p == 'OnlyOnExplicitRequest'
                                            else -2 )
    return ( -rank, d['factory'], d['source'], d['name'] )

def _query( facts, counts, *, physics, nthreads, quiet, progress ):
    #Query given factories (with or without physics), returning sorted list
    #of dicts:
    ntotal = sum( counts[f] for f in facts )
    extra = ( ( f'nthreads={nthreads}', ) if physics else ('cheap',) )
    data, ndone = [], 0
    from ._msg import _suppress_msgs_ctx
    msgctx = _suppress_msgs_ctx() if ( physics and quiet ) else _nullctx()
    with msgctx:
        for f in facts:
            #Smaller chunks if progress reporting is needed:
            nchunks = ( 1 if progress is None
                        else max( 1, -(-counts[f]//_progress_chunk) ) )
            for i in range( nchunks ):
                if progress is not None:
                    progress( ndone, ntotal )
                d = _q( 'browsedb', f, str(i), str(nchunks), *extra )
                data += d
                ndone += len(d)
    if progress is not None:
        progress( ndone, ntotal )
    data.sort( key = _sortkey )
    return data

class _Source:
    #Data shared by a DataBrowser and all browsers derived from it, in
    #particular caching physics properties loaded on demand (per factory).
    def __init__( self, factory, *, physics, nthreads, quiet, progress ):
        self.counts = list_factories()
        if factory is not None and factory not in self.counts:
            from .exceptions import NCBadInput
            raise NCBadInput(f'Unknown TextData factory: "{factory}"')
        self.facts = [ f for f,n in self.counts.items()
                       if n and ( factory is None or f == factory ) ]
        self.physics, self.nthreads = physics, nthreads
        self.quiet, self.progress = quiet, progress
        self.cheap = None
        self.loaded = {}
        self.loaded_facts = set()

    def entries( self ):
        if self.cheap is None:
            self.cheap = [ DataEntry( d, self ) for d in
                           _query( self.facts, self.counts, physics = False,
                                   nthreads = 1, quiet = True,
                                   progress = None ) ]
        return self.cheap

    def with_physics( self, entries ):
        if self.physics is False:
            from .exceptions import NCBadInput
            raise NCBadInput('Physics properties are not available, since'
                             ' the DataBrowser was created with'
                             ' physics=False')
        needed = []
        for e in entries:
            if e.factory not in self.loaded_facts and e.factory not in needed:
                needed.append( e.factory )
        if needed:
            for d in _query( needed, self.counts, physics = True,
                             nthreads = self.nthreads, quiet = self.quiet,
                             progress = self.progress ):
                self.loaded[d['fullkey']] = DataEntry(d)
            self.loaded_facts.update( needed )
        return [ self.loaded.get( e.fullkey, e ) for e in entries ]

def _truncate( s, n, truncate = True ):
    if not truncate or len(s) <= n:
        return s
    return s[:max(0,n-3)].rstrip() + '...'

def _strip_empty( lines ):
    lines = list( lines or [] )
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()
    return lines

def _noop( s ):
    return s

def _props_lines( entry, linewidth, indent, truncate ):
    pre = ' '*indent
    if entry.props is None:
        return [ _truncate( f'{pre}[could not load: {entry.error}]',
                            linewidth, truncate ) ]
    import textwrap
    parts = [ f'{k}={_fmt_prop(v)}' for k,v in entry.props.as_dict().items() ]
    return [ pre + ll for ll in textwrap.wrap( '  '.join(parts),
                                               width = linewidth - indent,
                                               break_long_words = False,
                                               break_on_hyphens = False ) ]

def _format_listing( b, *, comments, props, truncate, linewidth, hl ):
    hl = hl or _noop
    out = []
    groups = []
    for e in b:
        key = ( e.factory, e.source, e.priority )
        if not groups or groups[-1][0] != key:
            groups.append( ( key, [] ) )
        groups[-1][1].append( e )
    for (factname, source, priority), group in groups:
        n = len(group)
        src = f' ({source}, priority={priority})' if source else (
            f' (priority={priority})' )
        out.append(f'==> {n} entr{"y" if n==1 else "ies"} from "{factname}"'
                   f'{src}:')
        namew = min( 40, max( len(e.display_name) for e in group ) )
        for e in group:
            name = e.display_name
            extra = ' (hidden)' if e.hidden else ''
            descr = e.description
            #NB: Truncate before highlighting, so color codes do not count:
            padding = ' '*( max(0,namew-len(name)) )
            if descr and not comments:
                room = linewidth - 4 - max(namew,len(name)) - 2 - len(extra)
                line = ( f'    {hl(name)}{padding}  '
                         f'{hl(_truncate(descr,room,truncate))}{extra}' )
            else:
                line = f'    {hl(name)}{extra}'
            out.append( line.rstrip() )
            if b._search and not comments:
                for ll in e.matching_lines( b._search,
                                            regex = b._searchregex ):
                    out.append( '        | '
                                + hl(_truncate(ll,linewidth-10,truncate)) )
            if props:
                out += _props_lines( e, linewidth, 8, truncate )
            if comments:
                cl = _strip_empty( e.comments )
                out += [ f'        # {hl(ll)}'.rstrip() for ll in cl ]
                if cl:
                    out.append('')
    return ''.join( ll + '\n' for ll in out )

def _table_value( entry, col ):
    if col == 'description':
        return entry.description
    if entry.props is None:
        return '-'
    v = getattr( entry.props, col )
    if v is None:
        return '-'
    if isinstance( v, frozenset ):
        return ','.join( sorted(v) ) if v else '-'
    if isinstance( v, dict ):
        return ','.join( f'{k}:{x:g}' for k,x in v.items() ) if v else '-'
    if isinstance( v, float ):
        return '%g'%v
    return str(v)

def _table_json( entry, cols ):
    d = dict( name = entry.display_name )
    props = ( entry.props.as_dict( json_compatible = True )
              if entry.props is not None else None )
    for c in cols:
        if c == 'description':
            d[c] = entry.description
        else:
            d[c] = props[c] if props is not None else None
    return d

def _format_table( b, cols, *, truncate, linewidth ):
    #NB: Rows are never truncated (no data should be lost), except for the
    #free-text description column which is shortened to fit if possible.
    header = [ 'NAME' ] + [ c.upper() for c in cols ]
    rows = [ [ e.display_name ] + [ _table_value(e,c) for c in cols ]
             for e in b ]
    def width( k ):
        return max( len(r[k]) for r in rows + [header] )
    if 'description' in cols:
        kd = 1 + cols.index('description')
        other = sum( width(k) + 2 for k in range(len(header)) if k != kd )
        room = max( 20, linewidth - other )
        for r in rows:
            r[kd] = _truncate( r[kd], room, truncate )
    widths = [ width(k) for k in range(len(header)) ]
    def fmt( r ):
        return '  '.join( ( v.ljust(w) if k==0 or cols[k-1]=='description'
                            else v.rjust(w) )
                          for k,(v,w) in enumerate(zip(r,widths)) ).rstrip()
    return ''.join( fmt(r) + '\n' for r in [header] + rows )

def _csv_value( v ):
    #Full precision (repr) floats, empty for unavailable values:
    if v is None:
        return ''
    if isinstance( v, list ):
        return ','.join( str(e) for e in v )
    if isinstance( v, dict ):
        return ','.join( f'{k}:{x!r}' for k,x in v.items() )
    if isinstance( v, float ):
        return repr(v)
    return str(v)

def _format_csv( b, cols ):
    import csv
    import io
    buf = io.StringIO()
    w = csv.writer( buf, lineterminator = '\n' )
    w.writerow( [ 'name' ] + cols )
    for e in b:
        d = _table_json( e, cols )
        w.writerow( [ d['name'] ] + [ _csv_value(d[c]) for c in cols ] )
    return buf.getvalue()

def _format_html( b, cols ):
    import html
    def row( cells, tag ):
        return ( '<tr>' + ''.join( f'<{tag}>{html.escape(c)}</{tag}>'
                                   for c in cells ) + '</tr>\n' )
    out = [ '<table class="ncrystal-browse">\n<thead>\n',
            row( [ 'name' ] + cols, 'th' ), '</thead>\n<tbody>\n' ]
    for e in b:
        vals = [ ( '' if v == '-' else v )
                 for v in ( _table_value(e,c) for c in cols ) ]
        out.append( row( [ e.display_name ] + vals, 'td' ) )
    out.append( '</tbody>\n</table>\n' )
    return ''.join( out )

def _format_info( b, *, linewidth ):
    import textwrap
    out = []
    for n, e in enumerate( b ):
        if n:
            out.append('')
        name = e.display_name
        out.append(f'==> {name}')
        srcdescr = ( f'{e.source} (factory "{e.factory}",'
                     f' priority {e.priority})' )
        fields = [ ('Description', e.description or '-'),
                   ('Full key', e.fullkey),
                   ('Source', srcdescr),
                   ('Data type', e.datatype or '-'),
                   ('On-disk path', e.path or '-') ]
        if e.hidden:
            fields.append( ('Hidden', ( 'yes, by a higher priority entry'
                                        ' with the same name' ) ) )
        for k, v in fields:
            ll = textwrap.wrap( v, width = linewidth - 20,
                                break_long_words = False,
                                break_on_hyphens = False ) or ['']
            out.append(f'    {k:<14}: {ll[0]}')
            out += [ f'{"":20}{x}' for x in ll[1:] ]
        out.append('  Physics properties:')
        out += _props_lines( e, linewidth, 4, True )
        comments = _strip_empty( e.comments )
        if comments:
            out.append('  Header comments:')
            out += [ f'    # {ll}'.rstrip() for ll in comments ]
        out.append('  Usage examples:')
        pycmd = ( f'python3 -c \'import NCrystal as NC;'
                  f' NC.load("{name}").dump()\'' )
        for d, c in [ ( 'Plot cross sections', f'nctool "{name}"' ),
                      ( 'Show full content',
                        f'ncrystal browse -x "{name}"' ),
                      ( 'Load in Python', pycmd ) ]:
            out.append(f'    # {d}:')
            out.append(f'    {c}')
    return ''.join( ll + '\n' for ll in out )

def _suggestions( entries, pattern ):
    #Names of entries similar to the pattern:
    import difflib
    def stem( n ):
        return n.rsplit('.',1)[0].lower() if '.' in n else n.lower()
    #Compare with full names (without extension), and with their parts
    #(e.g. "diamond" in "C_sg227_Diamond.ncmat"):
    cands = {}
    for e in entries:
        st = stem(e.name)
        cands.setdefault( st, e.name )
        for part in st.split('_'):
            if len(part) >= 4:
                cands.setdefault( part, e.name )
    res = []
    for m in difflib.get_close_matches( stem(pattern), list(cands), n = 3,
                                        cutoff = 0.7 ):
        if cands[m] not in res:
            res.append( cands[m] )
    return res

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
    if isinstance( v, dict ):
        return '{' + ','.join( f'{k!r}:{_fmt_prop(x)}'
                               for k,x in v.items() ) + '}'
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
    #Expressions failing due to unavailable values (None) are considered
    #false (e.g. "sg > 200" or "max(debyetemps.values()) > 500"):
    try:
        return bool( fct( props ) )
    except TypeError:
        return False
    except AttributeError as e:
        if "'NoneType' object" in str(e):
            return False
        if not isinstance( fct, _WhereExpr ):
            raise
        from .exceptions import NCBadInput
        raise NCBadInput(f'Error evaluating where expression "{fct.expr}":'
                         f' {e}') from e
    except Exception as e:
        if not isinstance( fct, _WhereExpr ):
            raise
        from .exceptions import NCBadInput
        raise NCBadInput(f'Error evaluating where expression "{fct.expr}":'
                         f' {e}') from e
