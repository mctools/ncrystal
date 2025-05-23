
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

dbfilerelpath = 'devel/plugin_database.yml'

def short_description():
    return f'Extract information about plugins from {dbfilerelpath}'

def load_and_check_data():
    from .dirs import reporoot
    import yaml
    with reporoot.joinpath(dbfilerelpath).open() as fh:
        data =yaml.safe_load(fh)
    errmsg = _prepare_data_and_check_errors( data )
    if errmsg:
        raise SystemExit(f'ERROR invalid data in {dbfilerelpath} : {errmsg}')
    return data

def _add_piptargets_srcurls_etc( pluginname, info ):

    pypkgname = f'ncrystal-plugin-{pluginname}'
    repo_baseurl = 'https://github.com/%s'%info['github_repo_key']
    if not info['gitref'] and not info['repo_subdir']:
        src_url = repo_baseurl
    else:
        resolved_gitref = info['gitref'] or 'HEAD'
        src_url = repo_baseurl + f'/tree/{resolved_gitref}'
        if info['repo_subdir']:
            src_url += '/%s'%info['repo_subdir']

    if info['is_on_pypi']:
        pip_target = pypkgname
    else:
        pip_target = f'git+{repo_baseurl}'
        if info['gitref']:
            pip_target+='@%s'%info['gitref']
        if info['repo_subdir']:
            pip_target += '#subdirectory=%s'%info['repo_subdir']

    #shlex.quote assumes posix shells, and do not really work on windows,
    #where we have to use mslex.quote. However, it seems that simply adding
    #a " at each end is OK on both platforms in this very specific case:
    quoted_pip_target = ( f'"{pip_target}"'
                          if any( e in pip_target for e in ' @#+')
                          else  pip_target )

    pypi_url = ( f'https://pypi.org/project/{pypkgname}'
                 if info['is_on_pypi'] else '' )

    info['pypkgname'] = pypkgname
    info['src_url'] = src_url
    info['pip_target'] = pip_target
    info['quoted_pip_target'] = quoted_pip_target
    info['pypi_url'] = pypi_url

def _prepare_data_and_check_errors( data ):
    #Prepares data in place, returns an error message in case of issues
    if not isinstance(data,dict):
        return "top-level info must be a dictionary"
    known_keys = set(['github_repo_key','gitref','repo_subdir','description',
                      'disable_tests', 'is_on_pypi' ])
    required_keys = set(['github_repo_key','description'])
    boolean_keys = set(['disable_tests', 'is_on_pypi'])
    assert not ( boolean_keys - known_keys )
    assert not (required_keys-known_keys)
    import textwrap
    for k0, v0 in data.items():
        if not isinstance(k0,str):
            return "top-level keys must be strings"
        if not isinstance(v0,dict):
            return f'data base info for "{k0}" is not a dictionary'
        keys = set(v0.keys())
        missing = required_keys - keys
        if missing:
            return ( f'data base info for "{k0}" is '
                     f'missing required key: "{missing.pop()}"' )
        unknown = keys - known_keys
        if unknown:
            return ( f'data base info for "{k0}" has '
                     f'unknown key: "{unknown.pop()}"' )
        for k1, v1 in v0.items():
            _type,_type_str = str, 'string'
            if k1 in boolean_keys:
                _type,_type_str = bool, 'bool'
            if not isinstance(v1,_type):
                return f'value for {k0}.{k1}={repr(v1)} is not a {_type_str}'
        assert 'description' in v0
        v0['description'] = textwrap.fill(v0['description'],9999999)
        if v0['github_repo_key'].count('/')!=1:
            return f'invalid value for {k0}.github_repo_key'
        assert 'repourl' not in v0
        v0['repourl'] = 'https://github.com/%s'%v0['github_repo_key']
        #Add some defaults:
        for bk in boolean_keys:
            #Booleans should be marked as '1' if true and empty string '' else:
            v0[bk] = '1' if v0.get(bk) else ''
        #All known keys are set to empty strings if absent:
        for k in known_keys:
            if k not in v0:
                v0[k] = ''
        _add_piptargets_srcurls_etc(k0, v0)

def main( parser ):
    parser.init( short_description() )
    parser.add_argument(
        '-m','--mode', type = str, choices = ('list',
                                              'listjson',
                                              'json',
                                              'wiki',
                                              'pprint'),
        default = None,
        help="""Mode (default: pprint unless --extract is used)."""
    )

    parser.add_argument(
        '-e','--extract', nargs='+',
        help="""List of keys to extract values from."""
    )

    args = parser.parse_args()
    if args.extract:
        if args.mode is not None:
            parser.error('Incompatible flags')
        args.mode = 'extract'

    if not args.mode:
        parser.error('No mode selected')

    data = load_and_check_data()
    if args.mode=='pprint':
        import pprint
        pprint.pprint( data )
        return

    if args.mode=='wiki':
        print(_produce_wiki( data ),end='')
        return

    if args.mode=='json':
        import json
        print( json.dumps( data ) )
        return

    if args.mode=='list':
        for k,v in sorted(data.items()):
            print(k)
        return

    if args.mode=='listjson':
        import json
        print( json.dumps(sorted(data.keys())) )
        return

    assert args.mode == 'extract'
    keys = args.extract[:]
    while keys:
        k = keys.pop(0)
        if k not in data:
            if hasattr(data,'keys'):
                c = '", "'.join(sorted(data.keys()))
                parser.error(f'Key not found: "{k}" '
                             f'(valid choices would be "{c}")')
            else:
                parser.error(f'Key not found: "{k}" '
                             '(item at this position is not a dictionary)')
        data = data[k]
    print(data)

def _www_plugin_sort( e ):
    pluginname, info = e
    return (
        int(bool('Dummy' in pluginname)),
        -int(bool(info['is_on_pypi'])),
        int(bool(info['disable_tests'])),
        pluginname,
        info
        )

def _produce_wiki( data ):
    from builtins import print as _orig_print
    import io
    out_buf = io.StringIO()
    def print(*a,**kw):
        return _orig_print(*a, file=out_buf, *kw)
    def printdd(a):
        import textwrap
        print(textwrap.dedent(a))
    printdd("""
    # Curated list of plugins

    On this page we maintain a curated list of plugins, known to the NCrystal
    developers and gauged to be of sufficient interest for the wider
    community. If we are missing something, please [[Contact]] us to add another
    plugin to the list (or simply open a PR adding your plugin to the list we
    maintain in [this file](../blob/main/devel/plugin_database.yml)).

    For how to use plugins in general, consult the [[Plugins]] page.

    """)
    for pluginname, info in sorted(data.items(),
                                   key = _www_plugin_sort ):
        printdd(f"""
        ## {pluginname}

        {info['description']}
        """)
        if info['disable_tests']:
            print(f"* **WARNING:** `ncrystal-pluginmanager --test {pluginname}`"
                  " is currently failing for this plugin")
        print("* **Install with:** `pip install %s`"%info['quoted_pip_target'])
        if info['is_on_pypi']:
            print("* **PyPI URL:** [%s](%s)"%(info['pypi_url'],
                                              info['pypi_url']))
        print("* **Source URL:** [%s](%s)"%(info['src_url'],
                                            info['src_url']))
    return out_buf.getvalue().rstrip()+'\n'
