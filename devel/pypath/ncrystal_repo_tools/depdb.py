
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

dbfilerelpath = 'devel/dependency_database.yml'

def load_depdb():
    from .dirs import reporoot
    import yaml
    with reporoot.joinpath(dbfilerelpath).open() as fh:
        data = yaml.safe_load(fh)
    errmsg, data = _prepare_data_and_check_errors( data )
    if errmsg:
        raise SystemExit(f'ERROR invalid data in {dbfilerelpath} : {errmsg}')
    return data

def load_part2deplist_from_pyproject_toml( tomlfile ):
    from .toml import parse_toml
    t = parse_toml(tomlfile)
    part2deplist = {}
    part2deplist['BASE'] = t['project'].get('dependencies',[])
    for part,deplist in t['project'].get('optional-dependencies',{}).items():
        part2deplist[ part ] = deplist
    return part2deplist

def _handle_ADD( pkgname, part2deplist ):
    def has_add( deplist ):
        return any( d.startswith('ADD') for d in deplist )
    ready = dict( (k,v) for k,v in part2deplist.items() if not has_add(v) )
    pending = dict( (k,v) for k,v in part2deplist.items() if has_add(v) )
    i = 0
    while pending:
        i+=1
        if i==1000:
            return 'recursive ADD statements involving %s'%pending.keys()
        for partname,deplist in pending.items():
            for i in range(len(deplist)):
                if deplist[i].startswith('ADD'):
                    needed_part = deplist[i][len('ADD'):].strip()
                    if needed_part == partname:
                        raise SystemExit('ERROR: Self referencing ADD '
                                         f'statement involving {partname}')
                    if needed_part in ready:
                        deplist[i] = ''#to be removed later
                        deplist += ready[needed_part]
            if not has_add(deplist):
                #transfer copy to 'ready'
                ready[partname] = deplist[:]
                del deplist[:]#clear in place
        for partname in list(pending.keys()):
            if partname in ready:
                assert len(pending[partname]) == 0
                del pending[partname]
    assert set(part2deplist.keys()) == set(ready.keys())
    for k,v in ready.items():
        part2deplist[k] = v

_cmdcache=[{}]
def _handle_CMD( pkgname, part2deplist ):
    from .dirs import reporoot
    import subprocess
    import json
    import sys
    for part,deplist in sorted(part2deplist.items()):
        for i in range(len(deplist)):
            if not deplist[i].startswith('CMD'):
                continue
            cmd = []
            for e in deplist[i][len('CMD'):].split():
                if e=='<python>':
                    e = sys.executable
                elif e.startswith('<reporoot>/'):
                    e = str(reporoot.joinpath(*(e.split('/')[1:])))
                cmd.append(e)
            cmd = tuple(str(e) for e in cmd)
            if cmd not in _cmdcache[0]:
                rv = subprocess.run(cmd,capture_output=True)
                if rv.stderr or rv.returncode != 0:
                    import shlex
                    raise SystemExit('CMD failed: %s'%shlex.join(cmd))
                _cmdcache[0][cmd] = json.loads(rv.stdout.decode())
            deplist[i]=''#clear previous
            deplist += _cmdcache[0][cmd]

def append_version_var( s, val, operator, sep='' ):
    if not val:
        return s
    return s + sep + ( val if operator in val else f'{operator}={val}' )

def create_py_spec(info):
    pn = info['pypi_pkgname']
    min_v = info.get('min_version')
    max_v = info.get('max_version')
    assert not (min_v and max_v)
    pystd = info.get('in_pystdlib_from')
    s = append_version_var( pn, min_v, '>' )
    s = append_version_var( s, max_v, '<' )
    if pystd:
        s += f'; python_version < "{pystd}"'
    return s

def create_conda_spec(info):
    pn = info['pypi_pkgname']
    min_v = info.get('min_version')
    max_v = info.get('max_version')
    assert not (min_v and max_v)
    s = append_version_var( pn, min_v, '>', ' ' )
    s = append_version_var( s, max_v, '<', ' ' )
    return s

def _prepare_data_and_check_errors( rawdata ):
    from .dirs import reporoot
    #version = reporoot.joinpath('VERSION').read_text().strip()

    assert set(rawdata.keys())==set(['packages','dependencies'])
    #First pkg_deps:

    pkgdata = rawdata['packages']
    for pkgname,info in sorted(pkgdata.items()):
        if not info.get('deps'):
            info['deps'] = {}
        if not info['deps'].get('BASE'):
            info['deps']['BASE'] = []
        assert set(info.keys())==set(['tomlfiles','deps'])
        assert len(info['tomlfiles'])>=1
        part2deplist = info['deps']
        errmsg = ( _handle_CMD(pkgname,part2deplist)
                   or _handle_ADD(pkgname,part2deplist) )
        if errmsg:
            return errmsg, None
        for tf in info['tomlfiles']:
            if not reporoot.joinpath(*tf.split('/')).is_file():
                raise SystemExit(f'File not found: {tf}')

    def check_name_ok(p):
        if not ( isinstance(p,str)
                 and p
                 and p.replace('-','_').isidentifier() ):
            raise SystemExit(f'Invalid name : {p}')

    all_deps_used = set()

    for pkgname,info in pkgdata.items():
        part2deplist = info['deps']
        check_name_ok(pkgname)
        for part in list(part2deplist.keys()):
            check_name_ok(part)
            assert part.isidentifier()
            part2deplist[part] = sorted(set(e for e in part2deplist[part] if e))
            for dep in part2deplist[part]:
                all_deps_used.add(dep)
                check_name_ok(dep)
                if ( dep not in rawdata['dependencies']
                     and dep not in pkgdata ):
                    raise SystemExit(f'{pkgname}.{part} uses dependency'
                                     f' "{dep}" which is not defined in'
                                     ' dependencies: section and which is not'
                                     ' itself a package name.')

    #Then dependencies:
    for depname,info in rawdata['dependencies'].items():
        check_name_ok(depname)
        if depname not in all_deps_used:
            raise SystemExit('dependencies: section lists '
                             f'unused dependency "{depname}"')
        for k in ['import_check','pypi_pkgname','conda_pkgname']:
            if k not in info:
                info[k] = depname
        if 'min_version' in info and 'max_version' in info:
            raise SystemExit('dependencies: both min and max versions '
                             f'specified for "{depname}" (this is currently'
                             ' not supported by our infrastructure)')

    #add derived info:
    for depname,info in rawdata['dependencies'].items():
        info['py_spec'] = create_py_spec(info)
        info['conda_spec'] = create_conda_spec(info)

    return '',rawdata

def _produce_expected_deplists_for_env_files(db):
    out = {}
    def add(name,deplist):
        assert name not in out
        out[name] = set(deplist)
    db_ncmetapkgdeps = db['packages']['ncrystal-metapkg']['deps']
    for part,deplist in sorted(db_ncmetapkgdeps.items()):
        if part!='BASE':
            add( part, deplist )

    add( 'all_and_devel',
         db_ncmetapkgdeps['all'] + db_ncmetapkgdeps['devel'] )
    return out

def produce_expected_requirements_txt_files(db):
    def content( deplist ):
        lines = []
        for dep in set( deplist ):
            depinfo = db['dependencies'][dep]#assume ncrystal-metapkg has no
                                             #ncrystal-pkg debs outside of BASE
                                             #deps
            lines.append( depinfo['py_spec'] )
        return '\n'.join(sorted(lines)+[''])
    name2deplist = _produce_expected_deplists_for_env_files(db)
    return dict( (f'requirements_{n}.txt', content(d) )
                 for n,d in name2deplist.items() )

_conda_yml_init = """name: ncrystal_<<NAME>>
channels:
  - nodefaults
  - conda-forge
dependencies:
"""

def produce_expected_conda_yml_files(db):
    #For conda we always add whatever is needed for a development environment,
    #since people might wish to pip install plugins, or whatnot:
    deplist_extra = set(['python','pip','c-compiler',
                         'cxx-compiler','cmake','make'])
    def content( name, deplist ):
        lines = sorted(deplist_extra)[:]
        for dep in sorted(set( deplist )-deplist_extra ):
            depinfo = db['dependencies'][dep]#assume ncrystal-metapkg has no
                                             #ncrystal-pkg debs outside of BASE
                                             #deps
            lines.append(depinfo['conda_spec'])
        res = _conda_yml_init.replace('<<NAME>>',name)
        for line in lines:
            res += '  - %s\n'%line
        return res

    name2deplist = _produce_expected_deplists_for_env_files(db)
    return dict( (f'conda_{n}.yml', content(n,d) )
                 for n,d in name2deplist.items() )
