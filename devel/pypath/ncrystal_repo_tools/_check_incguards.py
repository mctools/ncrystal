
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

# Checks that all NCrystal headers have a correctly named and unique include
# guard. It also checks that NCrystal/NCrystal.hh contains correct include
# statements for all public api headers.

def get_first_two_lines( f ):
    return f.read_text().splitlines()[0:2]

def check_NCrystal_hh( content, incguards ):
    incguards = dict( (ig,include) for (include,ig) in incguards )
    lines = content.splitlines()
    assert lines[-1] == '#endif'
    lines = list( e.split('//',1)[0].strip() for e in lines[2:-1])
    lines = list( e for e in lines if e )
    i = 0
    while i < len(lines):
        assert lines[i].startswith('#ifndef ')
        ig = lines[i][len('#ifndef '):]
        if ig not in incguards:
            raise SystemExit(f'Forbidden include guard "{ig}" in NCrystal.hh')
        assert lines[i+1]=='#  include "%s"'%incguards[ig]
        assert lines[i+2]=='#endif'
        i += 3
        del incguards[ig]

    expected_missing = ( 'NCrystal/NCrystal.hh',#self
                         'NCrystal/NCPluginBoilerplate.hh',
                         'NCrystal/ncrystal.h',
                         'NCrystal/cinterface/ncrystal.h'
                        )

    for ig, include in incguards.items():
        if include in expected_missing:
            continue
        raise SystemExit(f'NCrystal misses correct include for {include}')
    assert len(incguards) == len(expected_missing)

def main():
    from .srciter import all_files_iter
    from .dirs import coreroot
    from .util import path_is_relative_to

    incroot = coreroot.joinpath('include')
    incrootinternal = coreroot.joinpath('include/NCrystal/internal')
    incguards_seen = set()
    incguards_for_nchh = set()
    f_NCrystal_hh = incroot/'NCrystal/NCrystal.hh'
    f_ncrystal_h = incroot/'NCrystal/cinterface/ncrystal.h'
    f_ncrystal_h_redirect = incroot/'NCrystal/ncrystal.h'

    for f in all_files_iter('*.h','*.hh',root=incroot):
        is_c = ( f.suffix == '.h' )
        ext = 'h' if is_c else 'hh'
        assert '.'+ext == f.suffix
        if is_c:
            assert f.name.startswith('nc')
        else:
            assert f.name.startswith('NC')
        if f.stem.lower().startswith('ncrystal'):
            shortname = f.stem
        else:
            shortname = f.stem[2:]#discard leaving NC
        ig = f'NCrystal_{shortname}_{ext}'
        if f.samefile(f_NCrystal_hh):
            ig = 'NCrystal_hh'
        elif f.samefile(f_ncrystal_h_redirect):
            ig = 'ncrystal_redirection_h'
        elif f.samefile(f_ncrystal_h):
            ig = 'ncrystal_h'

        assert ig not in incguards_seen, "duplicate include guard: %s"%repr(ig)
        incguards_seen.add(ig)
        l1, l2 = get_first_two_lines(f)
        assert l1 == f'#ifndef {ig}', f'Unexpected first line in {f}'
        assert l2 == f'#define {ig}'
        if not path_is_relative_to( f, incrootinternal ):
            incguards_for_nchh.add(
                (str(f.relative_to(incroot)).replace('\\','/'),ig)
            )

    check_NCrystal_hh( f_NCrystal_hh.read_text(),
                       incguards_for_nchh )
    print('All OK!')
