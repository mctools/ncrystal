
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


def short_description():
    return 'Analyse C++ component dependencies based on include statements'

def main( parser ):
    from .core_components import load_components
    from .extract_includes import get_include_staments_from_file as getinc
    from .dirs import coreroot, coreroot_include

    parser.init( """Analyse C++ code dependencies by looking for include
    statements in daughter components to header files in mother components.""" )

    parser.add_argument(
        '-m','--mother', nargs='+', required = True, metavar='COMPNAME',
        help="""Code component names (e.g. "core", "utils", ...) that you want
        to investigate usage of in daughter components."""
    )
    parser.add_argument(
        '-d','--daughter', nargs='+', required = True, metavar='COMPNAME',
        help="""Code component names (e.g. "core", "utils", ...) in which you
        want to investigate uses of mother components."""
    )
    parser.add_argument(
        '--hdr', action = 'store_true',
        help="""If this is enabled, only include statements in header files are
        listed."""
    )
    parser.add_argument(
        '--nonrecursive', action = 'store_true',
        help="""If this is enabled, includes inside mother header files are
        ignored."""
    )
    args = parser.parse_args()

    name2comp = load_components()
    for e in args.mother + args.daughter:
        if e not in name2comp:
            parser.error(f'Unknown component: {e}')


    if args.hdr:
        def get_files( comp ):
            return comp.hdrfiles
    else:
        def get_files( comp ):
            for f in comp.all_file_iter():
                yield f

    incstatements_prefixes = tuple(
        f'NCrystal/internal/{n}/' if name2comp[n].is_internal else f'NCrystal/{n}/'
        for n in args.mother
    )

    def get_relevant_includes( f ):
        for i in getinc(f):
            if any(i.startswith(e) for e in incstatements_prefixes):
                yield i

    #Expand mother includes (i.e. if a daughter file includes a mother hdrfile
    #which includes another mother hdrfile, that daughter file sees both of
    #these):
    motherincs = {}
    for m in ( [] if args.nonrecursive else args.mother ):
        c = name2comp[m]
        for f in sorted(c.hdrfiles):
            fn = str(f.relative_to(coreroot_include))
            assert fn not in motherincs
            ri = set()
            for i in get_relevant_includes(f):
                ri.add( i )
            motherincs[fn] = ri

    def expand_includes( includes ):
        done = set()
        todo = set( i for i in includes )
        while todo:
            i = todo.pop()
            for i2 in motherincs[i]:
                if i != i2 and i2 not in done and i2 not in todo:
                    todo.add( i2 )
            done.add(i)
        return done

    def get_relevant_includes_expanded( f ):
        return expand_includes( set( get_relevant_includes( f ) ) )

    def get_daughter_includes():
        i2f = {}
        for d in args.daughter:
            for f in get_files(name2comp[d]):
                for i in get_relevant_includes_expanded(f):
                    if i not in i2f:
                        i2f[i] = set()
                    i2f[i].add( str(f.relative_to(coreroot)) )
        for i,fs in sorted( (i,sorted(fs)) for i,fs in i2f.items()):
            yield i,fs

    included = set()
    for i,fs in get_daughter_includes():
        included.add( i )
        print(f'* {i} is included by:')
        print()
        for f in fs:
            print(f'     {f}')
        print()

    print("Header-files not included in these daughters:")
    for m in args.mother:
        c = name2comp[m]
        for f in sorted(c.hdrfiles):
            fn = str(f.relative_to(coreroot_include))
            if fn not in included:
                print('   ',fn)
