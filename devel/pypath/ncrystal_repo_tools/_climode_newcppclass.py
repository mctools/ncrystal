
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
    return 'Create new C++ file from skeleton class'

def main( parser ):
    parser.init( 'Create a new C++ class from a template.' )

    parser.add_argument(
        '--name', metavar='CamelCasedName',
        required = True,
        help="""Name of new class."""
    )
    parser.add_argument(
        '-c','--comp', metavar='COMPNAME',
        required = True,
        help="""NCrystal component name like "core", "utils", etc. to put class
        .hh/.cc files in."""
    )
    parser.add_argument(
        '-n','--new', choices = ('internal','public'), metavar='TYPE',
        help="""Allow --comp to specify an new component which is not yet
        present. TYPE must be either 'internal' or 'public'."""
    )


    args = parser.parse_args()

    if not args.comp or not args.comp.isidentifier():
        parser.error('Invalid component name')

    from .core_components import load_components
    from .dirs import coreroot, devpymoddir

    for f in coreroot.rglob(f'NC{args.name}.hh'):
        parser.error(f'ERROR: Chosen name conflicts with existing file: {f}')

    name2comp = load_components()

    files_to_create = []
    if args.comp not in name2comp:
        if not args.new:
            parser.error(f'Unknown component "{args.comp}" (supply '
                         '--new TYPE to create).')
        is_internal = args.new == 'internal'
        if is_internal:
            p_hdr = coreroot / f'include/NCrystal/internal/{args.comp}'
        else:
            p_hdr = coreroot / f'include/NCrystal/{args.comp}'
        p_src = coreroot / f'src/{args.comp}'
        files_to_create.append( (p_src.joinpath('dep.txt'), 'core\n') )
    else:
        comp = name2comp[args.name]
        p_hdr = comp.hdrdir
        p_src = comp.srcdir

    skelhh = devpymoddir.joinpath('templates/NCSkeleton.hh')
    skelcc = devpymoddir.joinpath('templates/NCSkeleton.cc')
    txt_hh = skelhh.read_text().replace('Skeleton',args.name)
    txt_cc = ( skelcc.read_text()
               .replace('Skeleton',args.name)
               .replace('"INCLUDEPATH/',
                        '"%s/'%p_hdr.relative_to(coreroot.joinpath('include'))))
    files_to_create.append( ( p_hdr.joinpath(f'NC{args.name}.hh'), txt_hh) )
    files_to_create.append( ( p_src.joinpath(f'NC{args.name}.cc'), txt_cc) )
    for f,txt in files_to_create:
        if f.exists():
            parser.error(f'ERROR: File already exists: {f}')

    for f,txt in files_to_create:
        print(f'Creating {f}')
        f.parent.mkdir(parents=True,exist_ok=True)
        f.write_text(txt)
