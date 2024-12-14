
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
    return 'Update C++ component dependencies based on include statements'

def main( parser ):
    parser.init( 'Update all ncrystal_core/src/dep.txt files based on'
                 ' include statements actually seen in the package.' )
    parser.add_argument(
        '-n','--dryrun', action = 'store_true',
        help='Show changes but do not modify anything'
    )
    args = parser.parse_args()

    from .extract_includes import get_include_staments_from_file
    from .core_components import load_components
    for c in load_components( init_deps = False ).values():
        deplist = set()
        for f in c.all_file_iter():
            for i in get_include_staments_from_file(f):
                if not i.startswith('NCrystal/'):
                    continue
                p = i.split('/')
                if len(p) == 4 and i.startswith('NCrystal/internal/'):
                    deplist.add( p[2] )
                elif len(p) == 3 and i.startswith('NCrystal/'):
                    deplist.add( p[1] )
        deptxt = '\n'.join(sorted(e for e in deplist if e != c.name)) + '\n'
        if c.depfile.read_text() != deptxt:
            if args.dryrun:
                print("Would update",c.depfile)
            else:
                print("Updating",c.depfile)
                c.depfile.write_text( deptxt )
