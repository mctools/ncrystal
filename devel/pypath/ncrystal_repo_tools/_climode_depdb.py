
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

def short_description():
    from .depdb import dbfilerelpath
    return f'Extract information about dependencies from {dbfilerelpath}'

def load():
    from .depdb import load_depdb
    return load_depdb()

def gathertestdeps():
    from .dirs import testroot
    needs = set()
    for s in testroot.joinpath('scripts').glob('*.py'):
        with s.open('rt') as fh:
            for line in fh:
                if line.startswith('# NEEDS:'):
                    needs.update(set(line[len('# NEEDS:'):].split()))
    return needs

def fixreqfiles():
    from .depdb import load_depdb
    from ._check_deps import check_env_files
    db = load_depdb()
    check_env_files( db, fix = True )

def main( parser ):
    parser.init( short_description() )

    parser.add_argument(
        '-m','--mode', type = str, choices = ('fixreqfiles',
                                              'testdeps2json',
                                              'json',
                                              'pprint'),
        default = None,
        help="""Mode (default: pprint)."""
    )

    args = parser.parse_args()

    if args.mode == 'fixreqfiles':
        fixreqfiles()
        return

    if args.mode=='testdeps2json':
        import json
        print( json.dumps(sorted(gathertestdeps())) )
        return

    if args.mode=='pprint':
        import pprint
        pprint.pprint( load() )
        return

    if args.mode=='json':
        import json
        print( json.dumps( load() ) )
        return

    parser.error('No mode selected')
