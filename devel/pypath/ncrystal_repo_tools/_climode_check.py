
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

def short_description():
    return 'Run static code checks (use before committing)'

def main( parser ):
    parser.init( """Launch all or some of the various source code checks
    (defined in <reporoot>/devel/pypath/ncrystal_repo_tools/_check_*.py). The
    default is to run all of them.""" )

    parser.add_argument(
        '-l','--list', action='store_true',
        help="""List all available checks and exit."""
    )
    parser.add_argument(
        'CHECK', nargs='*',
        help="""Run this check (runs all checks if nothing is explicitly
        requested)."""
    )
    parser.add_argument(
        '-n','--negate', action='store_true',
        help="""Will run all checks EXCEPT those listed."""
    )

    args = parser.parse_args()
    from . import check_runner
    all_checks = check_runner.get_available_checks_list()
    if args.list:
        for t in all_checks:
            print(t)
        return

    for c in args.CHECK:
        if c not in all_checks:
            candidates = [e for e in all_checks if ( c in e ) or ( e in c ) ]
            advice=''
            if len(candidates)==1:
                advice = ' (perhaps you meant "%s"?)'%candidates[0]
            raise SystemExit(f'Unknown check "{c}"{advice}. Run with --list '
                             'to see available checks.')

    check_list = args.CHECK or all_checks
    if args.negate and args.CHECK:
        check_list = [ c for c in all_checks if c not in args.CHECK ]

    done = set()
    for c in check_list:
        if c not in done:
            check_runner.run_check(c)
            done.add(c)
    print()
    print("DONE: Ran %s checks succesfully"%( len(check_list)
                                              if args.CHECK
                                              else '%i (all)'%len(check_list)))
