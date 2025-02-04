
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

"""

This __main__.py module is here to ensure that NCrystal commandline script
can be accessed even when not terminal scripts (entry points) are
installed. This is for instance done via:

$> python3 -mNCrystal nctool <args>
$> python3 -mNCrystal ncmat2cpp <args>
$> python3 -mNCrystal ncrystal_ncmat2cpp <args>

Note that the two last commands above run the same command line script.

"""

def _print_usage():
    from ._common import print
    from .cli import cli_tool_list
    import textwrap
    usagestr = textwrap.dedent("""
    Usage: provide name and arguments of NCrystal commandline-tool to run.

    Available tools:

    <<TOOLLIST>>
    Specifying -h or --help displays this message.""").lstrip()
    toolliststr = ''.join( f'       {t}\n' for t
                           in cli_tool_list( canonical_names = False ))
    print(usagestr.replace('<<TOOLLIST>>',toolliststr))

def _prepare():
    import sys
    from ._cliimpl import _resolve_cmd_and_import_climod as _resolve

    args = sys.argv[1:]

    initial_optargs = []
    for a in args:
        if a.startswith('-'):
            initial_optargs.append(a)
        else:
            break
    args = args[len(initial_optargs):]

    is_help_req = False
    opt_unblock = False
    bad_usage = False
    for a in initial_optargs:
        if a in ('-h','--help','/?'):
            is_help_req = True
        elif a == '--unblock':
            opt_unblock = True
        else:
            bad_usage = True

    if is_help_req or bad_usage or not args:
        _print_usage()
        raise SystemExit(0 if is_help_req else 1)

    return _resolve( args[0], args[1:] ), opt_unblock

def main():
    try:
        (climod, argv), opt_unblock = _prepare()
    except RuntimeError as e:
        from ._common import print
        _print_usage()
        print()
        raise SystemExit(f'ERROR: {e}')
    if opt_unblock:
        #Run like in pyapi (in particular letting all exceptions escape):
        from .cli import run
        run( *argv )
    else:
        climod.main( argv )

if __name__ == '__main__':
    main()
