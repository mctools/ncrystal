
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

from ._cliimpl import ( create_ArgumentParser,
                        cli_entry_point,
                        print )
import textwrap

def climod_metadata():
    return dict(
        displaygroup = 'misc',
        displayorder = 15,
        descr = "Send low-level queries for JSON data to NCrystal."
    )

def parseArgs( progname, args, return_parser=False ):
    descr="""Send low level queries for JSON data to the NCrystal C++ library.

    Queries consist of a series of strings and by default the returned JSON data
    will simply be decoded and printed to stdout (via the Python pprint module),
    but this behavour can be overridden by various flags.

    The actual queries and the JSON data returned is unless otherwise noted
    considered an implementation detail of NCrystal, and the format is NOT
    guaranteed to remain stable when new versions of NCrystal are released."""

    helpw = 58
    def hwrap(t):
        return textwrap.fill(' '.join(t.split()),width=helpw)
    parser = create_ArgumentParser(prog = progname,
                                   description=descr)
    parser.add_argument('query', type=str, nargs='+',metavar='STR',
                        help='The actual strings which make up the query.')
    parser.add_argument("-j","--json",action='store_true',
                        help=hwrap("""Output JSON rather than decoded Python
                        dictionary if destination is stdout."""))
    parser.add_argument("--output",'-o',default='',type=str,
                        help=hwrap("""Output file name in which to store
                        generated JSON (defaults to stdout)."""))
    parser.add_argument('--force','-f',action='store_true',
                        help=("Override output file if it already exists."))
    if return_parser:
        return parser
    args = parser.parse_args(args)
    return args

def create_argparser_for_sphinx( progname ):
    return parseArgs(progname,[],return_parser=True)

@cli_entry_point
def main( progname, args ):
    args = parseArgs( progname, args )
    from .misc import evaluate_query
    to_stdout = (not args.output or args.output=='stdout')
    unpack_json = to_stdout and not args.json
    res = evaluate_query( query = args.query,
                          unpack = unpack_json )
    if to_stdout:
        if unpack_json:
            from ._common import ncpprint
            ncpprint(res,do_sort=False)
        else:
            print(res)
    else:
        import pathlib
        p = pathlib.Path( args.output )
        if p.is_file() and not args.force:
            raise SystemExit(f'File already exists (use --force to overwrite): {p}')
        if not p.parent.is_dir():
            raise SystemExit(f'Output directory does not exist: {p.parent}')
        from ._common import write_text
        write_text(p,res)
        print(f"Wrote: {p}")

