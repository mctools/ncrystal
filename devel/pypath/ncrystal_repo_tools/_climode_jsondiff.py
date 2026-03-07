
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
    return "Diff two json files by diffing their pprint'ed output"

def get_gitfile_contents(commit_hash, file_path):
    import subprocess
    return subprocess.run(
        ['git', 'show', f'{commit_hash}:{file_path}'],
        capture_output=True, check=True, text=True
    ).stdout

def show_through_colordiff_or_print( content, force=False ):
    import subprocess
    import shutil
    import sys
    import os
    cmd = ( shutil.which('colordiff')
            if (force or os.isatty(sys.stdout.fileno()))
            else None )
    if force and not cmd:
        raise SystemExit( 'Error: --color specified by colordiff'
                          ' command not available' )
    if not cmd:
        print(content,end='')
        return
    with subprocess.Popen([cmd],stdin=subprocess.PIPE,text=True) as p:
        p.stdin.write(content)
        p.stdin.close()
        p.wait()

def main( parser ):
    import pathlib
    import json
    import pprint
    import difflib
    import io

    parser.init(short_description()
                +'. If the colordiff command is available,'
                ' it will be used to colorise the output')

    parser.add_argument( 'JSONFILE1',
                         help=("""First json file. If it has the value
                                  git:<tag>:<path> it will load contents of that
                                  file from that git tag.""" ))
    parser.add_argument( 'JSONFILE2',
                         nargs='?',
                         help=("""Second json file. Also supports git syntax
                                  like JSONFILE1. As a special use-case
                                  JSONFILE2 can be left out, meaning to produce
                                  a diff between the git:HEAD version and the
                                  current version of the single file specified.
                         """ ))
    parser.add_argument( '--nocolor', action='store_true',
                         help="""Never pipe output through colordiff.""" )
    parser.add_argument( '--color', action='store_true',
                         help="""Always pipe output through colordiff.""" )

    args = parser.parse_args()
    assert not ( args.color and args.nocolor )
    if args.JSONFILE2 is None:
        args.JSONFILE2 = args.JSONFILE1
        args.JSONFILE1 = f'git:HEAD:{args.JSONFILE1}'
    def load_json_and_pprint(f):
        if f.startswith('git:'):
            p = f.split(':',2)
            assert len(p)==3
            j = get_gitfile_contents(p[1], p[2])
            j = json.loads(j)
        else:
            f = pathlib.Path(f)
            assert f.exists()
            j = json.loads(f.read_text())
        os = io.StringIO()
        pp = pprint.PrettyPrinter(stream=os)
        pp.pprint(j)
        res = os.getvalue()
        os.close()
        return str(f), res.splitlines()
    n1, lines1 = load_json_and_pprint(args.JSONFILE1)
    n2, lines2 = load_json_and_pprint(args.JSONFILE2)

    diff = difflib.unified_diff( lines1, lines2, lineterm='',
                                 fromfile=n1, tofile=n2 )
    out = '\n'.join(diff)+'\n'
    if args.nocolor:
        print(out,end='')
    else:
        show_through_colordiff_or_print(out,force=args.color)
