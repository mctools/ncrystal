
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
    return 'Analyse version numbers in git tags'

def extract_version( pattern, tagstr ):
    import re
    mo = re.search( pattern.replace('X.Y.Z',
                                    r'([0-9]+\.[0-9]+\.[0-9]+)'),
                    tagstr )
    if mo and len(mo.groups())==1:
        return mo.groups()[0]

def _get_git_version_tag( pattern ):
    from .dirs import reporoot
    import subprocess
    import shutil
    if not shutil.which('git'):
        raise SystemExit('Error: git command not found')
    p = subprocess.run( ['git','tag','--points-at','HEAD'],
                        capture_output = True,
                        check = True,
                        cwd = reporoot )
    if p.returncode != 0 or p.stderr:
        raise SystemExit('Error: git command failed')
    versions = []
    for t in p.stdout.decode().splitlines():
        t = t.strip()
        v = extract_version(pattern,t) if t else None
        if v:
            versions.append(v)
    if not versions:
        raise SystemExit('Error: No tags matching pattern'
                         f' "{pattern}" on HEAD commit')
    if len(versions) > 1:
        raise SystemExit('Error: Multiple tags matching'
                         f' pattern "{pattern}" on HEAD commit')
    return versions[0]

def main( parser ):

    parser.init( """Analyse version numbers in git tags. This is mainly for CI
    usage. If any issues are found, exit with nonzero exitcode. Otherwise, print
    the extracted version.""" )

    parser.add_argument(
        '-t','--tag', metavar='TAG', default='',
        help="""Tag to analyse. If not provided, a "git describe" command will
        be invoked on the repo, requiring the tag to match the --tagpattern and
        to be on the current commit (HEAD)."""
    )
    parser.add_argument(
        '-p','--pattern', metavar='FMT', default='vX.Y.Z',
        help="""Tags pattern. The substring 'X.Y.Z' must be included where
        the three version numbers are expected."""
    )
    parser.add_argument(
        '--fail-if-devel', action='store_true',
        help="""If extracted version number has an odd patch number (the last
        digit), or a patch number >= 80, it is considered a release for internal
        development usage only. This is indicated by a non-zero exit code, and
        can be used to prevent automatic deployment scripts from proceeding."""
    )
    parser.add_argument(
        '--file-verify', metavar='FILE',
        help="""Extracts the version encoded in the provided file (a relative
        path to the repo root), and ends with non-zero exit code if it does not
        match that extracted from the tag"""
    )

    args = parser.parse_args()

    if 'X.Y.Z' not in args.pattern:
        parser.error('Pattern must contain the string "X.Y.Z"')

    if args.tag:
        version_str = extract_version( args.pattern, args.tag )
        if not version_str:
            raise SystemExit(f'Tag "{args.tag}" does not'
                             f' match pattern "{args.pattern}"')
    else:
        version_str = _get_git_version_tag(args.pattern)
        assert version_str

    if args.fail_if_devel:
        major,minor,patch = (int(e) for e in version_str.split('.'))
        if patch % 2 == 1 or patch >= 80:
            raise SystemExit(f'Abort: Version {version_str} '
                             'indicates a development version')

    if args.file_verify:
        from .dirs import reporoot
        f = reporoot.joinpath(args.file_verify)
        if not f.is_file():
            raise SystemExit(f'Error: File not found: {f}')
        from . import _check_versions as cv
        if f.name=='pyproject.toml':
            fv = cv.get_toml_version(f)
        elif f.name=='__init__.py':
            fv = cv.get_py_version(f)
        elif f.name=='VERSION':
            fv = cv.get_versionfile_version(f)
        else:
            raise SystemExit(f'Error: Unknown type of version file: {f}')
        if version_str != fv:
            raise SystemExit(f'Error: Tag version ({version_str}) does not'
                             f' match that extracted from file {f} ({fv})')

    print(version_str)
