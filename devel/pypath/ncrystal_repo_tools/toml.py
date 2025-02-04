
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

def can_parse_toml():
    import sys
    if sys.version_info < (3, 11):
        try:
            import tomli # noqa F401
            return True
        except ImportError:
            return False
    return True

def parse_toml(path,prefer_sys_exit = True):
    import sys
    if sys.version_info >= (3, 11):
        import tomllib as mod
    else:
        import tomli as mod
    import pathlib
    with pathlib.Path(path).open("rb") as f:
        try:
            data = mod.load(f)
        except mod.TOMLDecodeError as e:
            if not prefer_sys_exit:
                raise
            print()
            print('TOMLDecodeError:')
            print()
            print(f'  File: {path}')
            print()
            print(f'  Problem: {e}')
            print()
            raise SystemExit(1)

    return data
