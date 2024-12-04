
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

def can_parse_toml():
    import sys
    if sys.version_info < (3, 11):
        try:
            import tomli
            return True
        except ImportError:
            return False
    return True

def parse_toml(path):
    import sys
    if sys.version_info >= (3, 11):
        import tomllib as mod
    else:
        import tomli as mod
    import pathlib
    with pathlib.Path(path).open("rb") as f:
        data = mod.load(f)
    return data
