#!/usr/bin/env python3

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

# NEEDS: numpy

def main():
    import pathlib
    import textwrap
    pathlib.Path('fakepypath/endf_parserpy').mkdir(parents=True)
    pathlib.Path('fakepypath/endf_parserpy/__init__.py').write_text(
        textwrap.dedent(
            """
            raise ImportError('I am a bad endf_parserpy module')
            """
        )
    )
    import sys
    sys.path.insert(0,str(pathlib.Path('fakepypath').absolute()))
    from NCrystalDev.cli import run as ncrun
    from NCTestUtils.common import ensure_error
    with ensure_error(ImportError,'Could not import endf_parserpy'):
        ncrun('ncmat2endf','Al_sg225.ncmat;temp=350K')

if __name__ == '__main__':
    main()
