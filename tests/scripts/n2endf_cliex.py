#!/usr/bin/env python3

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

# NEEDS: numpy endf-parserpy

# Launch the examples from ncmat2endf --help text.

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.reprint_escaped_warnings # noqa F401
from NCTestUtils.ncmat2endf_utils import test_cli
import NCrystalDev._ncmat2endf_impl as ncmat2endf_impl
from NCrystalDev._cli_ncmat2endf import examples as cli_example_list

ncmat2endf_impl.unit_test_chop_vals[0] = True
ncmat2endf_impl.unit_test_abort_write[0] = 'dump'

def main():
    import copy
    for arglist in copy.deepcopy(cli_example_list):
        assert '.ncmat' in arglist[0]
        arglist[0] += ';vdoslux=0'
        test_cli(*arglist)

if __name__ == '__main__':
    main()
