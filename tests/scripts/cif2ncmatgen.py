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

# NEEDS: numpy spglib gemmi

# Tests cif2ncmat with synthetic CIF data (valid CIF data in various styles must
# reproduce the original crystal, invalid CIF data must be rejected). A more
# comprehensive version covering all space groups is in long_cif2ncmatgen.py.

import NCTestUtils.cifgentests

if __name__ == '__main__':
    NCTestUtils.cifgentests.run( ( 2, 14, 70, 141, 166, 194, 227 ),
                                 seed = 1234, invalid_for = ( 14, ) )
