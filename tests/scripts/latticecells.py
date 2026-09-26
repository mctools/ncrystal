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

# NEEDS: numpy spglib

# Validates the handling of unit cells for many synthetic cells, trigonal space
# groups in the -31m and -3m1 Laue classes (which were once confused), and all
# crystals in the standard data library. See also long_latticecells.py.

import NCTestUtils.enable_fpe  # noqa: F401
import NCTestUtils.latticecells as lc


def main():
    lc.test_synthetic()
    lc.test_stdlib()
    lc.test_spacegroups( ( 149, 150, 151, 152, 153, 154, 157, 156,
                           159, 158, 162, 164, 163, 165 ) )
    lc.test_rhombohedral_axes()

if __name__ == '__main__':
    main()
