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

# NEEDS: numpy

from NCTestUtils.env import ncsetenv
ncsetenv('DEBUG_VDOSREGULARISATION','1')
ncsetenv('DEBUG_PHONON','1')
import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC # noqa E402

def main():
    info = NC.createInfo('stdlib::BaO_sg225_BariumOxide.ncmat')
    di = info.findDynInfo('O')
    (emin,emax),density = di.vdosData()

main()
