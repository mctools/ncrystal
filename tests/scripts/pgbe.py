#!/usr/bin/env python3

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

# NEEDS: numpy

# Checks https://github.com/mctools/ncrystal/issues/171

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
import numpy as np
wl_x = np.arange(0.0,10.0,0.01)
pc_x = NC.createScatter("Be_sg194.ncmat;dir1=@crys_hkl:0,0,1@lab:0,1,0;dir2=@crys_hkl:1,0,0@lab:1,0,0;mos=1deg;dirtol=180deg;lcaxis=0,0,1")
pc_x.crossSection(NC.wl2ekin(wl_x),direction=(0,0,1))
wl_x = np.linspace(3.95,3.97,10000)
pc_x.crossSection(NC.wl2ekin(wl_x),direction=(0,0,1))
