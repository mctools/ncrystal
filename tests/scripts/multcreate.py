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

import NCrystal as NC

data_Au_sg225_nxs="""
space_group = 225
lattice_a = 4.07825
lattice_b = 4.07825
lattice_c = 4.07825
lattice_alpha = 90
lattice_beta  = 90
lattice_gamma = 90
debye_temp = 166.8082
[atoms]
add_atom = Au 7.63 0.43 98.65 196.966570114 0 0 0
"""

mat = NC.directMultiCreate( """NCMAT v1
@CELL
    lengths 4.07825 4.07825 4.07825
    angles 90. 90. 90.
@SPACEGROUP
    225
@ATOMPOSITIONS
    Au 0. 0.5 0.5
    Au 0. 0. 0.
    Au 0.5 0.5 0.
    Au 0.5 0. 0.5
@DEBYETEMPERATURE
    Au   167
""" )
print(mat.info.getDensity())

bla=NC.directMultiCreate(data_Au_sg225_nxs,dtype='nxs',doInfo=False,doScatter=False,doAbsorption=False)
print(bla)
bla=NC.directMultiCreate(data_Au_sg225_nxs,'temp=600K;dcutoff=1',dtype='nxs')
bla.info.dump()
#help(bla.scatter)
print(bla.scatter.crossSectionNonOriented(NC.wl2ekin(2.0)))
