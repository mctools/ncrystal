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

#test file with same element in different roles (here differing by MSD value)
import NCTestUtils.enable_fpe
import NCrystal as NC

NC.registerInMemoryFileData( "Al_sg225_multi.ncmat",
"""NCMAT v3
#Fake file, which should have coherent-elastic physics like the standard
#Al_sg225.ncmat. We add a fake and extreeeemely loosely bound Al element in the
#center of the cell.
#
# This tests how the same atom can play different roles in a file.
#
# NCRYSTALMATCFG[inelas=none;incoh_elas=0]
@CELL
 lengths 4.04958 4.04958 4.04958
 angles 90 90 90
@SPACEGROUP
  225
@ATOMPOSITIONS
  X1 0   1/2 1/2
  X1 0     0   0
  X1 1/2 1/2   0
  X1 1/2   0 1/2
  X2 1/2 1/2 1/2 #<--- adding fake loose atom in center... should not contribute to coherent-elastic cross section.
@DEBYETEMPERATURE
  X1 410.4
  X2 10
@ATOMDB
  #Standard Al has (in NCrystal <= v2.0.0 at least): Al 26.981538408u 3.449fm 0.0082b 0.231b
  #
  #Note that if we want similar per-atomic bragg cross sections as in
  #Al_sg225.ncmat, we would have to increase the coherent scat len with
  #sqrt(5/4):
  #  Al 26.981538408u 3.856099227198387fm 0.0082b 0.231b
  #by not having the line above, we get ~identical fsquared factors, but per-atomic
  #cross sections will differ by factor of 4/5.
  X1 is Al
  X2 is Al
""")


NC.createInfo('Al_sg225_multi.ncmat;dcutoff=0.6').dump()
