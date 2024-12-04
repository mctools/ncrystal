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

import NCTestUtils.enable_fpe
import NCrystalDev as NC

contents="""NCMAT v[[VERSION]]
#NB: No spacegroup since the symmetry of the atoms below is unknown (to me)
@CELL
  lengths 1 1 1
  angles 90 90 90
@ATOMPOSITIONS
  D 0   1/2 1/2
  D 0     0   0
  Al 1/2 1/2   0
  Al 1/2   0 1/2
@DEBYETEMPERATURE
  Al   410.4
  D   300.4
"""

for version in (1,2,3,4,5,6,7,8):
    print(f'Trying element D in NCMAT data with version {version}')
    NC.registerInMemoryFileData( "TestD.ncmat",contents.replace("[[VERSION]]",str(version)))
    try:
        sc=NC.createScatter("TestD.ncmat;dcutoff=0.5;inelas=0")
    except NC.NCBadInput as e:
        print ('   => Result in exception: ',str(e))
    else:
        xs=sc.crossSectionIsotropic(0.025)
        print (f'   => Create scatter with no problems (xs@25meV is {xs:g}barn).')
    print()
