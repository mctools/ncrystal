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

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
import math
NC.setDefaultRandomGenerator(None)
kToDeg = 57.295779513082322865

def performSomeWork(sc):
    #Some dummy usage of the scatter model.
    print("-"*61)
    wavelength=0.1
    for i in range(13):
        #for wavelength in [0.1+i*0.4 for i in range(0, 13)]:
        ekin=NC.wl2ekin(wavelength)
        xsect = sc.crossSectionIsotropic(ekin)
        ekin_final, mu = sc.sampleScatterIsotropic( ekin )
        delta_ekin = ekin_final - ekin
        scat_angle = math.acos(mu)
        print(f"Sigma(lambda={wavelength:g}Aa)={xsect:g}barn. "
              f"Scattering gave theta={scat_angle*kToDeg:g}degree"
              f" and dE={delta_ekin:g}eV")
        wavelength += 0.4

cfgstr = "Al_sg225.ncmat;dcutoff=0.5;inelas=0;incoh_elas=0;temp=25C"

#First try with the on-disk one:
performSomeWork(NC.createScatter(cfgstr))

#Now register similarly named file in memory, but with slightly larger lattice parameters:
NC.registerInMemoryFileData( "Al_sg225.ncmat",
"""NCMAT v2
@CELL
 lengths 5 5 5 #orig: 4.04958 4.04958 4.04958
 angles 90 90 90
@SPACEGROUP
  225
@ATOMPOSITIONS
  Al 0   1/2 1/2
  Al 0     0   0
  Al 1/2 1/2   0
  Al 1/2   0 1/2
@DEBYETEMPERATURE
  Al   410.4
""")

performSomeWork(NC.createScatter(cfgstr))

#Again, but with more usual lattice parameters:
NC.registerInMemoryFileData( "Al_sg225.ncmat",
"""NCMAT v2
@CELL
 lengths 4.04958 4.04958 4.04958
 angles 90 90 90
@SPACEGROUP
  225
@ATOMPOSITIONS
  Al 0   1/2 1/2
  Al 0     0   0
  Al 1/2 1/2   0
  Al 1/2   0 1/2
@DEBYETEMPERATURE
  Al   410.4
""")

performSomeWork(NC.createScatter(cfgstr))

#Now with in-file cfg:

#Again, but with more usual lattice parameters:
NC.registerInMemoryFileData( "Al_sg225.ncmat",
"""NCMAT v2
# NCRYSTALMATCFG[dcutoffup=1.0]
@CELL
 lengths 4.04958 4.04958 4.04958
 angles 90 90 90
@SPACEGROUP
  225
@ATOMPOSITIONS
  Al 0   1/2 1/2
  Al 0     0   0
  Al 1/2 1/2   0
  Al 1/2   0 1/2
@DEBYETEMPERATURE
  Al   410.4
""")

performSomeWork(NC.createScatter(cfgstr))

for fn in ('Al_sg225.ncmat','Ti_sg194.ncmat','something_that_does_not_exist.ncmat'):
    print(f"=======> {fn}:")
    try:
        td=NC.createTextData(fn)
    except NC.NCFileNotFound:
        td=None
    print(td.rawData if td else '<not-found>')
    print("<=======")


#Just wanted somewhere to test that vdosdebye is not possible without Debye
#temperature:
NC.registerInMemoryFileData( "bad.ncmat","""NCMAT v3
@DYNINFO
  element  Xe
  fraction 1
  type     vdosdebye
@DENSITY
  10 kg_per_m3
""")
try:
    NC.createInfo("virtual::bad.ncmat")
except NC.NCBadInput as e:
    print("Loading bad.ncmat gives expected error:",e)
else:
    raise SystemExit("Did not get expected error")
