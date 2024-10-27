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

import NCTestUtils.enable_fpe # noqa F401
import NCrystal as NC

def test1():
    cfgstr='stdlib::Au_sg225.ncmat;temp=600K;dcutoff=1Aa'
    print(f'Debugging "{cfgstr}":')
    from NCrystal.misc import detect_scattering_components
    for comp in detect_scattering_components(cfgstr):
        wl=2.0
        xs=NC.load(cfgstr+f';comp={comp}').scatter.xsect(wl=wl)
        print(f"  sigma_{comp} ( {wl:g} Aa ) = {xs:.14g} barn")

def test2():
    data_Au_sg225_ncmat="""NCMAT v4
@CELL
  cubic 4.07825
@SPACEGROUP
  225
@ATOMPOSITIONS
  Au   0 1/2 1/2
  Au   0   0   0
  Au 1/2 1/2   0
  Au 1/2   0 1/2
@DEBYETEMPERATURE
  Au  167
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

    mat2=NC.directMultiCreate( data_Au_sg225_ncmat,
                               dtype='ncmat',
                               doInfo=False,
                               doScatter=False,
                               doAbsorption=False )
    print(mat2)
    mat2=NC.directMultiCreate(data_Au_sg225_ncmat,'temp=600K;dcutoff=1')
    mat2.info.dump()
    print('%.14g'%mat2.scatter.crossSectionIsotropic(NC.wl2ekin(2.0)))

test1()
test2()
