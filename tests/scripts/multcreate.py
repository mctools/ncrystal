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

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC

def investigate_cross_sections(mat,lbl=None):
    import NCrystalDev.misc as ncmisc
    matsrc = ncmisc.MaterialSource(mat)
    if lbl:
        matsrc.set_plotlabel(lbl)

    print(f"Scattering cross sections of {matsrc.plotlabel}")
    from NCrystalDev.misc import detect_scattering_components
    for comp in detect_scattering_components(matsrc):
        loaded=matsrc.load(f'comp={comp}')
        for wl in ( 0.5, 2.0, 5.0):
            xs = loaded.scatter.xsect(wl=wl)
            #Inelastic cross sections vary a bit too much due to FP instability of
            #the VDOS expansion. But hopefully the vdosgn unit test will monitor the
            #issues for that adequately.
            fmtxs = f'{xs:.7g}' if comp=='inelas' else f'{xs:.14g}'
            print(f"  sigma_{comp} ( {wl:g} Aa ) = {fmtxs} barn")

def main():
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

    investigate_cross_sections(data_Au_sg225_ncmat,'Au_sg225.ncmat (inmem)')

if __name__=='__main__':
    main()
