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

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC

def validate_cfgstr(cfgstr):
    info = NC.createInfo(cfgstr)
    for i,di in enumerate(info.dyninfos):
        if not hasattr(di,'analyseVDOS'):
            continue
        lbl=di.atomData.displayLabel()
        print(f"Validating VDOS regularisation of {cfgstr} / {lbl}:")
        (emin,emax),density = di.vdosData()
        print(f"  egrid [meV]: {emin*1000:.14g} .. {emax*1000:.14g} ({len(density)} points)")
        assert len(density)>=7
        vals = (density[0],density[1],density[2],
                density[-3],density[-2],density[-1] )
        #Worst imprecision is in last few density entries, but not THE last. And
        #it is slightly worse for some of the hydrogren VDOS's. To make tests
        #run, and still retain as much as possible robustness, we vary the
        #output precisions accordingly:
        if lbl=='H':
            print("  density: [%.13g, %.13g, %.13g, .., %.9g, %.9g, %.14g]"%vals )
        else:
            print("  density: [%.13g, %.13g, %.13g, .., %.11g, %.11g, %.14g]"%vals )

def main():
    cfgstrs = [f.fullKey for f in NC.browseFiles(factory='stdlib')]
    for i,f in enumerate(sorted(cfgstrs)):
        validate_cfgstr(f)

main()
