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

import NCTestUtils.enable_fpe # noqa F401
import NCrystal as NC

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
        print("  density: [%.14g, %.14g, %.14g, .., %.14g, %.14g, %.14g]"%( density[0],
                                                                            density[1],
                                                                            density[2],
                                                                            density[-3],
                                                                            density[-2],
                                                                            density[-1] ) )
def main():
    cfgstrs = [f.fullKey for f in NC.browseFiles(factory='stdlib')]
    for i,f in enumerate(sorted(cfgstrs)):
        validate_cfgstr(f)


main()
