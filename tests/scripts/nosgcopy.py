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

import sys

def sec(title):
    sys.stdout.flush()
    print()
    print ('='*100)
    print ('='*100)
    print(f'====  {title.center(88)}  ====')
    print ('='*100)
    print ('='*100)
    print()
    sys.stdout.flush()

def main():
    data = """NCMAT v1
@CELL
    lengths 4.04958 4.04958 4.04958
    angles 90. 90. 90.
@ATOMPOSITIONS
    Al 0. 0.5 0.5
    Al 0. 0. 0.
    Al 0.5 0.5 0.
    Al 0.5 0. 0.5
@DEBYETEMPERATURE
    Al   410.3542
"""

    import NCrystal as NC
    NC.removeAllDataSources()
    NC.registerInMemoryFileData('Al_nosg.ncmat',data)
    NC.registerInMemoryFileData('Al.ncmat',data+'@SPACEGROUP\n    225\n')

    for verbose in (None,0,1,2,3,'sc'):
        for fn in ('Al_nosg.ncmat','Al.ncmat'):
            if verbose is None:
                sec(f'Contents of "{fn}"')
                print(NC.createTextData(fn).rawData)
                continue
            cfg=f'{fn};dcutoff=0.4'
            if verbose=='sc':
                sec(f'Scatter-dump of "{cfg}"')
                NC.createScatter(cfg).dump()
            else:
                sec(f'Info-dump of "{cfg}" with verbosity lvl {verbose}')
                NC.createInfo(cfg).dump(verbose=verbose)

if __name__ == '__main__':
    main()
