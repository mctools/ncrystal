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

# NEEDS: numpy

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.stabilise_ncpprint # noqa F401
from NCrystalDev._common import ncpprint
from NCrystalDev.misc import evaluate_query as ncquery

def main():
    query = ['sab', 'sglcell', '@40.0', '@158.28456094789', '@160.60010716438572', '@-2.2701714000978403e-14', '@-2.120250108329995e-17', '@0.02833179339754119', '@0.027532855700519738', '@0.027532855700519426', '@0.026748257342144074', '3', '12462']
    r = ncquery( query, huge_arrays=True)
    ncpprint(r)

if __name__ == '__main__':
    main()
