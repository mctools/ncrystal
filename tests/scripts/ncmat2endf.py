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

# fixme: ase is only here temporarily to allow scipy usage
# NEEDS: numpy ase endf-parserpy

import NCTestUtils.enable_fpe # noqa F401
from NCrystalDev.ncmat2endf import ncmat2endf, EndfParameters

def test( cfg, name, endf_parameters, **kwargs ):
    import pprint
    print()
    print()
    print()
    print('-'*80)
    if 'force_save' not in kwargs:
        kwargs['force_save'] = True
    kwargs['ncmat_fn']=cfg
    kwargs['name']=name
    kwargs['endf_parameters']=endf_parameters
    pprint.pprint(kwargs)
    res=ncmat2endf(**kwargs)
    print(res)
    # TODO: replace this with a proper test
    for endf_fn, frac in res:
        print(f"Created file {endf_fn} with fraction {frac}")
        with open(endf_fn) as f:
            lines = [next(f) for _ in range(100)]
        print("".join(lines))

endf_defaults = EndfParameters()
test('Al_sg225.ncmat', 'Al', endf_defaults, temperatures=[350])
test('Polyethylene_CH2.ncmat', 'CH2', endf_defaults, temperatures=[293.6, 350], mat_numbers={"C":37, "H": 38})
