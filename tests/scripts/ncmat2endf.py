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
from NCrystalDev.ncmat2endf import ncmat2endf
from NCrystalDev.exceptions import NCBadInput
from NCTestUtils.common import print_text_file_with_snipping
import math
import warnings

def reldiff( x, y ):
    if math.isinf(x):
        return ( float('inf') if ( not math.isinf(y) or
               ( ( x>0 ) != ( y>0 ) ) ) else 0.0 )
    return abs(x-y)/(max(1e-300,abs(x)+abs(y)))

def require_flteq( x, y ):
    def okfct( a, b ):
        return bool( reldiff( a, b ) < tol )
    tol = 1e-13
    if hasattr( x, '__len__' ):
        if ( not len(x) == len(y) or
             any( ( not okfct(a,b) ) for a,b in zip(x,y) )):
            raise RuntimeError('numpy flteq failed for arrays '
                              f'x={x} and y={y}!')
    elif not okfct(x,y):
        raise RuntimeError(f'require_flteq( x={x}, y={y} ) failed!')

def test( cfg, ref_teff=None, ref_parsed=None, **kwargs ):
    import pprint
    print()
    print()
    print()
    print('-'*80)
    if 'force_save' not in kwargs:
        kwargs['force_save'] = True
    kwargs['ncmat_cfg']=cfg
    pprint.pprint(kwargs)
    with warnings.catch_warnings():
        # Suppress warnings from ncmat2endf
        warnings.filterwarnings("ignore",category=UserWarning)
        res=ncmat2endf(**kwargs)
    for endf_fn, frac in res:
        print(f"Created file {endf_fn} with fraction {frac}")
        with open(endf_fn) as f:
            text = "".join(f.readlines())
            print_text_file_with_snipping(text,
                                          nstart=100,
                                          nend=25,
                                          prefix='endf>')
        if ref_teff:
            if endf_fn not in ref_teff.keys():
                raise RuntimeError(f'No reference Teff data for {endf_fn}')
            from endf_parserpy import EndfParser
            parser = EndfParser(cache_dir=False)
            endf_dic = parser.parsefile(endf_fn)
            teff = endf_dic[7][4]['teff0_table']['Teff0']
            print(teff, ref_teff[endf_fn])
            require_flteq(teff, ref_teff[endf_fn])
        if ref_parsed:
            if endf_fn not in ref_parsed.keys():
                raise RuntimeError( 'No reference parsed ENDF sections for '
                                   f'{endf_fn}')
            from endf_parserpy import EndfParser, list_parsed_sections
            parser = EndfParser(cache_dir=False)
            endf_dic = parser.parsefile(endf_fn)
            parsed = " ".join([" ".join(str(x) for x in _)
                               for _ in list_parsed_sections(endf_dic)])
            if parsed != ref_parsed[endf_fn]:
                raise RuntimeError( 'ENDF sections {parsed} expected but '
                                   f'sections {ref_parsed[endf_fn]} found')

def test_fail( e, *args, **kwargs ):
    try:
        test(*args,**kwargs)
    except NCBadInput as e:
        print("FAILED (as expected): %s"%e)
        return
    raise SystemExit('Did not fail as expected')

#
# Error handling tests
#

# Oriented materials not supported
test_fail( NCBadInput, 'Ge_sg227.ncmat;dcutoff=0.5;mos=40arcsec;'
           'dir1=@crys_hkl:5,1,1@lab:0,0,1;'
           'dir2=@crys_hkl:0,-1,1@lab:0,1,0')
# Wrong material number assignement
test_fail( NCBadInput, 'Al_sg225.ncmat;vdoslux=1', mat_numbers={"Ge":99})
# Negative temperatures
test_fail( NCBadInput, 'Al_sg225.ncmat;vdoslux=1', temperatures=[-100])
# Repeated temperatures
test_fail( NCBadInput, 'Al_sg225.ncmat;vdoslux=1', temperatures=[293.15])
# No inelastic
test_fail( NCBadInput, 'Al_sg225.ncmat;vdoslux=1;comp=coh_elas')
# Wrong elastic mode
test_fail( NCBadInput, 'Al_sg225.ncmat;vdoslux=1',
           elastic_mode='something wrong')

#
# Library production tests
#

test('Al_sg225.ncmat;vdoslux=1', material_name='Al',
     ref_teff={'tsl_Al.endf':[320.2258, 372.8392]},
     ref_parsed={'tsl_Al.endf':'0 0 1 451 7 2 7 4'},
     temperatures=[350])
test('Polyethylene_CH2.ncmat;vdoslux=1', material_name='CH2',
     ref_teff={'tsl_H_in_CH2.endf':[1208.094],
               'tsl_C_in_CH2.endf':[667.3864]},
     ref_parsed={'tsl_H_in_CH2.endf':'0 0 1 451 7 2 7 4',
                 'tsl_C_in_CH2.endf':'0 0 1 451 7 2 7 4'},
     mat_numbers={"C":37, "H": 38})

