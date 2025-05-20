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
# fixme: split into different tests
# NEEDS: numpy ase endf-parserpy

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.reprint_escaped_warnings # noqa F401
from NCrystalDev.ncmat2endf import ncmat2endf
from NCrystalDev.ncmat2endf import EndfMetaData
from NCrystalDev.exceptions import NCBadInput
from NCTestUtils.common import print_text_file_with_snipping
import NCrystalDev.cli as nc_cli
import shlex
import math

#FIXME: Add tests with:
import NCrystalDev._ncmat2endf_impl as ncmat2endf_impl
ncmat2endf_impl.unit_test_chop_svals[0] = True
#ncmat2endf_impl.is_unit_test[0] = True # just print final dict

def reldiff( x, y ):
    # FIXME: Add these functions to NCTestUtils
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

def test_cli( args ):
    if isinstance(args,str):
        args = shlex.split(args)
    hr=f"============= CLI >>{shlex.join(args)}<< ===================="
    print(hr)
    nc_cli.run('ncmat2endf',*args)
    print('='*len(hr))

def test( cfg, ref_teff=None,
          ref_parsed=None, ref_bragg_edges=None,
          **kwargs ):
    import pprint
    print()
    print()
    print()
    print('-'*80)
    if 'force_save' not in kwargs:
        kwargs['force_save'] = True
    kwargs['ncmat_cfg']=cfg
    pprint.pprint(kwargs)
    res = ncmat2endf(**kwargs)
    for endf_fn, frac in res:
        print(f"Created file {endf_fn} with fraction {frac}")
        with open(endf_fn) as f:
            text = "".join(f.readlines())
            print_text_file_with_snipping(text,
                                          nstart=140,
                                          nend=70,
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
        if ref_bragg_edges:
            if endf_fn not in ref_bragg_edges.keys():
                raise RuntimeError( 'No reference Bragg edges for '
                                   f'{endf_fn}')
            parser = EndfParser(cache_dir=False)
            endf_dic = parser.parsefile(endf_fn)
            ref_Eint, ref_S0 = ref_bragg_edges[endf_fn]
            Eint = tuple(endf_dic[7][2]['S_T0_table']['Eint'][:len(ref_Eint)])
            S0 = tuple(endf_dic[7][2]['S_T0_table']['S'][:len(ref_S0)])
            require_flteq(Eint, ref_Eint)
            require_flteq(S0, ref_S0)

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
metadata = EndfMetaData()
metadata.set_mat_numbers( {"Ge":99} )
test_fail( NCBadInput, 'Al_sg225.ncmat;vdoslux=1', endf_metadata=metadata)
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
# CLI tests
#

test_cli('"stdlib::Al_sg225.ncmat;temp=350K;vdoslux=1"'
         ' --verbosity 5 -m Al -f -e greater')

#
# Library production tests
#

ref_Eint = (0.003741252, 0.004988336, 0.009976672, 0.01371792, 0.01496501,
            0.01995334, 0.0236946, 0.02494168, 0.02993002, 0.03367127,
            0.03990669, 0.04364794, 0.04489502, 0.04988336, 0.05362461,
            0.05487169, 0.05986003, 0.06360128, 0.06484837, 0.0698367)
ref_S0 = (0.005134741, 0.00839204, 0.01258346, 0.01924357, 0.02131949,
          0.02254634, 0.02674178, 0.03073558, 0.03405288, 0.03793785,
          0.03912707, 0.04336451, 0.04591494, 0.04767616, 0.04925869,
          0.05078662, 0.05123031, 0.05363637, 0.05479998, 0.05684046)

test('Al_sg225.ncmat;vdoslux=1', material_name='Al',
     ref_teff={'tsl_Al.endf':[320.2258, 372.8392]},
     ref_parsed={'tsl_Al.endf':'0 0 1 451 7 2 7 4'},
     ref_bragg_edges={'tsl_Al.endf':(ref_Eint, ref_S0)},
     temperatures=[350], elastic_mode='scaled')

metadata = EndfMetaData()
metadata.set_mat_numbers( {"C":37, "H": 38} )
test('Polyethylene_CH2.ncmat;vdoslux=1', material_name='CH2',
     ref_teff={'tsl_H_in_CH2.endf':[1208.094],
               'tsl_C_in_CH2.endf':[667.3864]},
     ref_parsed={'tsl_H_in_CH2.endf':'0 0 1 451 7 2 7 4',
                 'tsl_C_in_CH2.endf':'0 0 1 451 7 2 7 4'},
     endf_metadata=metadata)

