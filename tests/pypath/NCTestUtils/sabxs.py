
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

# adaptation of .xs utils for inelastic sab
import NCrystalDev.core as nccore
from NCrystalDev._numpy import _np_geomspace, _np_linspace
from NCrystalDev.constants import constant_boltzmann, wl2ekin


def run( testgroup ):
    from .xs import XSMonitor
    def testlist_filtered():
        yield from test_list_gen( testgroup )
    mon = XSMonitor( refdatadir = f'sabxs_{testgroup}',
                     matloadfct = _load_fct,
                     egridgenfct = _egrid_fct,
                     testlistgenfct = testlist_filtered )
    mon.run()

_test_focus = ( 'Al_sg225.ncmat',
                'CaH2_sg62_CalciumHydride.ncmat',
                )
#                'Li2O_sg225_LithiumOxide.ncmat' )



def test_list_gen( testgroup ):
    assert testgroup in ('A','B','C','D')

    if testgroup=='B':
        yield from [
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=0;knllux=0',
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=0;knllux=1',
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=0;knllux=2',
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=0;knllux=3',
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=0;knllux=4',

            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=1;knllux=0',
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=1;knllux=1',
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=1;knllux=2',
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=1;knllux=3',
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=1;knllux=4',

            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2;knllux=0',
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2;knllux=1',
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2;knllux=2',
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2;knllux=3',
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2;knllux=4',

            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=3;knllux=0',
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=3;knllux=1',
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=3;knllux=2',
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=3;knllux=3',
            'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=3;knllux=4',
            ]

    if testgroup=='A':
        import NCTestUtils.enable_testdatapath # noqa F401
        yield from [
            'Li_from_Li2O.ncmat;vdoslux=0;knllux=0',
            'Li_from_Li2O.ncmat;vdoslux=0;knllux=1',
            'Li_from_Li2O.ncmat;vdoslux=0;knllux=2',
            'Li_from_Li2O.ncmat;vdoslux=0;knllux=3',
            'Li_from_Li2O.ncmat;vdoslux=0;knllux=4',

            'O_from_Li2O.ncmat;vdoslux=0;knllux=0',
            'O_from_Li2O.ncmat;vdoslux=0;knllux=1',
            'O_from_Li2O.ncmat;vdoslux=0;knllux=2',
            'O_from_Li2O.ncmat;vdoslux=0;knllux=3',
            'O_from_Li2O.ncmat;vdoslux=0;knllux=4',

            'Li_from_Li2O.ncmat;vdoslux=1;knllux=0',
            'Li_from_Li2O.ncmat;vdoslux=1;knllux=1',
            'Li_from_Li2O.ncmat;vdoslux=1;knllux=2',
            'Li_from_Li2O.ncmat;vdoslux=1;knllux=3',
            'Li_from_Li2O.ncmat;vdoslux=1;knllux=4',

            'O_from_Li2O.ncmat;vdoslux=1;knllux=0',
            'O_from_Li2O.ncmat;vdoslux=1;knllux=1',
            'O_from_Li2O.ncmat;vdoslux=1;knllux=2',
            'O_from_Li2O.ncmat;vdoslux=1;knllux=3',
            'O_from_Li2O.ncmat;vdoslux=1;knllux=4',

            'Li_from_Li2O.ncmat;vdoslux=2;knllux=0',
            'Li_from_Li2O.ncmat;vdoslux=2;knllux=1',
            'Li_from_Li2O.ncmat;vdoslux=2;knllux=2',
            'Li_from_Li2O.ncmat;vdoslux=2;knllux=3',
            'Li_from_Li2O.ncmat;vdoslux=2;knllux=4',

            'O_from_Li2O.ncmat;vdoslux=2;knllux=0',
            'O_from_Li2O.ncmat;vdoslux=2;knllux=1',
            'O_from_Li2O.ncmat;vdoslux=2;knllux=2',
            'O_from_Li2O.ncmat;vdoslux=2;knllux=3',
            'O_from_Li2O.ncmat;vdoslux=2;knllux=4',

            'Li_from_Li2O.ncmat;vdoslux=3;knllux=4',
            'Li_from_Li2O.ncmat;vdoslux=4;knllux=4',
            'O_from_Li2O.ncmat;vdoslux=3;knllux=4',
            'O_from_Li2O.ncmat;vdoslux=4;knllux=4',
        ]



    #FIXME:
    return


    #group A: all files in _test_focus with many configs.
    #group B: anything not in A + files starting with A..H
    #group C: anything not in A + files starting with K..N
    #group D: anything not in A + files starting with P..Z + solid::'s

    #Define list of cfgstrs to test (apart from common factors to be applied in
    #load, to keep filenames shorter)


    from NCrystalDev.datasrc import browseFiles
    is_A = testgroup == 'A'
    take_solids = ( testgroup == 'D' )
    if testgroup == 'B':
        letter_low, letter_up = 'A', 'H'
    elif testgroup == 'C':
        letter_low, letter_up = 'K', 'N'
    elif testgroup == 'D':
        letter_low, letter_up = 'P', 'Z'

    thinning_factor = 7 if is_A else 79
    vdoslux_vals = (0,1,2,3,4)
    vdoslux_vals = (4,)#0,1,2,3,4)#FIXME JUST A TEST
    knllux_vals = (0,1,2,3,4,5)

    i = 0
    for f in browseFiles():
        if f.factName == 'solid':
            if not take_solids:
                continue
        elif f.factName != 'stdlib':
            continue

        is_focus = f.factName=='stdlib' and f.name in _test_focus
        if is_focus != is_A:
            continue

        if ( not is_A
             and f.factName == 'stdlib'
             and not ( letter_low <= f.name[0] <= letter_up ) ):
                continue

        if f.factName == 'solid':
            vdoslux_vals_used = tuple( e for e in vdoslux_vals if e>=3 )
        else:
            vdoslux_vals_used = vdoslux_vals

        temp_vals = (10,None,1000)
        if '::Liquid' in f.fullKey:
            #workaround pre-generated liquid kernels
            temp_vals = (None,)
        if f.fullKey in ('stdlib::Polylactide_C3H4O2.ncmat',
                         'stdlib::AcrylicGlass_C5O2H8.ncmat'):
            #workaround some materials that can not reach 1000K
            temp_vals = (10,None,)

        if f.factName=='stdlib':
            testinfo = nccore.createInfo(f.fullKey+';vdoslux=0;knllux=0')
            if not any( hasattr(di,'loadKernel') for di in testinfo.dyninfos):
                continue

        for t in temp_vals:
            for vdoslux in vdoslux_vals_used:
                for knllux in knllux_vals:
                    c = f.fullKey
                    if c == 'stdlib::Li2O_sg225_LithiumOxide.ncmat':
                        if t!=10:
                            continue
                        import NCTestUtils.enable_testdatapath # noqa F401
                        c = f'Li2O_sg225_LithiumOxide_vdoslux{vdoslux}_temp10K.ncmat'

                    if t is not None:
                        c+=f';temp={t}'
                    c+=f';vdoslux={vdoslux}'
                    c+=f';knllux={knllux}'
                    if thinning_factor > 1:
                        i += 1
                        keep = ( (i-1)%thinning_factor == 0 )
                        if not keep:
                            continue
                    yield c

def _load_fct( cfgstr ):
    if 'knllux=' not in (''.join(cfgstr.split())):
        #NB: add knllux=3 if knllux not specified, to get through the migration
        #period where knllux=-1 is the default.
        cfgstr += ';knllux=3'
    #in any case, only load inelastic:
    cfgstr += ';comp=inelas'
    return nccore.load(cfgstr)

def _egrid_fct( mat ):
    kT = mat.info.getTemperature()*constant_boltzmann
    n = 10
    e = set(_np_linspace( kT*1e-6, kT*1000, n ))
    e |= set(_np_geomspace( kT*1e-6, kT*1000, n ))
    e |= set(_np_geomspace(1e-10, 1e3, n ) )
    e |= {wl2ekin(wl) for wl in _np_linspace(0.1,15.0,n)}
    return e

# stdlib::Pt_sg225.ncmat;vdoslux=0;knllux=4 are inconsistent at the reldiff=0.000964565
# stdlib::Al_sg225.ncmat;temp=1000;vdoslux=1;knllux=4 are inconsistent at the reldiff=0.000106591
# stdlib::Pt_sg225.ncmat;vdoslux=0;knllux=4 are inconsistent at the reldiff=0.000964565 level
# stdlib::Al_sg225.ncmat;temp=1000;vdoslux=1;knllux=4 are inconsistent at the reldiff=0.000106591
