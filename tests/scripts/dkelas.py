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

#Comprehensive reference test of direct-kernel (scatknl) material
#behaviour at the Info and process levels: state of matter, AtomInfo
#msd/Debye-temperature fields, dyninfo types, elastic/inelastic
#component breakdowns and cross sections, the analyser's turn-key
#Teff/msd values, and the (in)ability to load crystalline files whose
#only dynamics is a scatter kernel. This pins down the behaviour that
#the msd-based elastic support for scatknl materials is about to
#change, so that change will show up as a deliberate update of this
#log. Printed values are deliberately low-precision.

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.enable_testdatapath # noqa F401
import NCrystalDev as NC
from NCrystalDev.misc import evaluate_query as q

def fmt(x,prec=4):
    return 'absent' if x is None else f'{x:.{prec}g}'

def inspect_info( cfgstr ):
    print(f'==> Info inspection of {cfgstr}:')
    info = NC.createInfo(cfgstr)
    print(f'    stateOfMatter={info.stateOfMatter().name}'
          f' crystalline={info.isCrystalline()}'
          f' hasatominfo={info.hasAtomInfo()}')
    for di in info.dyninfos:
        print(f'    dyninfo {di.atomData.displayLabel()}:'
              f' {type(di).__name__} fraction={di.fraction:g}')
    for ai in ( info.atominfos if info.hasAtomInfo() else [] ):
        print(f'    atominfo {ai.atomData.displayLabel()}:'
              f' msd={fmt(ai.msd)} dt={fmt(ai.debyeTemperature)}')

def inspect_proc( cfgstr ):
    print(f'==> Process inspection of {cfgstr}:')
    m = NC.load(cfgstr)
    def combreakdown( s ):
        if s['isNull']:
            return 'NULL'
        sp = s.get('specific') or {}
        if 'components' in sp:
            return '+'.join( sorted( c[1]['name'] for c in sp['components'] ) )
        return s['name']
    print(f'    scatter comps: {combreakdown(m.scatter.getSummary())}')
    for comp in ('elas','inelas'):
        mc = NC.load( cfgstr + f';comp={comp}' )
        xs = ' '.join( f'{mc.scatter.crossSectionIsotropic(e):.4g}'
                       for e in (0.001, 0.0253, 1.0) )
        print(f'    comp={comp}: null={mc.scatter.isNull()} xs(3pts)={xs}')

def inspect_auto( cfgstr, lbl ):
    r = q(['sab','analyse',cfgstr,lbl])['sabanalyse']
    a = r['auto']
    print(f'==> analyser auto values for {cfgstr}//{lbl or "mono"}:'
          f' teff={fmt(a["teff"])} msd={fmt(a["msd"],2)}')

def crystalline_scatknl_probe():
    #A crystalline (unit-cell) file whose only dynamics is a scatter
    #kernel: splice stdlib Li2O structure sections onto the pre-expanded
    #10K kernel file, with and without @DEBYETEMPERATURE. NB: the
    #structure sections are nominally room-temperature, fine for a pure
    #wiring probe at 10K.
    struct = []
    take = False
    for line in NC.createTextData(
            'stdlib::Li2O_sg225_LithiumOxide.ncmat' ):
        if line.startswith('@'):
            take = line.startswith( ('@CELL','@SPACEGROUP',
                                     '@ATOMPOSITIONS') )
        if take:
            struct.append(line)
    knl = []
    take = False
    for line in NC.createTextData(
            'Li2O_sg225_LithiumOxide_sabsmall_temp10K.ncmat' ):
        if line.startswith('@DYNINFO'):
            take = True
        if take:
            knl.append(line)
    for name, extra in ( ('probe_cell_scatknl.ncmat', []),
                         ('probe_cell_scatknl_dt.ncmat',
                          ['@DEBYETEMPERATURE','  Li 430','  O 430']) ):
        NC.registerInMemoryFileData(
            name, '\n'.join( ['NCMAT v7'] + struct + extra + knl )+'\n' )
        print(f'==> load of virtual {name}:')
        try:
            inspect_info( f'virtual::{name}' )
            inspect_proc( f'virtual::{name}' )
        except NC.NCException as e:
            print(f'    refused: {type(e).__name__}: {e}')

def main():
    mats = [ ('stdlib::LiquidWaterH2O_T293.6K.ncmat','H'),
             ('stdlib::LiquidHeavyWaterD2O_T293.6K.ncmat','D'),
             ('Li2O_sg225_LithiumOxide_sabsmall_temp10K.ncmat','Li'),
             ('Li2O_sg225_LithiumOxide_sabsmall_temp10K.ncmat','O') ]
    for cfgstr in dict.fromkeys( c for c,_ in mats ):
        inspect_info( cfgstr )
        inspect_proc( cfgstr )
    for cfgstr, lbl in mats:
        inspect_auto( cfgstr, lbl )
    crystalline_scatknl_probe()

if __name__ == '__main__':
    main()
