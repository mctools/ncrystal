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
#component breakdowns and cross sections, and the analyser's turn-key
#Teff/msd values. Covers the msd-based elastic physics of solid
#direct-kernel materials: declared-solid structureless kernels gain
#incoherent-approximation elastic from the analyser msd (reserved for
#scatknl by the NCMAT v5 @STATEOFMATTER spec), Unknown-state and
#liquid materials do not, and crystalline files whose only dynamics is
#a scatter kernel still refuse to load without @DEBYETEMPERATURE
#(pending NCMAT v8). Printed values are deliberately low-precision.

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

def v8_gating_probe():
    #The v8 gating of automatic msd extraction: the same fixture
    #content restamped as NCMAT v8 regains the elastic physics that
    #v7 files must not have (for benzene with a warning naming the
    #H component whose msd is refused), and explicit v8 keywords both
    #enable refused components (msd for H) and override auto-detected
    #values (effective_temperature for C):
    btxt = NC.createTextData('benzene_solid_100K_sabsmall.ncmat').rawData
    for name in ( 'benzene_solid_100K_sabsmall.ncmat',
                  'C_from_benzene_solid_20K_sabsmall.ncmat' ):
        txt = NC.createTextData( name ).rawData
        NC.registerInMemoryFileData( 'v8_'+name,
                                     txt.replace('NCMAT v7','NCMAT v8',1) )
        inspect_proc( 'virtual::v8_'+name )
    parts = btxt.replace('NCMAT v7','NCMAT v8',1).split('@DYNINFO')
    for i, p in enumerate(parts):
        if 'element     H' in p:
            parts[i] = p + '  msd 0.9\n'
        if 'element     C' in p:
            parts[i] = p + '  effective_temperature 600\n'
    NC.registerInMemoryFileData( 'v8kw_benzene.ncmat',
                                 '@DYNINFO'.join(parts) )
    inspect_proc( 'virtual::v8kw_benzene.ncmat' )

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
    flat_O = [ '@DYNINFO', '  element     O', '  fraction    1/3',
               '  type        scatknl', '  temperature 10',
               '  alphagrid   0.01 1 2 3 4 5',
               '  betagrid    -5 -3 -1 1 3 5',
               '  sab_scaled  ' + '1.0 '*36 ]
    li_part = [ p for p in '\n'.join(knl).split('@DYNINFO')
                if 'element     Li' in p ]
    assert len(li_part) == 1
    knl_li_only = ( '@DYNINFO' + li_part[0] ).splitlines()
    for version, name, extra, theknl in (
            ('v7','probe_cell_scatknl.ncmat', [], knl),
            ('v7','probe_cell_scatknl_dt.ncmat',
             ['@DEBYETEMPERATURE','  Li 430','  O 430'], knl),
            #v8: loads WITHOUT Debye temperatures, the msd coming from
            #kernel analysis (feeding Bragg + elastic):
            ('v8','probe_cell_scatknl_v8.ncmat', [], knl),
            #v8 with an unanalysable (flat, liquid-like) O kernel:
            #per-atom load-time error naming O:
            ('v8','probe_cell_badO_v8.ncmat', [],
             knl_li_only + flat_O) ):
        NC.registerInMemoryFileData(
            name,
            '\n'.join( [f'NCMAT {version}'] + struct + extra + theknl )+'\n' )
        print(f'==> load of virtual {name}:')
        try:
            inspect_info( f'virtual::{name}' )
            inspect_proc( f'virtual::{name}' )
        except NC.NCException as e:
            print(f'    refused: {type(e).__name__}: {e}')

def stateofmatter_control_probe():
    #The elastic gain for solid scatknl materials keys on the declared
    #state of matter: the same kernel without @STATEOFMATTER (state
    #Unknown) must stay elastic-free:
    txt = NC.createTextData(
        'Li2O_sg225_LithiumOxide_sabsmall_temp10K.ncmat' ).rawData
    lines = txt.splitlines()
    i = lines.index('@STATEOFMATTER')
    NC.registerInMemoryFileData( 'probe_nostate.ncmat',
                                 '\n'.join( lines[:i]+lines[i+2:] )+'\n' )
    inspect_info( 'virtual::probe_nostate.ncmat' )
    inspect_proc( 'virtual::probe_nostate.ncmat' )

def main():
    mats = [ ('benzene_solid_100K_sabsmall.ncmat','C'),
             ('benzene_solid_100K_sabsmall.ncmat','H'),
             ('C_from_benzene_solid_20K_sabsmall.ncmat',''),
             ('stdlib::LiquidWaterH2O_T293.6K.ncmat','H'),
             ('stdlib::LiquidHeavyWaterD2O_T293.6K.ncmat','D'),
             ('Li2O_sg225_LithiumOxide_sabsmall_temp10K.ncmat','Li'),
             ('Li2O_sg225_LithiumOxide_sabsmall_temp10K.ncmat','O') ]
    for cfgstr in dict.fromkeys( c for c,_ in mats ):
        inspect_info( cfgstr )
        inspect_proc( cfgstr )
    for cfgstr, lbl in mats:
        inspect_auto( cfgstr, lbl )
    stateofmatter_control_probe()
    v8_gating_probe()
    crystalline_scatknl_probe()

if __name__ == '__main__':
    main()
