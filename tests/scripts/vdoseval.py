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

# NEEDS: mpmath numpy

import NCrystalDev as NC
#from NCrystalDev.vdos import analyseVDOS

import mpmath
mp = mpmath.mp
mp.dps = 200
mpf = mp.mpf
import NCTestUtils.enable_fpe # noqa F401

def integrate_vdos_mpmath(emin,emax,densities,
                          include_bins=True,
                          include_parabola=True):
    e0 = mpf(emin)
    e1 = mpf(emax)
    sumtrapez = mpf(0)
    assert len(densities)>=1
    if include_bins and len(densities)>=2:
        assert len(densities)>=2
        de = (e1-e0)/(len(densities)-1)
        half_de = de / mpf(2)
        for i in range( len(densities) - 1 ):
            contrib = half_de * ( mpf(densities[i]) + mpf(densities[i+1]) )
            sumtrapez += contrib
    #then we just want the integral of the parabola going from (0,0) through (e0,densities[0]):
    #d0 = mpf(densities[0])
    #integral( d0*(x/e0)**2 dx,0,e0 ) = (d0/(3*e0**2)) e0^3 = d0*e0/3
    contrib_parabola = ( mpf(densities[0]) * e0 / mpf(3)
                         if include_parabola
                         else mpf(0) )
    return sumtrapez + contrib_parabola

def integrate_vdos_ncrystal( emin,emax,densities ):
    from NCrystalDev.vdos import analyseVDOS
    d = analyseVDOS(emin=emin,
                    emax=emax,
                    density=densities,
                    temperature=100,#dummy
                    atom_mass_amu = 12.0,)#dummy
    print(d)

def extract_from_cfgstr(cfgstr, element = None):
    info = NC.createInfo(cfgstr)
    if element:
        if isinstance(element,int):
            di = info.dyninfos[element]
        else:
            di = info.findDynInfo(element)
    else:
        di = info.dyninfos[0]
    egrid,densities = di.vdosData()
    assert len(egrid)==2
    d = dict( emin = egrid[0],
              emax = egrid[1],
              densities = densities,
              ncvdoseval_integral = di.analyseVDOS()['integral'] )
    def calc_ref_integral(**kwargs):
        return integrate_vdos_mpmath(emin=d['emin'],
                                     emax=d['emax'],
                                     densities=d['densities'],
                                     **kwargs)
    d['ref_integral'] = calc_ref_integral()
    d['ref_integral_bins'] = calc_ref_integral(include_parabola=False)
    d['ref_integral_parabola'] = calc_ref_integral(include_bins=False)
    d['integral_vs_ref_deviation'] = float((d['ncvdoseval_integral']/d['ref_integral'])-mpf(1.0))
    #help(di)
    #raise SystemExit
    d['displayLabel'] = di.atomData.displayLabel()#fixme:
    return d

def validate_cfgstr(cfgstr):
    info = NC.createInfo(cfgstr)
    any_error = False
    for i,di in enumerate(info.dyninfos):
        if not hasattr(di,'analyseVDOS'):
            continue
        d = extract_from_cfgstr( cfgstr, i )
        lbl = d['displayLabel']
        print(f"Validating VDOS integral of {cfgstr} / {lbl}")
        err = d['integral_vs_ref_deviation']
        if err > 1.0e-14:
            print(f'High error in result: {err}')
            any_error = True
    return not any_error

def main():
    #def analyseVDOS(emin,emax,density,temperature,atom_mass_amu):
    #"""Analyse VDOS curve to extract mean-squared-displacements, Debye temperature,
    #effective temperature, gamma0 and integral. Input VDOS must be defined via
    #an array of density values, over an equidistant energy grid over [emin,emax]
    #(in eV). Additionally, it is required that emin>0, and a parabolic trend
    #towards (0,0) will be assumed for energies in [0,emin]. Units are kelvin and
    #eV where appropriate.
    #"""

    #integrate_vdos_mpmath

    #emin = 0.001
    #emax = 0.1
    #density = [ 0.1, 0.4, 0.2, 0.7 ]
    #print( integrate_vdos_mpmath( 0.1, 0.2, [2.0,2.0], include_parabola=True ) )
    #import numpy
    #print( integrate_vdos_ncrystal( 0.1, 0.2, numpy.asarray([2.0,2.0],dtype=float) ) )#fixme analyseVDOS should add numpy.asarray

    cfgstrs = [f.fullKey for f in NC.browseFiles(factory='stdlib')]
    saw_errors = False
    for i,f in enumerate(sorted(cfgstrs)):
        if i%5 == 0:#skip some for speedup
            if not validate_cfgstr(f):
                saw_errors = True
    if saw_errors:
        print("Errors detected!")
        raise SystemExit(1)


    #ref = extract_from_cfgstr('stdlib::KBr_sg225_PotassiumBromide.ncmat','K')



    #d['integral_vs_ref_deviation']


    #info = NC.createInfo('stdlib::KBr_sg225_PotassiumBromide.ncmat')
    #di=info.findDynInfo('K')
    #egrid,density = di.vdosData()
    #di.analyseVDOS()['integral']
    #import pprint
    #pprint.pprint(ref)


main()
