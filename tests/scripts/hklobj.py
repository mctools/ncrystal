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

import NCTestUtils.enable_fpe
import NCrystalDev as NC

print("test 1)")
mat = NC.load('Al_sg225.ncmat')

for e in mat.info.hklObjects():
    if 0.18 < e.d < 0.5:
        continue
    print(e)

import numpy
print("test 2)")

for idx,e in enumerate(mat.info.hklObjects()):
    assert e.hkl_type == NC.HKLInfoType.SymEqvGroup
    assert e.is_symequiv
    for a in e.h,e.k,e.l:
        assert a.dtype == numpy.int32
        assert isinstance(a,numpy.ndarray)
    assert type(e.d)==float # noqa E721
    assert type(e.f2)==float # noqa E721
    assert type(e.mult)==int # noqa E721
    assert type(e.hkl_label)==tuple # noqa E721
    assert type(e.hkl_label[0])==int # noqa E721
    assert type(e.hkl_label[1])==int # noqa E721
    assert type(e.hkl_label[2])==int # noqa E721
    assert e.f2 == e.fsquared
    assert e.d == e.dspacing
    assert e.multiplicity == e.mult

    if e.d < 0.5:
        if idx%20 != 0:
            continue

    print()
    print(f'{e.hkl_label=}')
    print(f'e.f2={e.f2:.12g}')
    print(f'e.fsquared={e.fsquared:.12g}')
    print(f'e.d={e.d:.14g}')
    print(f'e.dspacing={e.dspacing:.14g}')
    print(f'{e.h=}')
    print(f'{e.k=}')
    print(f'{e.l=}')
    print(f'{e.mult=}')
    print(f'{e.multiplicity=}')

print("test 3)")
info = mat.info
for e in info.hklObjects():
    #help( e );break #<- uncomment for usage info
    print()
    print( e )#<- a quick look
    print( e.hkl_label, e.mult, '%.14g'%e.d, '%.12g'%e.f2 )
    print( e.h, e.k, e.l )#all Miller indices as Numpy arrays.
    #Implement whatever selection logic suits you:
    if ( e.d < 1.5 ):
        break

print("test 4)")
gen = mat.info.hklObjects()
del mat
mat = None

for e in gen:
    print(e)

#some other functions:
def test5():
    print("test 5")
    mat = NC.load('Al_sg225.ncmat')
    wl=[0.5,1.0,1.5,1.8,2.0,5.0]
    print(mat.xsect(wl=wl))
    print(mat.macroscopic_xsect(wl=wl))
    print('%g'%mat.xsect(wl=1.8))
    print('%g'%mat.macroscopic_xsect(wl=1.8))
test5()
