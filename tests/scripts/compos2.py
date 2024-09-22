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

import NCrystal as NC

def compos2str( composition ):
    tl = lambda atomdata : ' [TopLevelAtomData]' if atomdata.isTopLevel() else ''
    return '[ %s ]'%(',\n  '.join('%g * %s%s'%(frac,str(atomdata),tl(atomdata)) for frac,atomdata in composition))

def test(cfgstr):
    i=NC.createInfo(cfgstr)
    c=i.getComposition()
    assert isinstance(c,list)
    assert len(c)>0
    print(f'==> Testing cfg "{cfgstr}"\n')
    print( compos2str(c) )
    print('\n')

test("phases<0.7*C_sg227_Diamond.ncmat&0.28*solid::V/6.12gcm3&0.02*Polyethylene_CH2.ncmat;comp=inelas>;temp=380K")
test("phases<1.0*C_sg227_Diamond.ncmat>;temp=380K")
test("C_sg227_Diamond.ncmat")
test("solid::B4C/2.52gcm3/B_is_0.95_B10_0.05_B11")
test("phases<0.5*solid::B4C/2.52gcm3/B_is_0.95_B10_0.05_B11&0.5*gasmix::BF3/2bar>")
