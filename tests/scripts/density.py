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

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
i1 = NC.createInfo("Ge_sg227.ncmat;temp=20K")
i2 = NC.createInfo("Ge_sg227.ncmat;density=1.5x;temp=20K")
print ( i1.uid )
print ( i2.uid )
assert i1.uid != i2.uid

print( i1._getUnderlyingUniqueID() )
print( i2._getUnderlyingUniqueID() )
assert i1._getUnderlyingUniqueID() == i2._getUnderlyingUniqueID()
assert i1._getUnderlyingUniqueID() == i1.uid

sc1 = NC.createScatter("Al_sg225.ncmat;temp=20K")
sc2 = NC.createScatter("Al_sg225.ncmat;density=1.5x;temp=20K")
print(sc1.uid)
print(sc2.uid)
assert sc1.uid==sc2.uid
sc1.dump()
sc2.dump()

i3 = NC.createInfo("phases<0.3*Ge_sg227.ncmat;density=1.5x&0.7*Ge_sg227.ncmat>;temp=20K")
assert i3.uid != i1.uid
assert i3.uid != i2.uid
assert len(i3.phases)==2
#print( i3.phases[1][1]._getUnderlyingUniqueID(), i2._getUnderlyingUniqueID(), i1._getUnderlyingUniqueID())
assert i3.phases[1][1]._getUnderlyingUniqueID() == i1._getUnderlyingUniqueID()
assert i3.phases[1][1]._getUnderlyingUniqueID() == i2._getUnderlyingUniqueID()
assert i3.phases[0][1]._getUnderlyingUniqueID() == i1._getUnderlyingUniqueID()
assert i3.phases[0][1]._getUnderlyingUniqueID() == i2._getUnderlyingUniqueID()

i4 = NC.createInfo("phases<0.3*Ge_sg227.ncmat;density=1.5x&0.7*Al_sg225.ncmat>;temp=20K")
assert i4.uid != i1.uid
assert i4.uid != i2.uid
assert len(i4.phases)==2
assert i4.phases[1][1]._getUnderlyingUniqueID() != i1._getUnderlyingUniqueID()
assert i4.phases[1][1]._getUnderlyingUniqueID() != i2._getUnderlyingUniqueID()
assert i4.phases[0][1]._getUnderlyingUniqueID() == i1._getUnderlyingUniqueID()
assert i4.phases[0][1]._getUnderlyingUniqueID() == i2._getUnderlyingUniqueID()
print()
i1.dump()
print()
i2.dump()
print()
i3.dump()
print()
i4.dump()
