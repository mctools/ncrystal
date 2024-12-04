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

#test NCrystal.test() + make sure it doesn't change RNG state

import NCTestUtils.enable_fpe
import NCrystalDev as NC

_rngstate1 = [99]
def rng1():
    global _rngstate1
    _rngstate1[0] = (_rngstate1[0]+1)%100
    return 0.01*(_rngstate1[0]+1)
_rngstate2 = [99]
def rng2():
    global _rngstate2
    _rngstate2[0] = (_rngstate2[0]+1)%100
    return 0.01*(_rngstate2[0]+1)

def nc_use_rng(n=1):
    #bkgd=none => PCBragg => 1rng/call:
    global _sc
    _sc=NC.createScatter("Al_sg225.ncmat;dcutoff=1.5;bkgd=none")
    [_sc.sampleScatterIsotropic(NC.wl2ekin(3.5)) for i in range(n)]
    print("  -> NCrystal consumed %i rngs"%n)

def print_state():
    print("RNG STATES: %02i %02i"%(_rngstate1[0],_rngstate2[0]))

NC.setDefaultRandomGenerator(rng1)
####################
print_state()
nc_use_rng()
print_state()
nc_use_rng(5)
print_state()
if True:#used to be: with NC.RandomCtxMgr(rng2):
    NC.setDefaultRandomGenerator(rng2)#ADDED THIS STATEMENT AFTER WE NO LONGER HAVE NC.RandomCtxMgr
    print("(change rng)")
    print_state()
    nc_use_rng(6)
    print_state()
    print("(change rng back)")
    NC.setDefaultRandomGenerator(rng1)#ADDED THIS STATEMENT AFTER WE NO LONGER HAVE NC.RandomCtxMgr
print_state()
nc_use_rng(6)
print_state()
NC.test()
print_state()
nc_use_rng(1)
print_state()
