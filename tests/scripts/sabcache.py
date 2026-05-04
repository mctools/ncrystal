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

# Verify that loading the same VDOS in different contexts will ultimately result
# in the same SABScatterHelper being used below the hood.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev.core as nccore
import NCTestUtils.enable_testdatapath # noqa F401
from NCTestUtils.sab import extract_scathelper_uids

def main():
    mat_al300k = nccore.load("stdlib::Al_sg225.ncmat;temp=300K")
    mat_alvdoslux2 = nccore.load("stdlib::Al_sg225.ncmat;temp=200K;vdoslux=2")
    mat_al = nccore.load("stdlib::Al_sg225.ncmat;temp=200K")
    mat_mp = nccore.load("""phases<
                              0.3 * stdlib::Al_sg225.ncmat;strain=0.1
                            & 0.2 * stdlib::Al_sg225.ncmat;strain=-0.1
                            & 0.5 * stdlib::Al_sg225.ncmat
                            >;temp=200K""")
    mat_al2o3 = nccore.load("stdlib::Al2O3_sg167_Corundum.ncmat;temp=200K")
    mat_al2o3_special = nccore.load("Al2O3_Alvdos_from_Al225.ncmat;temp=200K")

    su = extract_scathelper_uids
    uids_al300k = su(mat_al300k)
    uids_alvdoslux2 = su(mat_alvdoslux2)
    uids_al = su(mat_al)
    uids_mp = su(mat_mp)
    uids_al2o3 = su(mat_al2o3)
    uids_al2o3_special = su(mat_al2o3_special)

    #The three first Al have different conditions for SABs:
    assert len(uids_al300k) == 1
    assert len(uids_alvdoslux2) == 1
    assert len(uids_al) == 1
    uids_all_al = list(set(uids_al300k+uids_alvdoslux2+uids_al))
    assert len(set(uids_al300k+uids_alvdoslux2+uids_al)) == 3
    #All the multiphase materials should have the same SAB condition as in the
    #mat_al:
    assert len(uids_mp) == 1
    assert uids_mp == uids_al
    #Here should be two completely new ones:
    assert len(uids_al2o3) == 2
    assert set(uids_all_al) - set(uids_al2o3) == set(uids_all_al)
    #And in this weird material, the Al should match that from mat_al and the O
    #should match the O from mat_al2o3:
    assert len(uids_al2o3_special) == 2
    assert uids_al[0] in uids_al2o3_special
    uids_al2o3_special_O = list(set(uids_al2o3_special)-set(uids_al))
    assert len(uids_al2o3_special_O) == 1
    assert uids_al2o3_special_O[0] in uids_al2o3
    print("All OK")

if __name__ == '__main__':
    main()
