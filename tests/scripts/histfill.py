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

# NEEDS: numpy
import numpy as np
from NCrystalDev.hist import HistFiller1D

def main():
    hist_test = HistFiller1D( 20, -10.0, 10.0, title='test' )
    for i in range(4):
        hist_test.fill( 0.5 )
    for i in range(4):
        hist_test.fill( -0.5, 2.0 )
    hist_test.fill( [-5.0]*4, [2.0]*4 )
    hist_test.fill( (5.0,5.0,5.0,5.0,99999.0), np.ones(5)*2.0 )
    #hist_test.to_hist1d().plot()
    hist_test.to_hist1d().dump()

    print("===> CLONE EMPTY")
    hist_test2 = hist_test.clone_empty(title='test2')
    hist_test2.to_hist1d().dump()
    print("===> CLONE")
    hist_test3 = hist_test.clone(title='test3')
    hist_test3.to_hist1d().dump()
    print("===> RESET CLONE")
    hist_test3.reset()
    hist_test3.to_hist1d().dump()
    print("===> CLONE AGAIN (no title override)")
    hist_test4 = hist_test.clone()
    hist_test4.to_hist1d().dump()
    print("===> ADD CONTENTS OTHER")
    hist_test.add_contents(hist_test4)
    hist_test.to_hist1d().dump()

if __name__ == '__main__':
    main()
