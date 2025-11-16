
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2025 NCrystal developers                                   ##
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

import math

class TestRNG:
    #Very simplistic and unscientific RNG, exclusively to be used for unit test
    #reproducibility (we could have used NCrystal's own proper RNG, but we do
    #not currently expose it to python).
    def __init__( self, seed = 1234567891 ):
        self.__seed = seed
        self.__m = 2 ** 32
        self.__a = 1664525
        self.__c = 1013904223
        self.__state = seed
    def rand01( self ):
        self.__state = (self.__a * self.__state + self.__c) % self.__m
        return self.__state / self.__m
    def stdgauss( self ):
        return ( math.sqrt(-2.0 * math.log(self.rand01()))
                 * math.cos(2.0 * math.pi * self.rand01()) )
    def gauss( self, mu=0.0, sigma=1.0 ):
        return mu + sigma * self.stdgauss()
