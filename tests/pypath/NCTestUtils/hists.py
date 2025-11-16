
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

def _init_lib():
    from NCTestUtils.loadlib import Lib
    lib = Lib('hists')
    lib.dump()
    assert hasattr(lib,'nctest_hist1d_book')
    return lib
_lib = _init_lib()

class Hist1Dcpp:

    def __init__(self, nbins, xmin, xmax, *,
                 allow_weights = True, clamp_overflows = False ):
        self.__id = _lib.nctest_hist1d_book( nbins, xmin, xmax,
                                             1 if allow_weights else 0,
                                             0 if clamp_overflows else 1 )

    def fill( self, x, w = None ):
        if w is None:
            _lib.nctest_hist1d_fill( self.__id, float(x) )
        else:
            _lib.nctest_hist1d_fillw( self.__id, float(x), float(w) )

    def toJSON( self ):
        return _lib.nctest_hist1d_tojson( self.__id )

    def toPyHist( self ):
        data = self.toJSON()
        from NCrystalDev._hist import Hist1D
        return Hist1D(data)
