
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

import numpy as np
from NCrystalDev.hist import HistFiller1D

hist1d_nbins = 70
onedim_projections = [('a',lambda a,b:a),
                      ('b',lambda a,b:b),
                      ('amb',lambda a,b:(a-b)),
                      ('apb',lambda a,b:(a+b))]

class Samples:

    def __init__( self, name, hists_template = None ):
        self.name = name
        self.ntries = 0
        self.n = 0
        self.__hists = {}
        self.__hists_template = hists_template
        self.__data = None

    @property
    def alpha(self):
        assert self.__data is not None, "data not added (with keepdata=True)"
        return self.__data[0]

    @property
    def beta(self):
        assert self.__data is not None, "data not added (with keepdata=True)"
        return self.__data[1]

    def clone_empty( self, name ):
        assert self.__hists_template is not None, "fill before clone_empty"
        return Samples( name, self.__hists_template )

    def AR( self ):
        assert self.ntries >= self.n
        return ( self.n / self.ntries ) if self.ntries else None

    def __init_template( self, avals, bvals ):
        def _bi(x):
            if not len(x)>0:
                return (10,0.0,1.0)
            xmin,xmax,nbins = x.min(),x.max(),hist1d_nbins
            dx = max(xmax-xmin,1e-199)/nbins
            return ( nbins, xmin-dx, xmax+dx )
        self.__hists_template = []
        for title, fct in onedim_projections:
            self.__hists_template.append( (title,fct,
                                           HistFiller1D(_bi(fct(avals,bvals)),
                                                        title=title)))
    def add_data( self, avals, bvals, ntries = None, keepdata = False ):
        avals = np.asarray(avals,dtype=float)
        bvals = np.asarray(bvals,dtype=float)
        if keepdata:
            assert self.__data is None
            self.__data = (avals,bvals)
        if self.__hists_template is None:
            self.__init_template(avals, bvals)
        for title, fct, htemplate in self.__hists_template:
            h=htemplate.clone_empty()
            h.fill(fct(avals,bvals))
            if title not in self.__hists:
                self.__hists[title] = h
            else:
                self.__hists[title].add_contents( h )
        assert len(avals)==len(bvals)
        if ntries is not None:
            self.ntries += ntries
        self.n += len(avals)

    def create_hist( self, key ):
        return self.__hists[key].to_hist1d()
