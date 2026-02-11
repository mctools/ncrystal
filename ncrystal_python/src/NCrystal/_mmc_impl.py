
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

def results_check_compat_impl( _self, other, threshold, errfct ):
    if _self is other:
        return
    #Allow things like seed, nthreads, and size of statistics to fluctuate:
    volatile=[
        ('engine','cfgstr'),
        ('engine','decoded','seed'),
        ('engine','decoded','nthreads'),
        ('source','cfgstr'),
        ('source','decoded','n'),
    ]
    volatile = set(volatile)
    def cmp(d1,d2,keylist):
        k1=d1.keys()
        if k1!=d2.keys():
            return False
        for k in k1:
            keylist.append(k)
            block = tuple(keylist) in volatile
            if block:
                keylist.pop()
                continue
            v1 = d1[k]
            v2 = d2[k]
            if isinstance(v1,dict) and isinstance(v2,dict):
                if not cmp(v1,v2,keylist):
                    return False
            elif v1 != v2:
                return False
            keylist.pop()
        return True
    keylist=[]
    if not cmp( _self._raw_data()['input'],
                other._raw_data()['input'],
                keylist ):
        k='/'.join(keylist)
        return errfct(f'incompatible input values for "{k}"')

    tallies_s = dict( (t.name, t) for t in _self.tallies )
    tallies_o = dict( (t.name, t) for t in other.tallies )
    tally_names = sorted(tallies_s.keys())

    if ( tally_names != sorted(tallies_o.keys())
         or len(tally_names)!=len(tallies_s)
         or len(tally_names)!=len(tallies_o) ):
        #should already have been caught unless data format changed
        return errfct('incompatible tally list')

    nhists=sum( t._nhists for t in tallies_s.values() )
    actual_threshold = ( 0.0 if threshold is None else threshold ) / nhists

    if actual_threshold <= 0.0:
        return

    for tn in tally_names:
        t_s = tallies_s.get( tn )
        t_o = tallies_o.get( tn )
        if not t_s or not t_o:
            return errfct(f'unexpectedly tally missing: {tn}')
        hists_s = t_s.histograms
        hists_o = t_o.histograms
        hk = hists_s.keys()
        if hists_o.keys() != hk:
            return errfct(f'incompatible histograms available for tally {tn}')
        for hname in hk:
            h_s = hists_s[hname]
            h_o = hists_o[hname]
            if not h_s.check_compat( h_o,
                                     threshold = actual_threshold,
                                     force_norm = True ):
                return (f'incompatible histograms: {tn}/{hk}')
    return
