
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

from .exceptions import NCBadInput

_cache_tally_info = [None]
def tally_info():
    if _cache_tally_info[0] is not None:
        return _cache_tally_info[0]
    from .misc import evaluate_query as evalquery
    _cache_tally_info[0] =  evalquery(['mmc','cfgdoc','tally'], readonly = True)
    return _cache_tally_info[0]

def run( *, resclass, unpack,
         cfgstr, geomcfg, srccfg, scenario, enginecfg,
         callback, callback_options ):
    def check( v, sn):
        if not ( v is None or isinstance(v,str) ):
            raise NCBadInput(f'The {sn} parameter must be a string.')
    if cfgstr is None:
            raise NCBadInput('Missing required parameter: cfgstr.')
    enginecfg = '' if enginecfg is None else enginecfg
    check(cfgstr,'cfgstr')
    check(geomcfg,'geomcfg')
    check(srccfg,'srccfg')
    check(scenario,'scenario')
    check(enginecfg,'engiornecfg')
    check(callback_options,'callback_options')
    if unpack not in ('dict', 'json', 'dict_jsoncompat', 'object'):
        raise NCBadInput('Invalid value of unpack (must be "dict",'
                         ' "json", "dict_jsoncompat", or "object"):'
                         f' {repr(unpack)}')
    query = ['mmc','run', cfgstr]#, geomcfg, srccfg, enginecfg]
    n_geomsrc = ( ( 1 if geomcfg is not None else 0 )
                  + ( 1 if srccfg is not None else 0 ) )
    if n_geomsrc == 0 and scenario is None:
        raise NCBadInput('Missing required parameters for geometry and source'
                         '. Please supply either a scenario string,'
                         ' or both of geomcfg + srccfg strings.')
    if n_geomsrc > 0 and scenario is not None:
        raise NCBadInput('Inconsistent parameters. Do not supply geomcfg or'
                         ' srccfg when also supplying a scenario string.')
    if n_geomsrc == 2 and scenario is None:
        query += [ geomcfg, srccfg ]
    elif n_geomsrc == 0 and scenario is not None:
        query += [ scenario ]
    else:
        #Should have been caught above, but just as a safety we throw also here:
        raise NCBadInput('Inconsistent parameters.')
    query += [ enginecfg ]
    if callback:
        from ._chooks import _get_raw_cfcts
        _rawfct = _get_raw_cfcts()
        res = _rawfct['flexmmcrun']( query, callback, callback_options )
    else:
        if callback_options is not None:
            raise NCBadInput('Inconsistent parameters. Do not supply'
                             ' callback_options without a callback function.')
        from .misc import evaluate_query
        res = evaluate_query( query, unpack = False )

    if unpack == 'json':
        return res
    import json
    res = json.loads(res)
    if unpack == 'dict_jsoncompat':
        return res
    from .hist import Hist1D
    res = Hist1D.objectify_data(res)
    if unpack == 'dict':
        return res
    return resclass( res )

def results_check_compat_impl( _self, other, threshold, errfct ):
    if _self is other:
        return
    #Allow things like seed, nthreads, and size of statistics to fluctuate:
    volatile=[
        ('engine','cfgstr'),
        ('engine','decoded','seed'),
        ('engine','decoded','nthreads'),
        ('src','cfgstr'),
        ('src','decoded','n'),
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
