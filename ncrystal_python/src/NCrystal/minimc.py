
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

"""

Access the NCrystal "MiniMC" framework.

"""

__all__ = [
    'minimc_run',
    'minimc_decode_scenario',
    'minimc_run_scenario',
    'available_tallies',
    'MMCResults',
    'MMCTallyView',
]

#Wrapper objects for result dictionaries:

from ._mmc_utils import MMCResults, MMCTallyView

#Various API entry points:

def minimc_run( cfgstr, *, geomcfg, srccfg, enginecfg,
                callback = None, callback_options = None ):
    """Run the NCrystal MiniMC with the provided cfg-strings for material,
    geometry, source and engine. Returns result wrapped in an MMCResults object.

    Optionally, a callback function can be provided by expert users in order to
    access the full set of all tallied neutrons.
    """
    res = _minimc_raw( cfgstr = cfgstr, geomcfg = geomcfg,
                       srccfg = srccfg, enginecfg = enginecfg,
                       unpack = True, callback = callback,
                       callback_options = callback_options )
    return MMCResults(res)

def minimc_decode_scenario( cfgstr, scenario ):
    """Decodes geomcfg, srccfg, and enginecfg based on the provided material
    cfgstr and scenario string.

    FIXME point at documentation of scenario cfgs.

    The decoded values are returned in a dictionary, along with a normalised
    (FIXME check) version of the material cfgstr.

    """
    from ._common import json_query_cpplayer
    assert isinstance(cfgstr,str), "cfgstr parameter must be a string"
    assert isinstance(scenario,str), "scenario parameter must be a string"
    return json_query_cpplayer(['mmc','scenario',cfgstr, scenario] )

def minimc_run_scenario( cfgstr, scenario, *,
                         extra_engineopts = None,
                         callback = None,
                         callback_options = None ):
    """Convenience function which combines minimc_decode_scenario and minimc_run.

    First the minimc_decode_scenario(..) function is used to obtain geomcfg,
    srccfg, and enginecfg based on the provided material cfgstr and scenario
    string.

    Then the minimc_run(..) function is invoked with these parameters, and the
    result of that function call is returned.

    If extra_engineopts is not None, it must be a string which will be appended
    to the automatically generated engineopts.

    Any provided callback function or callback_options are simply passed along
    to the minimc_run call.

    """
    dec = minimc_decode_scenario( cfgstr, scenario )
    ec = dec['enginecfg']
    if extra_engineopts:
        ec += ';'
        ec += extra_engineopts
    return minimc_run(
        cfgstr = dec['cfgstr'],
        geomcfg = dec['geomcfg'],
        srccfg = dec['srccfg'],
        enginecfg = ec,
        callback = callback,
        callback_options = callback_options,
    )

def _minimc_raw( cfgstr, *,
                 geomcfg, srccfg, enginecfg,
                 unpack=False,
                 callback = None, callback_options = None ):
    """Raw invocation of the MiniMC engine via a low level query, accepting 4
    configuration strings and returning a JSON string with the results.

    Set unpack=True to decode the resulting JSON data and replace any Hist1D
    data inside it with Hist1D objects. Set unpack='json' to simply decode the
    JSON data.

    """
    assert isinstance(cfgstr, str), "cfgstr parameter must be a string"
    assert isinstance(geomcfg, str), "geomcfg parameter must be a string"
    assert isinstance(srccfg, str), "srccfg parameter must be a string"
    assert isinstance(enginecfg, str), "enginecfg parameter must be a string"
    if not callback_options is None:
        assert isinstance(callback_options, str), ("callback_options parameter"
                                                   " must be a string")

    query = ['mmc','run', cfgstr, geomcfg, srccfg, enginecfg]
    if callback:
        from ._chooks import _get_raw_cfcts
        _rawfct = _get_raw_cfcts()
        res = _rawfct['flexmmcrun']( query, callback, callback_options )
    else:
        from ._common import json_query_cpplayer
        res = json_query_cpplayer( query, unpack = False )

    if unpack:
        import json
        res = json.loads(res)
        if unpack != 'json':
            from .hist import Hist1D
            res = Hist1D.objectify_data(res)
    return res

#fixme: more documentation functions like the following (perhaps name all
#doc_...?):
_cache_availtallies=[None]
def available_tallies():
    """Returns dictionary with information about available MiniMC talies that
    can be used in the enginecfg string like enginecfg="tally=mu,de").
    """
    if _cache_availtallies[0] is None:
        from ._common import json_query_cpplayer
        from types import MappingProxyType#to make read-only
        res = MappingProxyType(
            dict( (k,tuple(v)) for k,v in
                  sorted(json_query_cpplayer(['mmc','tallylist']).items()) )
        )
        _cache_availtallies[0] = res
    return _cache_availtallies[0]
