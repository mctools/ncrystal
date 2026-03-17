
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

See also: https://github.com/mctools/ncrystal/wiki/minimc

"""

__all__ = [
    'run',
    'decode_scenario',
    'decode_cfgstr',
    'gen_doc',
    'tally_info',
    'MMCResults',
    'MMCTallyView',
]

#Wrapper objects for results:

from ._mmc_utils import MMCResults, MMCTallyView

#Various API entry points:

def run( cfgstr, *,
         geomcfg = None, srccfg = None, scenario = None,
         enginecfg = None,
         callback = None, callback_options = None,
         unpack = 'object' ):
    """Run the NCrystal MiniMC with the provided cfg-strings for material,
    geometry, source and engine. Returns result wrapped in an MMCResults object.

    Requires a material cfgstr, and either both of geomcfg and srccfg or a
    scenario string.

    Optionally, it is possible to specify an enginecfg string if the default
    settings are not suitable.

    Optionally, a callback function can be provided by expert users in order to
    access the full set of all tallied neutrons.

    The unpack parameter concerns the format of the returned result, and it must
    be either 'dict', 'json', 'dict_jsoncompat', or 'object'. The default option
    'object' cause the results to be wrapped in an MMCResults class, while the
    other options cause the data to be returned in a JSON string or a (possibly
    JSON compatible) dictionary.

    More details at: https://github.com/mctools/ncrystal/wiki/minimc

    """
    from ._mmc_impl import run
    return run( resclass = MMCResults,
                cfgstr = cfgstr, geomcfg = geomcfg,
                srccfg = srccfg, scenario = scenario,
                enginecfg = enginecfg, callback = callback,
                callback_options = callback_options,
                unpack = unpack )

def decode_scenario( cfgstr, scenario ):
    """Decodes the provided material cfgstr and scenario string in order to
    generate a geomcfg and srccfg.

    The decoded values are returned in a dictionary with keys 'geomcfg', and
    'srccfg'.

    For more details about scenario strings, refer to:

      https://github.com/mctools/ncrystal/wiki/minimc_scenario

    Or alternatively invoke the gen_doc function from this module:

      NCrystal.minimc.gen_doc("scenario")

    """
    from .misc import evaluate_query
    assert isinstance(cfgstr,str), "cfgstr parameter must be a string"
    assert isinstance(scenario,str), "scenario parameter must be a string"
    return evaluate_query(['mmc','scenario',cfgstr, scenario] )

def decode_cfgstr( cfgstr, cfgtype ):
    """Decode a MiniMC cfg-string, whose type must be "src", "geom", or
    "engine".

    The decoded values are returned in a dictionary, along with a normalised
    version of the cfg-string itself.
    """
    from .misc import evaluate_query
    if not isinstance(cfgstr,str):
        from .exceptions import NCBadInput
        raise NCBadInput('cfgstr parameter must be a string')
    if not ( isinstance(cfgtype,str) and cfgtype in ('src','geom','engine') ):
        from .exceptions import NCBadInput
        raise NCBadInput('cfgtype parameter must be "src", "geom", or "engine"')
    return evaluate_query(['mmc','inspectcfg',cfgtype, cfgstr])

def gen_doc( subject, mode = None ):
    """Produce reference documentation concerning MiniMC. The single required
    argument specifies the desired documentation subject, and must be one of
    "geom", "src", "engine", or "scenario".

    By default the documentation is simply printed, but the optional mode
    keyword can be used to modify this:

      mode='print' : (default) print the information.
      mode='lines' : return information as list of strings, each
                     representing a single line.
      mode='txt'   : return information as single string.
      mode='dict'  : return information unformatted and in a
                     dictionary (might be empty if not supported).

    Note that documentation is also available at:

      https://github.com/mctools/ncrystal/wiki/minimc

    """
    from ._mmc_doc import gen_doc_impl
    return gen_doc_impl( subject = subject,
                         mode='print' if mode is None else mode )

def tally_info():
    """Returns dictionary with information about available MiniMC talies that
    can be used in the enginecfg string like enginecfg="tally=mu,de").
    """
    from ._mmc_impl import tally_info
    return tally_info()
