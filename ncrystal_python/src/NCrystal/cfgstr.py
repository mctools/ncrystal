
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

"""

Utilities for accessing NCrystal's database of atomic data, with information
about atomic masses, scattering lengths, etc. Also contains a few other related
utilities, like a list of all element names.

"""

def normaliseCfg(cfgstr):
    """Returns normalised version of cfgstr. This is done behind the scenes by
       loading the specified cfg-string into a C++ MatCfg object and then
       re-encoding it as a string.
    """
    from ._chooks import _get_raw_cfcts
    return _get_raw_cfcts()['nc_normcfgstr'](cfgstr)

def decodeCfg(cfgstr,*,asJSONStr=False):
    """Decodes cfg-string and returns as Python data structure (a dictionary). The
       format of this data structure should be mostly self-evident by
       inspection, and is not guaranteed to stay the same across NCrystal
       versions. If asJSONStr=true, the data structure will be returned as a
       JSON-encoded string, instead of a Python dictionary.
    """
    from ._chooks import _get_raw_cfcts
    _js = _get_raw_cfcts()['nc_cfgstr2json'](cfgstr)
    if asJSONStr:
        return _js
    import json
    return json.loads(_js)

def generateCfgStrDoc( mode = "print" ):
    """Generates documentation about the available cfg-str variables. Mode can
    either be 'print' (print detailed explanation to stdout), 'txt_full' (return
    detailed explanations as string), 'txt_short' (return short explanations as
    string), 'json' (return json-encoded string), or 'python' (return python
    data structures).
    """
    modemap={'print':0,'txt_full':0,'txt_short':1,'json':2,'python':2}
    modeint = modemap.get(mode,None)
    if modeint is None:
        from .exceptions import NCBadInput
        raise NCBadInput('mode must be one of %s'%list(sorted(modemap.keys())))
    from ._chooks import _get_raw_cfcts
    _=_get_raw_cfcts()['nc_gencfgdoc'](modeint)
    if mode == 'print':
        print(_)
    elif mode == 'python':
        import json
        return json.loads(_)
    else:
        return _

def decodecfg_vdoslux(cfgstr):
    """Extract vdoslux value from cfgstr."""
    from ._chooks import _get_raw_cfcts,_str2cstr
    return int(_get_raw_cfcts()['ncrystal_decodecfg_vdoslux'](_str2cstr(cfgstr)))
