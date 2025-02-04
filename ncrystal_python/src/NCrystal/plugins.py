
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

Utilities related to plugins and dynamic factories.

"""

def hasFactory(name):
    """Check if a factory of a given name exists"""
    from ._chooks import _get_raw_cfcts,_str2cstr
    return bool(_get_raw_cfcts()['ncrystal_has_factory'](_str2cstr(name)))

def browsePlugins(dump=False):

    """Return list of plugins [(pluginname,filename,plugintype),...].

    If the dump flag is set to True, the list will not be returned. Instead it
    will be printed to stdout.
    """
    from ._chooks import _get_raw_cfcts
    ll=_get_raw_cfcts()['ncrystal_get_pluginlist']()
    if not dump:
        return ll
    from ._common import print
    print('NCrystal has %i plugins loaded.'%len(ll))
    for i in range(len(ll)):
        pluginname, filename, plugintype = ll[i]
        print('==> %s (%s%s)'%(pluginname,plugintype,
                             ' from %s'%filename if filename else ''))
