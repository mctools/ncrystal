
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

_prefix = [None]
def ncsetenv(k,v):
    """Tests should use this to set NCrystal-environment variables,
    since depending on the build system they might use different prefixes.
    """

    assert not k.startswith('NCRYSTAL')
    import os
    if _prefix[0] is None:
        from .modeinfo import is_simplebuild_mode
        if is_simplebuild_mode():
            #In simplebuild devel mode we run with a namespace of 'dev' and the
            #NCrystal lib is compiled with NCRYSTAL_NAMESPACED_ENVVARS:
            _prefix[0] = 'NCRYSTALDEV_'
        else:
            #In CTests we don't use NCRYSTAL_NAMESPACED_ENVVARS (and normally
            #not even a namespace):
            _prefix[0] = 'NCRYSTAL_'
    k = _prefix[0] + k
    if v is None:
        if k in os.environ:
            del os.environ[k]
    else:
        os.environ[k] = str(v)
