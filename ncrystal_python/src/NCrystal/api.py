
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

Meta-module providing the most commonly needed public API functions and classes
from NCrystal in a single module. It can be used as:

import NCrystal.api as NC

Which will for now do the same as "import NCrystal as NC". However, it might be
that we will eventually modify the default behaviour to not include anything
when merely doing "import NCrystal", so the "import NCrystal.api as NC" will be
more stable in the long rum.

"""

#NB: reduce imported symbols here a bit in a future release (possibly by
#wrapping the removed function and placing in obsolete.py);
from .exceptions import * # noqa F403
from .core import * # noqa F403
from .datasrc import * # noqa F403
from .constants import wl2ekin, ekin2wl, ekin2ksq, wl2k, wl2ksq, constant_boltzmann # noqa F401
from .atomdata import atomDB, iterateAtomDB # noqa F401
from .cfgstr import normaliseCfg, decodeCfg, generateCfgStrDoc # noqa F401
from .ncmat import NCMATComposer, formatVectorForNCMAT # noqa F401
from .plugins import hasFactory, browsePlugins # noqa F401
from ._testimpl import * # noqa F403
from .vdos import createVDOSDebye, debyeIsotropicMSD, PhononDOSAnalyser, debyeTempFromIsotropicMSD, analyseVDOS # noqa F401
from .obsolete import * # noqa F403

#Some modules are left out on purpose (due to esoteric usage or non-standard
#dependencies that most users might not need):
#
# from .cifutils import *
# from .misc import *
# from .mcstasutils import *
#
