"""

Meta-module providing the most commonly needed public API functions and classes
from NCrystal in a single module. It can be used as:

import NCrystal.api as NC

Which will for now do the same as "import NCrystal as NC". However, it might be
that we will eventually modify the default behaviour to not include anything
when merely doing "import NCrystal", so the "import NCrystal.api as NC" will be
more stable in the long rum.

"""

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2023 NCrystal developers                                   ##
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

#NB: reduce imported symbols here a bit in a future release (possibly by
#wrapping the removed function and placing in obsolete.py);
from .exceptions import *
from .core import *
from .datasrc import *
from .constants import wl2ekin, ekin2wl, ekin2ksq, wl2k, wl2ksq, constant_boltzmann #TODO: only wl2ekin, ekin2wl
from .atomdata import atomDB, iterateAtomDB
from .cfgstr import normaliseCfg, decodeCfg, generateCfgStrDoc
from .ncmat import NCMATComposer, formatVectorForNCMAT
from .plugins import hasFactory, browsePlugins
from ._testimpl import test
from .vdos import createVDOSDebye, debyeIsotropicMSD, PhononDOSAnalyser, debyeTempFromIsotropicMSD, analyseVDOS
from .obsolete import *

#Some modules are left out on purpose (due to esoteric usage or non-standard
#dependencies that most users might not need):
#
# from .cifutils import *
# from .misc import *
# from .mcstasutils import *
#
