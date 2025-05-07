
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

Python module for using the NCrystal library for thermal neutron transport in
crystals and other materials.

Please find more information about NCrystal at the website:

   https://mctools.github.io/ncrystal/

In particular, a small example using the NCrystal python module can be found at:

   https://github.com/mctools/ncrystal/blob/HEAD/examples/ncrystal_example_py

A substantial effort went into developing NCrystal. If you use it for your work,
we would appreciate it if you would use the following reference in your work:

  X.-X. Cai and T. Kittelmann, NCrystal: A library for thermal neutron
  transport, Computer Physics Communications 246 (2020) 106851,
  https://doi.org/10.1016/j.cpc.2019.07.015

For work benefitting from our inelastic physics, we furthermore request that you
additionally also use the following reference in your work:

  X.-X. Cai, T. Kittelmann, et. al., "Rejection-based sampling of inelastic
  neutron scattering", Journal of Computational Physics 380 (2019) 400-407,
  https://doi.org/10.1016/j.jcp.2018.11.043

For detailed usage conditions and licensing of this open source project, see:

   https://github.com/mctools/ncrystal/blob/master/NOTICE
   https://github.com/mctools/ncrystal/blob/master/LICENSE

"""

#NB: Synchronize meta-data below with fields in setup.py+template_setup.py.in meta data:
__license__ = "Apache 2.0, http://www.apache.org/licenses/LICENSE-2.0"
__version__ = '4.1.8'
__status__ = "Production"
__author__ = "NCrystal developers (Thomas Kittelmann, Xiao Xiao Cai)"
__copyright__ = "Copyright 2015-2024 %s"%__author__
__maintainer__ = __author__
__email__ = "ncrystal-developers@cern.ch"

import sys as _sys
import os as _os

#Place f-string here to catch python <3.6 in a more obvious way than a syntax error below:
f'NCrystal does not work with Python2 (or Python3 < v3.8)' #noqa F541
_minpyversion=(3,8,0)

pyversion = _sys.version_info[0:3]
if pyversion < _minpyversion:
    raise RuntimeError('Unsupported python version %i.%i.%i detected (needs %i.%i.%i or later).'%(pyversion+_minpyversion))

#NB: The following env var can NOT be namespaced. E.g. it will always be
#NCRYSTAL_SLIMPYINIT and never e.g. NCRYSTAL<namespacehere>_SLIMPYINIT:

if not _os.environ.get('NCRYSTAL_SLIMPYINIT'):
    from .api import * # noqa F403

###################################
#Same as NCRYSTAL_VERSION macro (same as the get_version_num() function from .core):
version_tuple = tuple( int(i) for i in __version__.split('.') )#introduced in release 3.6.0
version_num = sum(int(i)*j for i,j in zip(version_tuple,(1000000,1000,1)))#Introduced in release 0.9.1
