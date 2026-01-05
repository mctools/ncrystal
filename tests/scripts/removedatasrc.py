#!/usr/bin/env python3

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

#If removeAllDataSources does not have an "ensurePluginsLoaded" call internally,
#the code below will load the wrong file from the local directory:

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
import pathlib

pathlib.Path('Al_sg225.ncmat').write_text("INVALID FILE DO NOT LOAD THIS")
NC.removeAllDataSources()
NC.enableStandardDataLibrary()
NC.createInfo("Al_sg225.ncmat;dcutoff=-1")
