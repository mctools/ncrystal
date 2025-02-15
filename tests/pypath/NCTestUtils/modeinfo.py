
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

_mode_is_sbld = [None]
def is_simplebuild_mode():
    """Whether or not using simplebuild development mode"""
    if _mode_is_sbld[0] is None:
        import pathlib
        _mode_is_sbld[0] = ( pathlib.Path(__file__).parent
                             .joinpath('_is_simplebuild.py').is_file() )
    return _mode_is_sbld[0]

def is_ncrystalverify_mode():
    """Whether or not exported as part of ncrystal-verify package"""
    return False
