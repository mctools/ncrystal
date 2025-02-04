
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

def simplebuild_bundle_list():
    import pathlib
    devel_sbld_dir = ( pathlib.Path(__file__).resolve().absolute()
                       .parent.parent.parent.parent )
    sbldcfg = devel_sbld_dir.joinpath( 'simplebuild.cfg' )
    if ( not sbldcfg.is_file()
         or not devel_sbld_dir.joinpath( 'sbgenerate.py' ).is_file() ):
        raise SystemExit('simple-build-ncrystaldynamic pypkg only works'
                         ' if installed in editable mode')

    return [ sbldcfg ]
