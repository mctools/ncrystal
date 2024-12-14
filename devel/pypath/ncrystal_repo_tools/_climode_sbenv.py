
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2024 NCrystal developers                                   ##
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

def short_description():
    return 'Run command in simplebuild environment'

def main( parser ):
    parser.init( '''Run sbenv associated with the simplebuild setup in
    <reporoot>/devel/simplebuild. This is mainly useful to run a command in that
    environment without modifying the current environment.''' )
    #Fixme many more options here, for now hardcoding below
    args = parser.get_raw_args()
    #Special: Do not actually init or use the parser!
    #args = parser.parse_args()

    import shutil
    import subprocess
    import os
    from .dirs import reporoot

    sb = shutil.which('sbenv')
    if not sb:
        raise SystemExit('ERROR: Please install simple-build-system')


    env = os.environ.copy()
    env['SIMPLEBUILD_CFG'] = str(reporoot.joinpath('devel',
                                                   'simplebuild',
                                                   'simplebuild.cfg'))
    ev = subprocess.run([sb]+args, env = env )
    raise SystemExit(ev.returncode)
