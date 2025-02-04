
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


from ._cliimpl import cli_entry_point

def create_argparser_for_sphinx( progname ):
    raise RuntimeError('Do not call create_argparser_for_sphinx'
                       ' for ncrystal-config')

@cli_entry_point
def main( progname, arglist ):
    import subprocess
    import shutil
    cmd = shutil.which('ncrystal-config')
    assert cmd, 'ncrystal-config command not found!'
    rv = subprocess.run( [cmd]+arglist[:] )
    raise SystemExit(rv.returncode)
