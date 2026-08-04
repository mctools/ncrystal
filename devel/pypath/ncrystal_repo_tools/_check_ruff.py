
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

def main():
    import shutil
    import subprocess

    from .srciter import all_files_iter
    ruff = shutil.which('ruff')
    if not ruff:
        raise SystemExit('ERROR: ruff command not available')
    #TODO: Work on these introduced with ruff 0.16.1:
    ignore=('UP031,C408,C401,RUF059,SIM102,C400,SIM101,C405,SIM118,C402,'
            'N999,RUF015,C403,C419,B018,PLC0206,PERF102,FURB188')
    #For ruff 0.15.20 (needed for FreeBSD):
    ignore += ',E402,E721,E741,F403,E743'

    rv = subprocess.run(['ruff','check','--ignore',ignore]
                        + list(all_files_iter('py')),
                        check = False)
    if rv.returncode!=0:
        raise SystemExit(1)

if __name__=='__main__':
    main()
