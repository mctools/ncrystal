
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

    from .dirs import reporoot
    from .srciter import all_files_iter
    ruff = shutil.which('ruff')
    if not ruff:
        raise SystemExit('ERROR: ruff command not available')
    #Default in ruff 0.15.20 (used on FreeBSD) but not in 0.16:
    extsel = 'E402,E721,E741,F403,E743'

    #Empty example file, can not carry a noqa comment (path relative to cwd):
    pfi = ( 'lint.per-file-ignores = {"examples/plugin_dataonly/src/'
            'ncrystal_plugin_DummyDataPlugin/__init__.py" = ["N999"]}' )
    rv = subprocess.run(['ruff','check','--extend-select',extsel,
                         '--config',pfi]
                        + list(all_files_iter('py')),
                        check = False, cwd = reporoot )
    if rv.returncode!=0:
        raise SystemExit(1)

if __name__=='__main__':
    main()
