#!/usr/bin/env python3

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

# NEEDS: toml

def main():
    from NCTestUtils.dirs import ncrystal_repo_root
    import sys
    import subprocess
    ncdevtool = ncrystal_repo_root.joinpath('devel','bin','ncdevtool')
    rv = subprocess.run( [ sys.executable or 'python3', '-BI',
                           ncdevtool, 'check', '-n','fixme' ] )
    if rv.returncode != 0:
        raise SystemExit("Check failed")
    print("All checks passed")

if __name__ == '__main__':
    main()
