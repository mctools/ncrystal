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

import NCrystalDev._clientry as ncclientry
import shlex

def ncrystalcmd(*args):
    argv = ['/some/where/ncrystal']+[str(e) for e in args]
    print()
    print('='*80)
    print("==> Invoking: %s"%(shlex.join(argv)))
    print('='*80)
    errmsg=None
    try:
        ncclientry.main(argv)
    except SystemExit as e:
        se=str(e)
        if se not in ('','0'):
            errmsg = se
    print('='*80)
    if errmsg is None:
        print("==> Ended with no error")
    else:
        print("==> Ended with SystemExit(%s)"%se)
    print('='*80)

def main():
    ncrystalcmd()
    ncrystalcmd('-l')
    ncrystalcmd('--list')
    ncrystalcmd('-h')
    ncrystalcmd('--help')
    ncrystalcmd('nctool','--cfg','Al_sg225.ncmat ;temp=200K')

if __name__ == '__main__':
    main()
