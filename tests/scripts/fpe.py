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

import NCTestUtils.enable_fpe
from NCTestUtils.loadlib import Lib
from multiprocessing import Process, freeze_support
import platform
import sys
import os

lib_misc = Lib('misc')
div = lib_misc.nctest_divide_args

assert div(1.0,2.0) == 0.5
assert div(10.0,2.0) == 5.0

def worker_div( *div_args ):
    #Redirect stderr to stdout, to avoid a test failure due to stderr output. Use
    #dup2+fileno's, to make it work for the compiled layer as well:
    os.dup2(sys.stdout.fileno(),
            sys.stderr.fileno())
    a,b = div_args
    div( a, b )

def div_in_subproc( a, b ):
    print(f"-----------> Division {a}/{b} starting")
    expect_error = (not b)
    if expect_error:
        print("About to trigger FPE: ------------------>",flush=True)
    p = Process(target=worker_div, args=(a,b))
    p.start()
    p.join()
    ok = p.exitcode==0
    assert ok == (not expect_error)
    print(f"-----------> Division {a}/{b} ended as expected")

def main():
    div_in_subproc( 1.0,2.0)
    div_in_subproc( 10.0,2.0 )
    if platform.system().lower() not in ('windows','darwin'):
        div_in_subproc( 1.0, 0.0 )

if __name__ == '__main__':
    freeze_support()
    main()
