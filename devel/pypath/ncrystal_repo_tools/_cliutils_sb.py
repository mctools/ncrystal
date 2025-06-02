
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

"""Module implementing the various simplebuild modes."""

def _cfgfilename( is_debug ):
    return 'simplebuild_debug.cfg' if is_debug else 'simplebuild_reldbg.cfg'

def short_description_sb( mode, is_debug ):
    d = dict(
        sbenv = 'Run command in simplebuild environment',
        sb = 'Build and test code with simplebuild',
        sbrun =  'Build code and run command in simplebuild environment',
    )[mode]
    return '%s (using %s, initial --long for more tests).'%(d,_cfgfilename( is_debug ))

def _find_sbcmd( cmdname ):
    import shutil
    cmd = shutil.which(cmdname)
    if not cmd:
        raise SystemExit('ERROR: Please install simple-build-system')
    return cmd

def _invoke( cmdname, args, env, block = False):
    import subprocess
    cmd = _find_sbcmd( cmdname )
    ev = subprocess.run([cmd]+args, env = env,
                        capture_output = block )
    if ev.returncode != 0:
        if block:
            print( ev.stderr.decode() )
        raise SystemExit(ev.returncode)

def mainsb( mode, is_debug, parser ):
    import os
    from .dirs import reporoot
    args = parser.get_raw_args()
    allow_long_tests = False
    if args and args[0] == '--long':
        allow_long_tests = True
        args = args[1:]

    sbcfg = reporoot.joinpath('devel','simplebuild','cfgs',
                              _cfgfilename( is_debug ))
    env = os.environ.copy()
    env['SIMPLEBUILD_CFG'] = str(sbcfg)
    if allow_long_tests:
        env['NCDEVSBL_ALLOW_LONG_TESTS'] = '1'

    if mode == 'sb':
        _invoke( 'unwrapped_simplebuild', args, env )
    elif mode == 'sbenv':
        _invoke( 'sbenv', args, env )
    elif mode == 'sbrun':
        _invoke( 'unwrapped_simplebuild', ['--quiet'], env, block = True )
        _invoke( 'sbenv', args, env )
    raise SystemExit(0)
