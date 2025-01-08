
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

def _is_exceptional( f, incstatement ):
    if incstatement == 'NCrystal/ncapi.h':
        return True
    from .dirs import coreroot
    from .util import path_is_relative_to
    if not path_is_relative_to( f, coreroot ):
        return False
    frel = str(f.relative_to( coreroot )).replace('\\','/')
    #if ( frel == 'include/NCrystal/plugins/NCPluginBoilerplate.hh'
    #     and incstatement == 'NCrystal/NCrystal.hh' ):
    #    return True
    if ( frel == 'src/utils/NCCFileUtils.cc'
         and incstatement == 'NCCFileUtils.h' ):
        return True

_include_extractor = [None]
def get_include_staments_from_file( path, *,
                                    ignore_exceptional = True ):
    if _include_extractor[0] is None:
        import re
        pattern = b'^\\s*#\\s*include\\s*"\\s*(([a-zA-Z0-9_/.]+))\\s*"'
        _include_extractor[0] = re.compile(pattern).match

    extractor = _include_extractor[0]

    import platform
    if platform.system() != 'Windows':
        #For efficiency, initial dig through file using the grep command:
        import subprocess
        rv = subprocess.run( ['grep','.*#.*include.*"..*"',
                              str(path.absolute())],
                             capture_output = True )
        #grep exit code of 1 simply indicates no hits
        if rv.returncode not in (0,1) or rv.stderr:
            raise RuntimeError('grep command failed')
        if rv.returncode == 1:
            return set()
        lines = rv.stdout.splitlines()
    else:
        #No grep on Windows, use slower fall back:
        lines = []
        for line in path.read_text('utf8').splitlines():
            if 'include' in line:
                lines.append(line.encode('utf8'))
    res = []
    for line in lines:
        v = extractor(line)
        if v:
            v = v.groups()[1].decode('utf8')
            if not ( ignore_exceptional and _is_exceptional( path, v ) ):
                res.append( v )
    return set(res)

def get_nccomp_include_statements( f, *, ignore_list = None ):
    #Iterate over (incstatement,compname_of_inc_statement)
    incs = get_include_staments_from_file( f )
    res = set()
    ignore_list
    for i in incs:
        comp = None
        if i.startswith('NCrystal/internal/'):
            comp = i.split('/')[2]
        elif i.startswith('NCrystal/'):
            comp = i.split('/')[1]
        if comp and ( not ignore_list or comp not in ignore_list ):
            res.add( ( i, comp ) )
    return res
