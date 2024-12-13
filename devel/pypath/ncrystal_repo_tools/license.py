
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

_licenseblurb_data = """
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
""".lstrip()

def licenseblurb_data_with_doubledash_comments():
    return _licenseblurb_data

_cache_hc = [None]
def licenseblurb_data_with_hash_comments():
    if _cache_hc[0] is not None:
        return _cache_hc[0]
    assert _licenseblurb_data.startswith('///////')
    assert _licenseblurb_data.endswith('///////\n')
    out = []
    for line in _licenseblurb_data.splitlines():
        #replace / with # but only in the frame
        if line.count('/')>40:
            line = line.replace('/','#')
        else:
            assert line.startswith('//') and line.endswith('//')
            line = '##' + line[2:-2] + '##'
        out.append(line)
    out = '\n'.join(out)+'\n'
    assert out.startswith('#######')
    assert out.endswith('#######\n')
    _cache_hc[0] = out
    return _cache_hc[0]

_cache_cc = [None]
def licenseblurb_data_with_ansic_comments():
    if _cache_cc[0] is not None:
        return _cache_cc[0]
    assert _licenseblurb_data.startswith('///////')
    assert _licenseblurb_data.endswith('///////\n')
    out = []
    for line in _licenseblurb_data.splitlines():
        #replace / with # but only in the frame
        if line.count('/')>40:
            assert line[0]=='/'
            assert line[-1]=='/'
            line = line.replace('/','*')
            line = '/' + line[1:-1] + '/'
        else:
            assert line.startswith('//') and line.endswith('//')
            line = '/*' + line[2:-2] + '*/'
        out.append(line)
    out = '\n'.join(out)+'\n'
    assert out.startswith('/******')
    assert out.endswith('******/\n')
    _cache_cc[0] = out
    return _cache_cc[0]
