
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

# For checking some file that must be kept completely synchronised.

def get_content( relpath ):
    from .dirs import reporoot
    print(f'    {relpath}')
    return reporoot.joinpath(relpath).read_text()

def check_same( reffile, *otherfiles ):
    assert len(otherfiles) > 0
    print('  Checking for same contents:')
    ref = get_content( reffile )
    for o in otherfiles:
        if get_content(o) != ref:
            print()
            raise SystemExit(f'ERROR: Content of {o} and {reffile} differs')

def main():
    check_same( 'README.md',
                'ncrystal_core/README.md',
                'ncrystal_core/empty_pypkg/README.md',
                'ncrystal_python/README.md',
                'ncrystal_metapkg/README.md',
                'ncrystal_verify/README.md',
               )

    check_same( 'LICENSE',
                'ncrystal_core/LICENSE',
                'ncrystal_core/empty_pypkg/LICENSE',
                'ncrystal_python/LICENSE',
                'ncrystal_metapkg/LICENSE',
                'ncrystal_pypluginmgr/LICENSE',
                'ncrystal_verify/LICENSE',
                'examples/plugin/LICENSE',
                'examples/plugin_dataonly/LICENSE')

    check_same( 'examples/ncrystal_example_cpp.cc',
                'examples/downstream_cmake/main.cc' )

    #NB: For ncrystal-verify we combine the more robust
    #ncrystal_core/app_test/main.cc with the CMake file from the downstream
    #example.
    check_same( 'ncrystal_core/app_test/main.cc',
                'ncrystal_verify/extra/data/downstream_cmake_main.cc',
                'tests/src/app_cpp/main.cc',
               )
    check_same( 'examples/downstream_cmake/CMakeLists.txt',
                'ncrystal_verify/extra/data/downstream_cmake_CMakeLists.txt' )

if __name__=='__main__':
    main()
