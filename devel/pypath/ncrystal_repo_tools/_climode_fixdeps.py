
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

def short_description():
    return 'Update C++ component dependencies based on include statements'

def main( parser ):
    parser.init( 'Update all ncrystal_core/src/dep.txt files based on'
                 ' include statements actually seen in the package.' )
    parser.add_argument(
        '-n','--dryrun', action = 'store_true',
        help='Show changes but do not modify anything'
    )
    args = parser.parse_args()
    from .core_components import fix_deps
    fix_deps( dryrun = args.dryrun )
