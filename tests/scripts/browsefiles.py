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
import NCrystalDev as NC

#Testing presence of data "files", as well as their sorting.

# More robust testing:
NC.removeAllDataSources()
NC.enableStandardDataLibrary()

files_default_order = {}
files_sorted_order = {}

print()
print('========= In default order ===========')
print()
for f in NC.browseFiles():
    print(f.fullKey)
    if f.factName not in files_default_order:
        files_default_order[f.factName] = []
    files_default_order[f.factName].append(f.fullKey)

print()
print('========= Sorted order ===========')
print()
for f in sorted(NC.browseFiles()):
    print(f.fullKey)
    if f.factName not in files_sorted_order:
        files_sorted_order[f.factName] = []
    files_sorted_order[f.factName].append(f.fullKey)

print()
if files_default_order == files_sorted_order:
    print("All good: Consistent sorting order observed!")
else:
    print("ERROR: Inconsistent sorting order observed!!")
    print()
    import pprint
    pprint.pprint(files_default_order)
    print()
    pprint.pprint(files_sorted_order)
    raise SystemExit("ERROR: Inconsistent sorting order observed!!")

