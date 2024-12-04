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

def otherphasescfgstr(cfgstr):
    return (f"""NCMAT v6
@DYNINFO
  element H
  fraction 1
  type sterile
@DENSITY
  1 g_per_cm3
@OTHERPHASES
  0.001 {cfgstr}#comment
""",)
def bad_otherphasescfgstr(cfgstr):
    return otherphasescfgstr(cfgstr)[0]

#utf8_filename
a_ring_char = b'\xc3\xa5'.decode('utf-8')#Funny danish letter (a with ring on top)
NC.registerInMemoryFileData('utf8_filename.ncmat',NC.createTextData('void.ncmat').rawData)
NC.registerInMemoryFileData(f'utf8_filen{a_ring_char}me.ncmat',NC.createTextData('void.ncmat').rawData)


ll=[
    #UTF8 chars never allowed in ncmat data, even if normally allowed in filenames in cfg strings:
    bad_otherphasescfgstr(f'utf8_filen{a_ring_char}me.ncmat'),
    #no newlines:
    bad_otherphasescfgstr(f'utf8_filen{a_ring_char}me.ncmat ; \n dcutoff = 0.4'),
    bad_otherphasescfgstr('utf8_filename.ncmat ; \n dcutoff = 0.4'),
    otherphasescfgstr('utf8_filename.ncmat ; dcutoff = 0.4'),
    otherphasescfgstr('Al_sg225.ncmat;dcutoff=0.4'),
    otherphasescfgstr('Al_sg225.ncmat ; dcutoff =   0.4'),
    bad_otherphasescfgstr('Al_sg225.ncmat dcutoff=0.8'),
    bad_otherphasescfgstr('Al_sg225.ncmat;dcutoff=0.4 dcutoffup=2.0'),
    otherphasescfgstr('Al_sg225.ncmat;atomdb=H   is   D'),
    otherphasescfgstr('  Al_sg225.ncmat \t \t ; \t\t\tatomdb\t\t = H:: \t\t:  is:  :: D  '),
    otherphasescfgstr('phases<0.4*Al_sg225.ncmat;dcutoffup=2.0&0.6*void.ncmat> ;dcutoff=0.4'),
    bad_otherphasescfgstr('phases<0.4*Al_sg225.ncmat;dcutoffup=2.0&0.6*void.ncmat> dcutoff=0.4'),#missing ';' after '>'
    otherphasescfgstr('phases<0.4*Al_sg225.ncmat;dcutoffup=2.0&0.6*void.ncmat> ;dcutoff=0.4'),
    otherphasescfgstr('phases<0.4*Al_sg225.ncmat;dcutoffup=2.0&0.6*void.ncmat > ;dcutoff=0.4'),
    otherphasescfgstr('phases<0.4*Al_sg225.ncmat&0.6*void.ncmat>'),
    otherphasescfgstr(' phases   <\t\t0.4 *   Al_sg225.ncmat\t\t&0.6*void.ncmat >'),
    bad_otherphasescfgstr(''),
    #otherphasescfgstr(''),
]
for e in ll:
    if isinstance(e,tuple):
        bad = False
        ncmatdata=e[0]
    else:
        bad = True
        ncmatdata=e
    print("=====> TRYING TO LOAD NCMAT DATA:\n\n")
    print(ncmatdata)
    print("\n<=====\n\n\n")
    if not bad:
        NC.directMultiCreate(ncmatdata)
    else:
        got_err = True
        try:
            NC.directMultiCreate(ncmatdata)
            got_err = False
        except (NC.NCBadInput,NC.NCFileNotFound) as e:
            print("Got expected error: "+str(e))
        if not got_err:
            raise SystemExit('Did not end in error as expected!')
        #NC.clearCaches()
