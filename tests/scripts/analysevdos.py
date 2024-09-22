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

# NEEDS: numpy

import NCrystal as NC

itest = 0
for vdoslux in (0,3):
  for temp in (0.001,300,1e6):
    for ncfile in sorted(NC.browseFiles(factory='stdlib')):
      if 'Liquid' in ncfile.name:
        continue
      itest += 1
      if not itest%43==0:
        continue
      cfgstr=f'{ncfile.name};temp={temp};vdoslux={vdoslux}'
      info=NC.createInfo(cfgstr)
      for di in info.dyninfos:
        if hasattr(di,'analyseVDOS'):
          print(f'==> {cfgstr}//{di.atomData.displayLabel()} => ',end='')
          print('; '.join('%s=%.12g'%(k,v) for k,v in sorted(di.analyseVDOS().items())))
