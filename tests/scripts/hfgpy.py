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

import NCTestUtils.enable_fpe # noqa F401
from NCrystalDev.hfg2ncmat import hfg2ncmat
from NCrystalDev.ncmat import NCMATComposer

def test(spec,formula,**kwargs):
    import pprint
    print()
    print()
    print()
    print('-'*80)
    if 'title' not in kwargs:
        kwargs['title'] = 'Dummy title'
    if 'density' not in kwargs:
        kwargs['density'] = 123.45
    kwargs['spec']=spec
    kwargs['formula']=formula
    pprint.pprint(kwargs)
    res=hfg2ncmat(**kwargs)
    print(res)

def test_fail( *args, **kwargs ):
    from NCrystalDev.exceptions import NCBadInput
    try:
        test(*args,**kwargs)
    except NCBadInput as e:
        print("FAILED (as expected): %s"%e)
        return
    raise SystemExit('Did not fail as expected')

test('1xCHali+1xCHaro+1xCH2+1xCH3+1xNH+1xNH2+1xNH3+1xOH+1xSH','H15')
test('1xCHali+1xCHaro+1xCH2+1xCH3+1xNH+1xNH2+1xNH3+1xOH+1xSH','H15C4N3OS')
test('5xCHaro+1xCHali+1xCH2','C8H8',density=0.99)
test_fail('1xCH2','CH2s')

c = NCMATComposer.from_hfg( '5xCHaro+1xCHali+1xCH2',
                            'C8H8',
                            density=0.99,
                            title='polystyrene' )
#c.inspect()
print(c.create_ncmat())
