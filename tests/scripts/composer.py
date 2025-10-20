#!/usr/bin/env python3

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

import NCrystalDev as NC
from NCTestUtils.common import ensure_error
import pathlib

def main():
    print('\n\n  ================> He Ne gas example\n\n')

    c_gas = NC.NCMATComposer()
    c_gas.set_dyninfo_freegas( 'He', fraction = 0.91 )
    c_gas.set_dyninfo_freegas( 'Ne', fraction = 0.09 )
    c_gas.set_state_of_matter( 'gas' )#to be expressive, this line not strictly needed
    c_gas.set_density(0.233642,'kg/m3')

    c_gas_ncmat = c_gas()
    print(c_gas_ncmat)

    print("\nc_gas.load(doScatter=False,doAbsorption=False).dump() ::\n")
    c_gas.load(doScatter=False,doAbsorption=False).dump()
    print("\nc_gas.load(doInfo=False,doAbsorption=False).dump() ::\n")
    c_gas.load(doInfo=False,doAbsorption=False).dump()
    print("\nc_gas.load(doInfo=False,doScatter=False).dump() ::\n")
    c_gas.load(doInfo=False,doScatter=False).dump()
    print("\nc_gas.load().dump() ::\n")
    c_gas.load().dump()
    c_gas.write('mysillygas.ncmat')
    assert pathlib.Path('mysillygas.ncmat').read_text() == c_gas_ncmat

    scatter = NC.createScatter("mysillygas.ncmat;temp=50K")
    scatter.dump()

    c_gas.register_as('mysillygas2.ncmat')
    scatter = NC.createScatter("mysillygas2.ncmat;temp=50K")
    scatter.dump()


    print('\n\n  ================> silly PE example\n\n')

    c_PE = NC.NCMATComposer()
    #help(c_PE.set_dyninfo_msd)
    c_PE.set_dyninfo_msd('H',msd=0.022,temperature=50,fraction=2/3)
    c_PE.set_dyninfo_msd('C',msd=0.0076,temperature=50,fraction=1/3)
    c_PE.set_density(0.92,'g/cm3')
    print(c_PE.create_ncmat())


    print('\n\n  ================> U example\n\n')

    c_u = NC.NCMATComposer("solid::U/1gcm3")#note: can also init from a cfg-string!
    c_u.set_composition('U','U232')
    with ensure_error(NC.NCBadInput,
                      'Atom with label "U232" is unknown. If it is a valid'
                      ' isotope which is simply missing in NCrystal\'s internal'
                      ' database you must define it yourself.'):
        c_u.load()


    c_u.update_atomdb('U232','232.04u -12.3fm 4.5b 6.789b')#NB: dummy data values!
    print(c_u())
    c_u.load().dump()
    #alternative syntax:
    print(c_u())
    c_u.update_atomdb('U232',mass=232.04, coh_scat_len=-12.3, incoh_xs=4.5, abs_xs=6.789)
    c_u.load().dump()


if __name__ == '__main__':
    main()
