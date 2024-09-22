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

import NCrystal as NC

def test(cfgstr,bad=False,**fcargs):
    print(f'CFGSTR: >>>{cfgstr}<<<')
    def name(v):
        return f'FCT::{v.__name__}' if '__name__' in dir(v) else str(v)
    print(f'ARGS: >>>{", ".join(["%s=%s"%(k,name(v)) for k,v in sorted(fcargs.items())])}<<<')
    info = NC.createInfo( cfgstr )
    try:
        fc = info.getFlattenedComposition( **fcargs )
    except NC.NCCalcError as e:
        if bad:
            print( f'  -> Got expected failure: >>>{e}<<<' )
            return
        else:
            raise
    if bad:
        raise SystemExit('  -> Did not fail as expected!')
    print(f'  -> flat compos: >>>{fc}<<<')

def testBad(cfgstr,**fcargs):
    return test(cfgstr,bad=True,**fcargs)


test('Al_sg225.ncmat')
test('Al2O3_sg167_Corundum.ncmat')
test('Al2O3_sg167_Corundum.ncmat;atomdb=Al is O')
test('Al2O3_sg167_Corundum.ncmat;atomdb=Al is Al27')

testBad('Al2O3_sg167_Corundum.ncmat;atomdb=Al is O18')

def simpleDB( Z ):
    #Just a few elements here:
    if Z == 8:
        return [ (16, 0.9975), (17, 0.0004), (18, 0.0021) ]
    if Z == 13:
        return [ (27, 1.0) ]
    if Z==92:
        return [ (238,0.9), (235,0.1) ]
test('Al2O3_sg167_Corundum.ncmat;atomdb=Al is O18',naturalAbundProvider=simpleDB)

def simpleDB2( Z ):
    #silly, but easy to calculate. Missing Al.:
    if Z == 8:
        return [ (16, 0.5), (17, 0.5) ]

test('Al2O3_sg167_Corundum.ncmat;atomdb=Al is O18',naturalAbundProvider=simpleDB2)

test('Al2O3_sg167_Corundum.ncmat;atomdb=O is Al27',naturalAbundProvider=simpleDB)
testBad('Al2O3_sg167_Corundum.ncmat;atomdb=O is Al27',naturalAbundProvider=simpleDB2)#missing Al

testBad('UO2_sg225_UraniumDioxide.ncmat;atomdb=O is 0.4 U235 0.5 U238 0.1 U234')
test('UO2_sg225_UraniumDioxide.ncmat;atomdb=O is 0.4 U235 0.5 U238 0.1 U234',naturalAbundProvider=simpleDB)

testBad('Al2O3_sg167_Corundum.ncmat',preferNaturalElements=False)
testBad('Al2O3_sg167_Corundum.ncmat',preferNaturalElements=False,naturalAbundProvider=simpleDB2)
test('Al2O3_sg167_Corundum.ncmat',preferNaturalElements=False,naturalAbundProvider=simpleDB)

def simpleDB3( Z ):
    #like DB2 but adds component with frac=0:
    if Z == 8:
        return [ (16, 0.5), (17, 0.5), (18,0.0) ]

test('Al2O3_sg167_Corundum.ncmat;atomdb=Al is O16',naturalAbundProvider=simpleDB3)

def simpleDB2_inexact( Z ):
    #like DB2 but fractions do not add exactly to unity - but almost:
    if Z == 8:
        return [ (16, 0.499999999999), (17, 0.499999999999) ]

def simpleDB2_bad( Z ):
    #like DB2 but fractions do not add to unity:
    if Z == 8:
        return [ (16, 0.5), (17, 0.4999) ]

test('Al2O3_sg167_Corundum.ncmat;atomdb=Al is O18',naturalAbundProvider=simpleDB2_inexact)
testBad('Al2O3_sg167_Corundum.ncmat;atomdb=Al is O18',naturalAbundProvider=simpleDB2_bad)
