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

import NCTestUtils.enable_fpe # noqa F401
from NCrystalDev.cfgstr import decodeCfg, normaliseCfg
from NCTestUtils.common import ensure_error
from NCrystalDev.exceptions import NCBadInput
import pathlib
import pprint

def t( extncfg_val, expected_norm_val ):
    print(f'\n---> Testing {extncfg_val}:')

    base = "dummy.ncmat;extn="
    cfgstr = f'{base}{extncfg_val}'
    d = decodeCfg(cfgstr)
    assert d['format'] == 'NCrystal-MatCfg-v1'
    assert len(d['pars'])==1 and d['pars'][0][0]=='extn'
    extn_json_dict = d['pars'][0][1]
    print("extn decoded json:")
    pprint.pprint(extn_json_dict)
    nc = normaliseCfg(cfgstr)
    assert nc.startswith(base)
    nc = nc[len(base):]
    print(f'Normalisation of extn="{extncfg_val}" gave "{nc}"')
    assert expected_norm_val is not None
    if nc != expected_norm_val:
        raise SystemExit('Expected normalisation of extn cfg'
                         f' "{extncfg_val}" to yield'
                         f' "{expected_norm_val}" but gave "{nc}"')

def main():
    pathlib.Path('dummy.ncmat').touch()
    #fixme: test "1.0mu" as well.
    t("2.0","2.0")#fixme: "2.0" or "2" ?
    t("2","2")#fixme: "2.0" or "2" ?
    t("0.1","0.1")
    t("0.1Aa","0.1")
    t("0.1 m","0.1m")

    with ensure_error( NCBadInput,
                       'Invalid domain size (length) length'
                       ' value in cfg-string: "0.1m2"' ):
        t("0.1m2",None)

    t("1e-6m/mdl:sabine","1e-6m/mdl:sabine")
    t("1e-6m","1e-6m")
    #t("1e-6m/mdl:default","1e-6m")

    with ensure_error( NCBadInput,
                       'Syntax error in extinction cfg "2.0/mdl:default":'
                       ' Unknown model "default".'):
        t("2.0/mdl:default",None)

    t("2.0/mdl:sabine/corr:1","2.0/mdl:sabine")
    t("2.0/mdl:sabine/corr:0","2.0/mdl:sabine/corr:0")
    t("2.0/mdl:sabine/tilt:rec","2.0/mdl:sabine")
    t("2.0/mdl:sabine/tilt:tri","2.0/mdl:sabine/tilt:tri")
    t("2.00000000000000000000000  /   mdl:sabine  / tilt:tri ","2/mdl:sabine/tilt:tri")

    t("3.0mu/mdl:bc","3.0mu/mdl:bc")
    t("3.0Aa/mdl:bc","3.0/mdl:bc")
    t("3.0Aa/mdl:bc/yp:cls","3.0/mdl:bc/yp:cls")
    t("3.0Aa/mdl:bc/yp:ucls","3.0/mdl:bc/yp:ucls")
    t("3.0Aa/mdl:bc/yp:lux","3.0/mdl:bc")#the default

    with ensure_error( NCBadInput,
                       'Syntax error in extinction cfg "2.0/mdl:bc/yp:bla": '
                       'Value of "yp" must be "lux", "cls", or "ucls".' ):
        t("2.0/mdl:bc/yp:bla",None)

    with ensure_error( NCBadInput,
                       'Syntax error in extinction cfg "2.0/mdl:bc/foo:bar": '
                       'Model "bc" does not support a "foo" parameter.' ):
        t("2.0/mdl:bc/foo:bar",None)

if __name__ == '__main__':
    main()
