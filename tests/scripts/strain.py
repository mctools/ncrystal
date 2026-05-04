#!/usr/bin/env python3

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2026 NCrystal developers                                   ##
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

import NCrystalDev.core as nccore
from NCTestUtils.common import ensure_error
from NCrystalDev.exceptions import NCBadInput

from NCTestUtils.sab import extract_scathelper_uids

def test( cfgstr ):
    print('\n'*3)
    print('-------------- Testing %s --------------'%repr(cfgstr))
    m = nccore.load(cfgstr)
    m_bragg = nccore.load(cfgstr+';comp=bragg')
    m_nonbragg = nccore.load(cfgstr+';bragg=0')

    print("--> scat uid = %i"%m.scatter.uid)
    print("--> nonbragg scat uid = %i"%m_nonbragg.scatter.uid)
    print("--> bragg scat uid = %i"%m_bragg.scatter.uid)
    m.info.dump()
    wls = [3.0,3.5,4.0,4.5,5.0,5.5,6.0]
    for wl in wls:
        print('  Bragg XS(%3gAa) = %g barn'%
              (wl,m_bragg.scatter.xsect( wl=wl )))
    for wl in wls:
        print('  NonBragg XS(%3gAa) = %g barn'%
              (wl,m_nonbragg.scatter.xsect( wl=wl )))
    scathelper_uids_base = extract_scathelper_uids( m.scatter )
    scathelper_uids_nonbragg = extract_scathelper_uids( m_nonbragg.scatter )
    scathelper_uids_bragg = extract_scathelper_uids( m_bragg.scatter )
    assert len(scathelper_uids_base) == 1
    assert len(scathelper_uids_nonbragg) == 1
    assert len(scathelper_uids_bragg) == 0
    assert scathelper_uids_base[0] == scathelper_uids_nonbragg[0]

    return dict( numdens = m.info.numberdensity,
                 braggthr = m.info.braggthreshold,
                 scat_uid_bragg = m_bragg.scatter.uid,
                 scat_uid_nonbragg = m_nonbragg.scatter.uid,
                 scathelper_uid = scathelper_uids_base[0],
                 nonbragg_xsvals = [ m_nonbragg.scatter.xsect( wl=wl )
                                     for wl in wls ] )


def assert_equal(a,b):
    assert abs(a-b)<1e-6
def assert_equal_lists(a,b):
    assert len(a)==len(b)
    for aa,bb in zip(a,b):
        assert_equal(aa,bb)

def main():
    basecfg = 'stdlib::Al_sg225.ncmat;dcutoff=1.5'
    base = test(basecfg)
    strainp =test(basecfg+';strain=0.2')
    strainm =test(basecfg+';strain=-0.2')

    assert_equal(base['numdens']/1.2**3,strainp['numdens'])
    assert_equal(base['braggthr']*1.2,strainp['braggthr'])
    assert_equal(base['numdens']/0.8**3,strainm['numdens'])
    assert_equal(base['braggthr']*0.8,strainm['braggthr'])

    #TODO: Uncomment these checks as well!
    #print('A',base['nonbragg_xsvals'])
    #print('B',strainp['nonbragg_xsvals'])
    assert_equal_lists(base['nonbragg_xsvals'],strainp['nonbragg_xsvals'])
    assert_equal_lists(base['nonbragg_xsvals'],strainm['nonbragg_xsvals'])
    assert base['nonbragg_xsvals'] == strainp['nonbragg_xsvals']
    assert base['nonbragg_xsvals'] == strainm['nonbragg_xsvals']
    #assert base['scat_uid_nonbragg'] == strainp['scat_uid_nonbragg']
    #assert base['scat_uid_nonbragg'] == strainm['scat_uid_nonbragg']
    assert base['scat_uid_bragg'] != strainp['scat_uid_bragg']
    assert base['scat_uid_bragg'] != strainm['scat_uid_bragg']
    assert strainp['scat_uid_bragg'] != strainm['scat_uid_bragg']
    assert base['scathelper_uid'] == strainp['scathelper_uid']
    assert base['scathelper_uid'] == strainm['scathelper_uid']

    with ensure_error(NCBadInput,
                      'strain must be in the interval [-0.5,0.5]'):
        nccore.createInfo(basecfg+';strain=0.6')
    with ensure_error(NCBadInput,
                      'strain must be in the interval [-0.5,0.5]'):
        nccore.createInfo(basecfg+';strain=-10.0')
    with ensure_error(NCBadInput,
                      'Syntax error - invalid value "none" '
                      'provided for parameter "strain"'):
        nccore.createInfo(basecfg+';strain=none')
    with ensure_error(NCBadInput,
                      'Syntax error - invalid value "" '
                      'provided for parameter "strain"'):
        nccore.createInfo(basecfg+';strain =  ')



if __name__ == '__main__':
    main()
