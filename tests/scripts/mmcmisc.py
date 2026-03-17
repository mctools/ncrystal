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

# NEEDS: numpy

def main():
    from NCTestUtils.minimc_ref import main_minimc_unittest as m

    # Test ignoremiss parameter:

    for ignoremiss in [ True, False ]:
        ec='seed=4;nthreads=2;nscatlimit=0'
        if ignoremiss:
            ec += ';ignoremiss=1'
            #beam radius twice of geom radius => 75% misses.
        res = m( cfgstr='void.ncmat',
                 geomcfg='sphere;r=0.05',
                 srccfg='circular;n=1e6;z=-1;wl=1.8;r=0.1',
                 key='<auto>_%s'%('ignmiss' if ignoremiss else 'noignmiss'),
                 tally='nscat,nscat_uw',
                 extra_enginecfg=ec )
        #res.tally('nscat').dump()
        #res.tally('nscat_uw').dump()
        h1 = res.tally('nscat').hist_total
        h2 = res.tally('nscat_uw').hist_total.clone().set_title(h1.title)
        assert h1 == h2
        #print( h1.content )
        assert h1.minfilled == ( 0.0 if ignoremiss else -1.0 )
        assert h1.maxfilled == 0.0
        if ignoremiss:
            assert h1.content[0] == 0.0
        else:
            assert abs( h1.content[0]/1e6 - 0.75 ) < 0.01
        assert abs( h1.content[1]/1e6 - 0.25 ) < 0.01

    #Test absorption:

    from NCTestUtils.minimc_ref import main_minimc_unittest_stdsphere as msph

    for absorption in [ True, False ]:
        ec='nscatlimit=0'#no scattering, pure attenuation
        ec += ';tallybins=-' #less bins
        if not absorption:
            ec += ';absorption=0'
        res = msph( cfgstr='solid::B10/1gcm3',
                    illuminate_uniformly = False,
                    sphere_diam_meter = 0.01,
                    neutron_energy='4.0Aa',
                    key='<auto>_%s'%('abs' if absorption else 'noabs'),
                    extra_enginecfg = ec,
                    tally = 'mu' )
        h = res.tally('mu').hist_total
        #no scattering, mu=1 always:
        assert h.minfilled == 1.0
        assert h.maxfilled == 1.0
        nstat = res.setup['src']['decoded']['n']
        assert nstat == 1e5
        if absorption:
            #value of 5.26e-37 found in test run with nstat=1e5
            assert abs( h.integral/nstat - 5.26e-37/1e5 ) < 0.01
        else:
            print(h.integral)
            # no absorption, all weights delivered at mu=1:
            assert abs(h.integral/nstat-1.0)<0.01

    #Test seed (with nthreads=1 results should be reproducible):
    from NCrystalDev.minimc import run as mmcrun

    def runseed( seed ):
        return mmcrun( cfgstr='Al_sg225.ncmat',
                       geomcfg='slab;dz=0.05',
                       srccfg='constant;n=1e5;z=-1;wl=1.8',
                       enginecfg = f'tally=mu;nthreads=1;seed={seed}'
                      ).tally('mu')
    tallymu_1 = runseed( 12345 )
    tallymu_2 = runseed( 123456 )
    tallymu_3 = runseed( 12345 )
    h1 = tallymu_1.hist_total
    h2 = tallymu_2.hist_total
    h3 = tallymu_3.hist_total
    h1.check_compat(h2,check=True)
    h1.check_compat(h3,check=True)
    h1.check_compat(h3,threshold=0.9999,check=True)

    assert tallymu_1 == tallymu_3
    assert h1 == h3
    assert ( h1.content == h3.content ).all()
    assert ( h1.content != h2.content ).all()
    assert not ( h1 == h2 )
    assert not ( tallymu_1 == tallymu_2 )


if __name__ == '__main__':
    main()

