#!/usr/bin/env python3
import NCrystalDev._mmc as ncmmc
#import NCrystalDev._mmc_utils as ncmmc_utils
import pprint
pprint.pprint(ncmmc.available_tallies())
pprint.pprint(ncmmc.available_tallies())

#fixme unit test both "comp=" and "comp=bragg" (and beyond bragg threshold)

cfg = dict( cfgstr="Al_sg225.ncmat;temp=500K",
            geomcfg="cyl;r=0.1;dy=10",
            #geomcfg="box;dx=0.1;dy=0.1;dz=0.1",
            #geomcfg="sphere;r=0.4",
            #geomcfg="box;dx=0.4;dy=1.0;dz=1.0",
            srccfg="constant;wl=1.8;z=-1;n=1e5",
            #NCrystal: TKTEST distToCylExit. xyz=-0.0324204 5.71694 -0.0945987
            #srccfg="circular;x=0.0375266;y=6.93176;z=-1.5;r=2;wl=1.8;n=1e5",
            #srccfg="circular;wl=1.8;;x=0.097;z=-1.5;r=10.004;n=1e5", ##THIS ONE
            #srccfg="constant;wl=1.8;x=0.13;z=-1;n=1e2",
            enginecfg=(
                #";absorption=1"
                #";nthreads=1"
                #";seed=123456789"
                #";ignoremiss=1"
                #";beamdirx=0;beamdiry=0;beamdirz=1"
                ";tally=cosmu,nscat"
                #";nscatlimit=1"
                #";tallybins=e:1000:0.0:0.4"
                ";roulette=0.1,0.01,2"
            ) )

res = ncmmc.run_minimc(**cfg)
print("ENDED OK")
pprint.pprint(res._raw_data())
res.dump(prefix='>>>')
t = res.tally('nscat')
t.plot()
#t.plot( max_nbins = 100 )
