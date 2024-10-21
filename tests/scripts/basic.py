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
from NCTestUtils.printnumpy import format_numpy_1darray_asfloat as npfmt
import NCrystal as NC
import pathlib
import os
import sys
import pprint
import math
import numpy as np
from numpy import set_printoptions as np_setprintopts

np_setprintopts(infstr='inf')#test reproducibility

#NC.disableCaching()
NC.setDefaultRandomGenerator(None)#better test reproducibility


#conversions
for f in 'ekin2wl wl2ekin ekin2ksq wl2k wl2ksq'.split():
    for val in 0.5, 1.0, 2.0, 0.0, float('inf'):
        print(f'Conversion: NC.{f}({val}) = {getattr(NC,f)(val):.14g}')

#non-oriented scatters:

datafile = "Al_sg225.ncmat"
info = NC.createInfo(datafile)

print('get_build_namespace=%s'%NC.get_build_namespace())

def myprint(*args):
    sys.stdout.flush()
    sys.stderr.flush()
    print(*args)
    sys.stdout.flush()
    sys.stderr.flush()


myprint("Temperature : %s"%('%g'%info.getTemperature() if info.hasTemperature() else '<n/a>'))
myprint("DebyeTemp : %s"%('yes' if info.hasDebyeTemperature() else 'no'))
myprint("XSectAbsorp : %g"%info.getXSectAbsorption())
myprint("XSectFree : %g"%info.getXSectFree())
myprint("Density : %g"%info.getDensity())
myprint("Structure : %s"%(', '.join('%s=%.12g'%(k,v) for k,v in sorted(info.getStructureInfo().items())) if info.hasStructureInfo() else '<n/a>'))

calc_pcbragg = NC.createScatter(datafile+";inelas=0;incoh_elas=0")
xs = calc_pcbragg.crossSectionIsotropic(NC.wl2ekin(4.0))#FIXME obsolete method
myprint("Aluminium %s diffraction cross-section @ 4.0Aa: %g barn"%(calc_pcbragg.name,xs))

wls=np.linspace(1.0,4.0,4)
myprint("Aluminium %s diffraction cross-section @ [%s]Aa: [%s] barn"%(calc_pcbragg.name,
                                                                      (', '.join('%g'%e for e in wls)),
                                                                      (', '.join('%g'%e for e in calc_pcbragg.xsect(wl=wls)))))


calc_bkgd = NC.createScatter(datafile+";bragg=false")
xs = calc_bkgd.crossSectionIsotropic(NC.wl2ekin(4.0))
myprint("Aluminium Bkgd cross-section @ 4.0Aa: %g barn"%xs)

calc_scomp = NC.createScatter(datafile)
xs = calc_scomp.crossSectionIsotropic(NC.wl2ekin(4.0))
myprint("Aluminium Bkgd+PCBragg cross-section @ 4.0Aa: %g barn"%xs)

xs_abs = NC.createAbsorption(datafile).crossSectionIsotropic(NC.wl2ekin(4.0))
myprint("Aluminium absorption cross-section @ 4.0Aa: %g barn"%xs_abs)

#single crystal Ge (also testing using python multi-line string for complex settings):
calc_scatfactsc = NC.createScatter("""
                                   Ge_sg227.ncmat;
                                   mos=0.001;
                                   dir1=@crys_hkl:5,1,1@lab:1,0,0;
                                   dir2=@crys_hkl:0,0,1@lab:0,1,1;
                                   dirtol=45deg""")

xs = calc_scatfactsc.crossSection( NC.wl2ekin(4.0), (1,0.15,0.7) )
myprint("Germanium ScatFactSC cross-section @ 4.0Aa and (1.0, 0.15, 0.7): %g barn"%xs)
myprint("... same with xsect method: %g barn"%calc_scatfactsc.xsect( wl=4.0, direction=(1,0.15,0.7) ))
myprint("... same with xsect method2: %g barn"%calc_scatfactsc.xsect( NC.wl2ekin(4.0), direction=(1,0.15,0.7) ))
myprint("... same with xsect method3: %g barn"%calc_scatfactsc.xsect( direction=(1,0.15,0.7),ekin=NC.wl2ekin(4.0) ) )



myprint("A few generated scattering angles using the full 3d vector interface:")
_iii = [0]
def test_genscatter(scatter,wl,indir):
    global _iii
    _iii[0] += 1
    select = _iii[0]%4
    if select == 0:
        #tuple(ekin_final,direction_final)
        out_ekin, outdir = scatter.sampleScatter(NC.wl2ekin(wl),indir)
    elif select == 1:
        out_ekin, outdir = scatter.scatter(NC.wl2ekin(wl),indir)
    elif select == 2:
        out_ekin, outdir = scatter.scatter(ekin=NC.wl2ekin(wl),direction=indir)
    else:
        assert select == 3
        out_ekin, outdir = scatter.scatter(wl=wl,direction=indir)
    delta_e = out_ekin - NC.wl2ekin(wl)

    x,y,z=indir
    a,b,c=outdir
    angle = math.acos((x*a+y*b+z*c)/math.sqrt((x*x+y*y+z*z)*(a*a+b*b+c*c)))
    indir_str = '(%g, %g, %g)'%indir#like this rather than str(indir) for python 2.6 test reproducibility
    myprint("scattering at %g Aa from %s happens at %g deg and delta-E %g eV"%(wl,indir_str,angle*180/math.pi,delta_e))

for i in range(20):
    test_genscatter(calc_scatfactsc,4.0, (1,0.15,0.7))
for i in range(20):
    test_genscatter(calc_scatfactsc,4.0, (1,0.15,0.6))
#angles,de = calc_bkgd.generateScatteringNonOriented(NC.wl2ekin(4.0),repeat=10000)
ekin_final, mu = calc_bkgd.sampleScatterIsotropic(NC.wl2ekin(4.0),repeat=10000)
angles, de = np.vectorize(math.acos)(mu), ekin_final - NC.wl2ekin(4.0)

myprint(npfmt(angles*(180/math.pi)))
myprint(npfmt(de*1000))

ekin_final, mu = calc_bkgd.scatter(wl=4.0,repeat=10000)
angles,de = np.vectorize(math.acos)(mu), ( ekin_final-NC.wl2ekin(4.0) )
myprint(npfmt(angles*(180/math.pi)))
myprint(npfmt(de*1000))

#vanishing processes
for fn in ('Be_sg194.ncmat','Al_sg225.ncmat'):
    for part in ('',";inelas=0;incoh_elas=0",';coh_elas=0'):
        cfg=fn+';dcutoff=10.0'+part
        myprint('===> Testing potentially vanishing scatter: "%s":'%cfg)
        sc = NC.createScatter(cfg)
        myprint('    xs@1Aa : %g'%sc.crossSectionIsotropic(NC.wl2ekin(1.0)))
        myprint('    xs@5Aa : %g'%sc.crossSectionIsotropic(NC.wl2ekin(5.0)))
        #ang,de=sc.generateScatteringNonOriented(NC.wl2ekin(1.0),repeat=10000)
        ekin_final,mu=sc.sampleScatterIsotropic(NC.wl2ekin(1.0),repeat=10000)
        ang, de = np.vectorize(math.acos)(mu), ekin_final - NC.wl2ekin(1.0)
        myprint('    mean scat@1Aa : %.2g rad, de=%.2geV'%(ang.mean(),de.mean()))

myprint('===> Testing exceptions')

#test catch all ncrystal exceptions:
try:
    NC.createInfo("nosuchfile.ncmat")
except NC.NCException as e:
    myprint("Caught NCException!:")
    myprint('  Type    : ',e.__class__.__name__)
    myprint('  Message : ',e.message)

#test catch specifically a given type:
try:
    NC.createInfo("nosuchfile.ncmat")
except NC.NCFileNotFound as e:
    myprint("Caught NCFileNotFound!:")
    myprint('  Type    : ',e.__class__.__name__)
    myprint('  Message : ',e.message)

#or the wrong type:
try:
    caught = False
    try:
        NC.createInfo("nosuchfile.ncmat")
    except NC.NCCalcError as e:
        caught = True
        myprint("Caught NCCalcError!:")
        myprint('  Type    : ',e.__class__.__name__)
        myprint('  Message : ',e.message)
except Exception:
    pass
assert not caught
myprint("Did not catch exception! (as expected)")

###extract packingfactor:
##def _testpf(cfgstr):
##    print('decodepackingfactor("%s") = %g'%(cfgstr,NC.decodecfg_packfact(cfgstr)))
##_testpf("Al_sg225.ncmat;dcutoff=0.5")
##_testpf("Al_sg225.ncmat;packfact=0.235;dcutoff=0.5")
###_testpf("Al_sg225.ncmat;packfact=0.235;dcutoff=0.5;packfact=0.6;")

#ensure proper memory cleanup by releasing internal default-assigned random
#generator:
NC.setDefaultRandomGenerator(None)

myprint("wl2ekin:")
myprint('%.12g'%NC.wl2ekin(4.0))
myprint(NC.wl2ekin(0.0))
myprint(NC.wl2ekin(float('inf')))
myprint(npfmt(NC.wl2ekin(np.asarray([4.0,0.0,float('inf')]))))
myprint("ekin2wl:")
myprint('%.12g'%NC.ekin2wl(0.025))
myprint(NC.ekin2wl(0.0))
myprint(NC.ekin2wl(float('inf')))
myprint(npfmt(NC.ekin2wl(np.asarray([0.025,0.0,float('inf')]))))

#Integrated unit test:
NC.test()


#getBraggThreshold:
#same file, different temp => no Info cache but should be identical bragg threshold:
print("\n---> testing fast Bragg threshold detection:")
sc1 = NC.createInfo('NaCl_sg225_SodiumChloride.ncmat;temp=123K')
sc2 = NC.createInfo('NaCl_sg225_SodiumChloride.ncmat;temp=234K')
fast_bt = sc1.braggthreshold
print("  Fast Bragg threshold(1) %.7g"%fast_bt)
alt_bt = 2.0 * next(sc2.hklList())[4]#twice the dspacing of the first hkl entry
print("  Full init bragg threshold(1) %.7g"%alt_bt)
assert abs(fast_bt-alt_bt)<1e-3
assert abs(fast_bt-alt_bt)<1e-14
assert abs(sc2.braggthreshold-alt_bt)<1e-14
print()
#test paths with unicode characters:

dirname=u'unicodedir_test\u4500abc'
testdir=os.path.join(os.getcwd(),dirname)
try:
    os.makedirs(testdir)
except UnicodeEncodeError:
    #file system encoding does not allow our special character - avoid test
    #failure by falling back to ascii name:
    dirname=u'asciidir_testXabc'
    testdir=os.path.join(os.getcwd(),dirname)
    os.makedirs(testdir)

al_sg225_content = NC.createTextData('stdlib::Al_sg225.ncmat').rawData
(pathlib.Path(testdir)/'Al_sg225.ncmat').write_text(al_sg225_content)

def testLoadOK(cfgstr,expectBadInput = False):
    print(f"Test loading of :>>>{cfgstr}<<<")
    try:
        NC.createInfo(cfgstr)
    except NC.NCBadInput as e:
        _='expected' if expectBadInput else 'unexpected'
        myprint(f"Caught {_} NCBadInput!: {e.message}")
        if not expectBadInput:
            raise SystemExit('Unexpected failure')
        return
    if expectBadInput:
        raise SystemExit('Did not fail as expected')

def testLoadFails(cfgstr):
    return testLoadOK(cfgstr,expectBadInput = True)

specialfn=os.path.join(dirname,'Al_sg225.ncmat')
testLoadOK(specialfn)
testLoadFails(f'{specialfn};inelas=test\u4500abc')
testLoadOK(f'phases<1.0*{specialfn}>')
testLoadOK(f'phases<1.0*{specialfn}>')
testLoadOK(f'phases<0.1*{specialfn}&0.9*void.ncmat>')

info = NC.createInfo(os.path.join(dirname,'Al_sg225.ncmat'))
assert info.hasTemperature() and abs(info.getTemperature()-293.15)<1e-10

for temp in [1,300,10000]:
    for mass in [1.007, 10.0, 300.0]:
        for debye_temp in [5., 400., 10000.]:
            msd = NC.debyeIsotropicMSD( debye_temperature=debye_temp, temperature=temp, mass=mass )
            print(f'  debyeIsotropicMSD( TDebye={debye_temp}, T={temp}, mass={mass} ) = {msd}')
            debye_temp_calc = NC.debyeTempFromIsotropicMSD(msd=msd, temperature=temp, mass=mass)
            print(f'  ...going back via debyeTempFromIsotropicMSD gives TDebye={debye_temp_calc}')
            assert debye_temp_calc>0.99*debye_temp
            assert abs(debye_temp-debye_temp_calc)/debye_temp < 0.0001


#test tht json encoding of cfg docs can be decoded by the python json module:
cfgaspy = NC.generateCfgStrDoc('python')
pprint.pprint(cfgaspy)

def get_par_doc(parname):
    for e in cfgaspy:
        for p in e['parameters']:
            if p['name']==parname:
                return p
    raise RuntimeError(f'No such parameter {parname}')
assert isinstance( get_par_doc('dcutoff')['default_value'], float )
assert isinstance( get_par_doc('dcutoffup')['default_value'], float )
assert isinstance( get_par_doc('mosprec')['default_value'], float )
assert isinstance( get_par_doc('vdoslux')['default_value'], int )
assert isinstance( get_par_doc('mosprec')['default_value_str'], str )
assert isinstance( get_par_doc('dcutoff')['default_value_str'], str )
assert isinstance( get_par_doc('dcutoffup')['default_value_str'], str )
assert isinstance( get_par_doc('vdoslux')['default_value_str'], str )
