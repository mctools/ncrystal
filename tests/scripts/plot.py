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

# Dedicated test for various plotting code, to increase test coverage.

from NCrystalDev.ncmat import NCMATComposer
import NCrystalDev.plot as ncplot
import NCrystalDev.core as nccore
import NCrystalDev.vdos as ncvdos
from NCrystalDev.misc import AnyVDOS
from NCTestUtils.env import ncsetenv
from NCrystalDev._numpy import _np_linspace
from NCTestUtils.common import ensure_error
import NCTestUtils.reprint_escaped_warnings # noqa F401
from NCrystalDev.exceptions import NCBadInput, NCLogicError

def main(do_plot):
    if not do_plot:
        ncsetenv('FAKEPYPLOT','1')
    test_plotcontext()
    test_vdos()
    test_xsects()
    test_knl()

def test_plotcontext():
    print(">>>> PlotContext test #1")
    with ensure_error(NCBadInput,
                      'explicit do_show not supported with plot_context'):
        ncplot.plot_xsect('void.ncmat',
                          do_show=True,
                          plot_context=ncplot.PlotContext())
    print(">>>> PlotContext test #2")

    with ensure_error(NCLogicError,
                      'Do not use PlotContext object after it has already been'
                      ' consumed. It is consumed when its .finalise() method is'
                      ' invoked or when it is passed directly as the value of a'
                      ' plot_context argument (you most likely want to pass it'
                      ' as plot_context=obj.subcontext() instead in this'
                      ' case).'):
        pctx = ncplot.PlotContext()
        ncplot.plot_xsect('void.ncmat',plot_context=pctx)
        print("Now call pctx.finalise")
        pctx.finalise()

    print(">>>> PlotContext test #3")
    pctx = ncplot.PlotContext()
    ncplot.plot_xsect('void.ncmat',plot_context=pctx.subcontext())
    pctx.finalise()

def test_vdos():
    egrid = _np_linspace(100.0,400.0,17)*1e-3#eV
    vdos1 = (egrid,(egrid/0.2)**5*400+500)
    vdos2 = (egrid,(egrid/0.3)**3*(-400)+1000)
    anyvdos1 = AnyVDOS( vdos1, label='vdos1' )
    anyvdos2 = AnyVDOS( vdos2, label='vdos2' )
    print(">>>> ncplot.plot_vdos #1")
    ncplot.plot_vdos( vdos1, vdos2, unit='meV' )
    def testlblfct( lbl ):
        print(lbl)
        return lbl
    print(">>>> ncplot.plot_vdos #2")
    ncplot.plot_vdos( vdos1, anyvdos2, unit='meV',
                      labelfct=testlblfct,
                      show_orig_data = True )

    print(">>>> Test PhononDOSAnalyser")
    dosana = ncvdos.PhononDOSAnalyser( [ anyvdos1, anyvdos2 ] )
    print(">>>> PhononDOSAnalyser.plot #1")
    dosana.plot()
    print(">>>> PhononDOSAnalyser.plot #1b")
    dosana.plot(logy=True)
    print(">>>> PhononDOSAnalyser.plot_cutoff_effects #1")
    dosana.plot_cutoff_effects( [0.1,0.2,0.3] )
    print(">>>> PhononDOSAnalyser.plot_cutoff_effects #1b")
    c = NCMATComposer('stdlib::Al2O3_sg167_Corundum.ncmat')
    dosana.plot_cutoff_effects_on_xsects( c,[0.1,0.2,0.3],
                                          lblmap={'vdos1':'Al','vdos2':'O'})
    print(">>>> PhononDOSAnalyser.plot_cutoff_effects #2")
    dosana.plot_cutoff_effects( [0.1,0.2,0.3], gn = 1, masses=[10.0,20.0] )
    dosana = dosana.apply_cutoff( 0.2 )
    print(">>>> PhononDOSAnalyser.plot #2 (after apply_cutoff)")
    dosana.plot()

    mat = NCMATComposer()
    mat.set_density(2.0,'g/cm3')
    mat.set_dyninfo_vdos('Al', *dosana.dos('vdos2'), fraction=0.3 )
    mat.set_dyninfo_vdos('O', *vdos1, fraction=0.7 )
    info = mat.loadInfo()
    print(">>>> Plotting di with show_orig_data=True")
    for di in info.dyninfos:
        print(">>>> di.plot_vdos %s"%di.atomData.displayLabel())
        di.plot_vdos(show_orig_data=True)

def test_xsects():
    #TODO: Way more tests (also, the ylim/xlim logic is not sound inside the
    #function):
    print(">>>> plot_xsect #1")
    ncplot.plot_xsect('stdlib::Al_sg225.ncmat;comp=incoh_elas')
    print(">>>> plot_xsect #2 (macro)")
    ncplot.plot_xsect('stdlib::Al_sg225.ncmat;comp=incoh_elas;density=2x',
                      xsmode='macroscopic')
    print(">>>> plot_xsect #3 (void ekin)")
    ncplot.plot_xsect('void.ncmat',logy=True,mode='ekin')
    print(">>>> plot_xsect #3 (ymin>ymax)")
    ncplot.plot_xsect('void.ncmat',ymin=5.0,ymax=0.1)
    print(">>>> plot_xsect #4 (no procs)")
    with ensure_error(NCBadInput,
                      'Can not produce plots for material source which'
                      ' contains no physics processes'):
        ncplot.plot_xsect(
            NCMATComposer('stdlib::Al_sg225.ncmat').load(doScatter=False,
                                                         doAbsorption=False)
        )
    print(">>>> plot_xsect #5 (only absorption 1/2)")
    ncplot.plot_xsect(
        NCMATComposer('stdlib::Al_sg225.ncmat').load(doInfo=False,
                                                     doScatter=False)
    )
    print(">>>> plot_xsect #6 (only absorption 2/2)")
    ncplot.plot_xsect('stdlib::Al_sg225.ncmat',
                      show_scattering=False)
    print(">>>> plot_xsect #7 (only scatter 1/2)")
    ncplot.plot_xsect(
        NCMATComposer('stdlib::Al_sg225.ncmat').load(doInfo=False,
                                                     doAbsorption=False)
    )
    print(">>>> plot_xsect #8 (only scatter 2/2)")
    ncplot.plot_xsect('stdlib::Al_sg225.ncmat',
                      show_absorption=False)

    print(">>>> plot_xsect #9 (macroscopic fail)")
    with ensure_error(NCBadInput,
                      'Can not produce macroscopic cross section plots for'
                      ' material source which contains no NCrystal.Info '
                      'object (and therefore no material density).'):
        ncplot.plot_xsect(
            NCMATComposer('stdlib::Al_sg225.ncmat').load(doInfo=False),
            xsmode='macroscopic'
        )

    print(">>>> plot_xsect #10 (bad plotcontext)")
    with ensure_error(NCBadInput,
                      'explicit do_show not supported with plot_context'):
        ncplot.plot_xsect('stdlib::Al_sg225.ncmat',
                          do_show=True,
                          plot_context=ncplot.PlotContext())
    print(">>>> plot_xsect #11 (bad arg)")
    with ensure_error(NCBadInput,
                      'Unsupported argument: foobar=123'):
        ncplot.plot_xsect('stdlib::Al_sg225.ncmat',
                          foobar=123)

    print(">>>> plot_xsects #1")
    ncplot.plot_xsects('stdlib::Al_sg225.ncmat;comp=incoh_elas',
                       'stdlib::Be_sg194.ncmat;comp=incoh_elas')
    print(">>>> plot_xsects #2")
    ncplot.plot_xsects('void.ncmat',title='foo')
    print('estimate_longest_interesting_wavelength multiphase = %g'
          %ncplot._estimate_longest_interesting_wavelength(
              nccore.createInfo(
                  'phases<0.1*stdlib::Al_sg225.ncmat'
                  '      &0.5*stdlib::MgF2_sg136_MagnesiumFlouride.ncmat'
                  '      &0.4*gasmix::air>')
          ))



def test_knl():
    mat = NCMATComposer('stdlib::Al_sg225.ncmat')
    mat._unofficial_vdos2sab_ignore( order_low=2, order_high=9999 )
    info = mat.loadInfo()
    di = info.dyninfos[0]
    print(">>>> di.plot_Gn(1)")
    di.plot_Gn(1)
    print(">>>> di.plot_knl")
    di.plot_knl(vdoslux=0,phasespace_curves=[0.001,0.1])

if __name__ == '__main__':
    import sys
    main( do_plot = '--plot' in sys.argv[1:] )
