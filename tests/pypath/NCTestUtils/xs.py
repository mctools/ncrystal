
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

# Idea: Test stability of cross-sections, both across platforms and code
# releases. So fail on even epsilon change, but make it easy to plot the ones
# with the biggest changes and auto-update the reference files.

import numpy as np
from NCrystalDev._numpy import _np_geomspace, _np_linspace

import NCTestUtils.dirs as testdirs
import NCTestUtils.enable_fpe  # noqa F401


class XSMonitor:

    def __init__( self,
                  refdatadir,
                  matloadfct,
                  egridgenfct,
                  testlistgenfct,
                  test_rdtol = 1e-5,
                  reffile_format_evals = '%.7g',
                  reffile_format_xsvals = '%.7g' ):
        ( self.__do_plot,
          self.__do_update,
          self.__test_select ) = _parse_sysargv()
        self.__test_rdtol = test_rdtol
        self.__reffile_format_evals = reffile_format_evals
        self.__reffile_format_xsvals = reffile_format_xsvals
        self.__matloadfct = matloadfct
        self.__egridgenfct = egridgenfct
        self.__testlist = list( testlistgenfct() )
        if not self.__do_plot:
            from NCTestUtils.env import ncsetenv
            ncsetenv('FAKEPYPLOT','1')
        self.__refdir = _init_refdir( self.__do_update,
                                      refdatadir,
                                      self.__testlist )

    def reffile( self, teststr ):
        return self.__refdir.joinpath(_reffile_bn(teststr))

    def save_ref( self, teststr, evals, xsvals ):
        fmtstr = f'{self.__reffile_format_evals} {self.__reffile_format_xsvals}'
        o = '\n'.join(fmtstr%(e,x) for e,x in zip(evals,xsvals))
        o += '\n'
        f = self.reffile(teststr)
        f.write_text(o)
        print(f'Wrote: {f}')

    def load_ref( self, teststr ):
        f = self.reffile(teststr)
        data = np.loadtxt(f, dtype=float)
        assert data.ndim == 2 and data.shape[1] == 2
        evals, xsvals = data[:, 0], data[:, 1]
        assert len(evals) == len(xsvals)
        return evals.copy(), xsvals.copy()

    def run( self ):
        badtests = set()
        nused = 0
        for teststr in self.__testlist:
            testkey = testdirs.encode_safe(teststr)
            if self.__test_select and testkey not in self.__test_select:
                print(f"=============> SKIPPING test {testkey}")
                continue
            nused += 1
            print(f"=============> Launching test {testkey}")
            mat = self.__matloadfct(teststr)
            evals = self.__egridgenfct( mat )
            #collapse, sort, and reduce precision to that which will be stored
            #in reffiles:
            def unique_and_sort(a):
                a = set( a )
                a = np.asarray(list(a),dtype=float)
                a.sort()
                return a
            evals = [ float(self.__reffile_format_evals%v) for v in evals ]
            evals = unique_and_sort( evals )
            xsvals = mat.scatter.xsect(evals)
            if self.__do_update:
                self.save_ref(teststr,evals,xsvals)
                continue
            ref_evals, ref_xsvals = self.load_ref(teststr)
            if _reldiff(evals,ref_evals).max() > self.__test_rdtol:
                badtests.add((None,testkey))
                print(f"ERROR: reference e-grid for {testkey} is inconsistent."
                      " Developers: If expected, --update after investigating"
                      f" with: --plot {testkey} ")
            else:
                rda = _reldiff( xsvals, ref_xsvals )
                rd = rda.max()
                if rd>self.__test_rdtol:
                    badtests.add((rd,testkey))
                    print(f"ERROR: reference cross sections for {testkey} are"
                          f" inconsistent at the reldiff={rd:g} level."
                          " Developers: If expected, consider --update after"
                          f" investigating with: --plot {testkey} ")
                    print('Data dump:')
                    print('               Eref                  E            XSref               XS  XSreldiff')
                    def fmt(v, n=10):
                        return f"{v:.{n}g}".rjust(n+6)#+7 if can be negative
                    for i in range(len(xsvals)):
                        ref_estr = fmt(ref_evals[i],12)
                        estr = fmt(evals[i],12)
                        if estr==ref_estr:
                            estr = '            <same>'
                        rdval = rda[i]
                        rdstr = fmt(rdval,4)
                        if rdval>self.__test_rdtol:
                            rdstr+=' <-- problem'
                        print(f' {ref_estr}'
                              f' {estr}'
                              f' {fmt(ref_xsvals[i],10)}'
                              f' {fmt(xsvals[i],10)}'
                              f' {rdstr}')

            if self.__do_plot:
                evals_lux = set(evals)
                evals_lux |= set(_np_linspace(evals[0],evals[-1],5000))
                evals_lux |= set(_np_geomspace(evals[0],evals[-1],5000))
                evals_lux = unique_and_sort(evals_lux)
                import matplotlib.pyplot as plt
                fig, axs = plt.subplots(2, 1, sharex=True)
                fig.subplots_adjust(hspace=0)
                ax, axdiff = axs
                ax.set_title(teststr)
                ax.plot(evals_lux,mat.scatter.xsect(evals_lux),'-')
                ax.plot(evals,xsvals,'o',alpha=0.5,label='current')
                ax.plot(ref_evals,ref_xsvals,'d',alpha=0.5,label='ref')
                ax.legend().set_draggable(True)
                ax.loglog()
                ax.grid()
                ax.set_ylabel('XS [barn/atom]')
                axdiff.plot( evals, _reldiff( xsvals, ref_xsvals ),
                             label='observed difference')
                axdiff.loglog()
                axdiff.grid()
                axdiff.axhline(self.__test_rdtol,color='red',ls=':',
                               label='test tolerance')
                axdiff.set_ylabel('XS relative difference')
                axdiff.set_xlabel('Neutron energy [eV]')
                axdiff.legend().set_draggable(True)
                plt.show()

        if not nused:
            raise SystemExit('ERROR: No tests were run!')

        if self.__test_select and len(self.__test_select)!=nused:
            raise SystemExit('ERROR: Test selection selected one or more'
                             ' non-existent tests!')

        if not badtests:
            if not self.__do_update:
                print("ALL OK")
            return#all ok
        print()
        print("ERROR - ISSUES DETECTED: ")
        print()
        for rd, testkey in sorted(badtests):
            descr = ( f'xs rel diff: {rd:g}' if rd is not None
                      else 'inconsistent egrid' )
            print(f'  Issue with test ({descr}): {testkey}')
        print()
        print('Remember in a dev env you can rerun with --plot or --update. '
              'You can also specify one or more test names to only run those tests')
        raise SystemExit(1)

def _reffile_bn( teststr ):
    return testdirs.encode_safe(teststr)+'.txt'

def _parse_sysargv(args=None):
    import sys
    if args is None:
        args = sys.argv[1:]
    do_plot, do_update = False, False
    while '--plot' in args:
        args.remove('--plot')
        do_plot = True
    while '--update' in args:
        args.remove('--update')
        do_update = True
    return do_plot, do_update, set(args)

def _reldiff( a,b ):
    if a.size != b.size:
        return np.inf
    return np.abs(a - b) / (0.5*(np.abs(b)+np.abs(a))+np.finfo(float).eps)

def _init_refdir(do_update,dirname,testlist):
    refdir = testdirs.get_named_test_data_dir(dirname,
                                              for_updates = do_update)
    #Reminder to remove any unused files:
    present = {f.name for f in refdir.glob('*.txt')}
    expected = {_reffile_bn(e) for e in testlist}

    if present-expected:
        print('ERROR: Excess files found (please remove):')
        for e in present-expected:
            print(f'  {refdir.joinpath(e)}')
        raise SystemExit(1)
    if not do_update and expected-present:
        print('ERROR: Some reference files missing'
              ' (rerun with --update to generate):')
        cmd = '  --update'
        for e in expected-present:
            print(f'  {e}')
            assert e.endswith('.txt')
            cmd += ' ' + str(e)[:-4]
        print()
        print("To update JUST those missing, run with:")
        print()
        print(cmd)
        raise SystemExit(1)
    return refdir
