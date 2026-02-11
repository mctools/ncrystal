
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

""" Obsolete interface. Use the NCrystal.minimc module instead."""

__all__ = ['quick_diffraction_pattern','runsim_diffraction_pattern']

def runsim_diffraction_pattern( *a, **kwargs ):
    """Obsolete function. Calling it now will result in an error."""
    #Fixme: include URL to new documentation page:
    from .exceptions import NCException
    raise NCException('The runsim_diffraction_pattern(..) function is'
                      ' obsolete. Please migrate your code to use the'
                      ' minimc_run(..) function instead')

__cache_qdpwarn=[True]
def quick_diffraction_pattern( cfgstr, *,
                               neutron_energy,
                               material_thickness,
                               nstat = 'auto',
                               nthreads = 'auto',
                               suppress_obsoletion_warning = False ):

    from .minimc import minimc_run_scenario
    migratemsg = ('Please migrate your code to use the'
                  ' NCrystal.minimc.minimc_run(..) or'
                  ' NCrystal.minimc.minimc_run_scenario(..) functions instead.')
    warnmsg = ('The quick_diffraction_pattern(..) function is obsolete.'
               f' {migratemsg}')

    if __cache_qdpwarn[0] and not suppress_obsoletion_warning:
        __cache_qdpwarn[0] = False
        from ._common import warn
        warn(warnmsg)

    assert ';' not in neutron_energy
    assert ';' not in material_thickness
    neutron_energy = ''.join(neutron_energy.split())
    material_thickness = ''.join(material_thickness.split())
    scenario_cfg = f'{neutron_energy} pencil on {material_thickness} sphere'
    enginecfg = f';nthreads={nthreads};tally=mu;tallybins=mu:1800:0:180'

    def simfct( n, cfgstr ):
        import time
        t0 = time.time()
        res = minimc_run_scenario( cfgstr,
                                   scenario_cfg + f' {n} times',
                                   extra_engineopts = enginecfg )
        t1 = time.time()
        return t1-t0, res

    if nstat is None or nstat=='auto':
        for nstat in [1e4,1e5,1e6,1e7]:
            t,res = simfct(nstat,cfgstr)
            #Usually, end within a second in total, but in worst cases, up to
            #10seconds:
            if t>0.1 and nstat >= 1e6:
                break
            if t>1.0:
                break
    else:
        t,res=simfct(nstat,cfgstr)

    class MiniMCObsoleteResult:
        """Backwards compatible results class for the obsolete
        quick_diffraction_pattern() function. The NCrystal.minimc.minimc_run(..)
        or NCrystal.minimc.minimc_run_scenario(..) functions should be used
        instead.
        """
        def __init__(self,res):
            self.__res = res

        def plot_breakdown(self,rebin_factor=1,logy=False):
            self.__res.tally('mu').plot(rebin_factor=rebin_factor,
                                        logy=logy)

        def __getattr__(self, name):
            if not name.startswith('_'):
                raise RuntimeError(f'"{name}" is not available on results from'
                                   ' the obsolete quick_diffraction_pattern'
                                   ' function. Only .plot_breakdown(..) is.'
                                   f' {migratemsg}')

    return MiniMCObsoleteResult(res)
