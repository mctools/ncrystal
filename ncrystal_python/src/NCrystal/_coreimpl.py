
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

"""

Internal implementation details for functionality in core.py

"""

__all__ = []

class divdos_methods:

    """Common methods shared between DI_VDOS and DI_VDOSDebye"""

    @staticmethod
    def _plot_Gn():
        def plot_Gn( _self, n, nmax=None, without_xsect = False, **plotkwargs ):
            """Calculate and plot Sjolander's Gn function of n'th order. If
            without_xsect is True, the curves do not include a factor of the
            bound scattering cross section. Internally, the plots are produced
            with NCrystal.plot.plot_vdos_Gn, and any plotkwargs are passed along
            to that function
            """
            assert 1 <= n <= 99999 and ( nmax is None or n<=nmax<=99999 )
            n2 = int( nmax if nmax is not None else n )
            from .vdos import extractGn
            from .plot import plot_vdos_Gn
            def f(_n):
                return extractGn( n=_n, vdos = _self, mass_amu = _self.atomData.averageMassAMU(),
                                  temperature = _self.temperature, expand_egrid = True,
                                  scatxs = 1.0 if without_xsect else _self.atomData.scatteringXS() )
            plot_vdos_Gn( [ (list(f(i))+[i]) for i in range(int(n),n2+1) ], **plotkwargs )
        return plot_Gn

    @staticmethod
    def _extract_Gn():
        def extract_Gn( _self, n, *, expand_egrid = True, without_xsect = False ):
            """Calculate Sjolander's Gn function of n'th order. If expand_egrid
            returns (egrid,Gn_values), otherwise ((emin,emax),Gn_values). If
            without_xsect is True, the result will not be multiplied by the
            bound scattering cross section.
            """
            assert 1 <= n <= 99999
            from .vdos import extractGn as _extgn
            return _extgn( _self, n=int(n),
                           mass_amu=_self.atomData.averageMassAMU(),
                           temperature=_self.temperature,
                           scatxs = 1.0 if without_xsect else _self.atomData.scatteringXS(),
                           expand_egrid = expand_egrid )
        return extract_Gn

    @staticmethod
    def _extract_custom_knl():
        def extract_custom_knl( _self,
                                vdoslux = 3,
                                order_weight_fct = None,
                                without_xsect = False,
                                plot = False,
                                **plotkwargs ):
            """Extract (and optionally plot) a scattering kernel based on the
            VDOS curve. See the NCrystal.vdos.extractKnl function for a
            description of the arguments, the only difference being that the
            mass_amu, temperature, and scatxs parameters do not need to be
            provided here, and the without_xsects parameter can be used to avoid
            having the kernel include a factor of the bound scattering cross
            section.
            """
            from .vdos import extractKnl
            return extractKnl( _self,
                               vdoslux = vdoslux,
                               mass_amu = _self.atomData.averageMassAMU(),
                               temperature = _self.temperature,
                               scatxs = 1.0 if without_xsect else _self.atomData.scatteringXS(),
                               order_weight_fct = order_weight_fct,
                               plot = plot, **plotkwargs )

        return extract_custom_knl

    @staticmethod
    def _plot_vdos():
        def plot_vdos( _self, **kwargs ):
            """Plot VDOS using the NCrystal.plot.plot_vdos function. Any
            kwargs are simply passed along."""
            from .plot import plot_vdos
            plot_vdos( _self, **kwargs )
        return plot_vdos

    @staticmethod
    def _plot_knl():
        def plot_knl( self, vdoslux = 3, **kwargs ):
            """Plot the scattering kernel using the NCrystal.plot.plot_knl
               function. Any kwargs are simply passed along."""
            from .plot import plot_knl
            plot_knl( self._loadKernel( vdoslux=vdoslux ), **( kwargs or {}) )
        return plot_knl
