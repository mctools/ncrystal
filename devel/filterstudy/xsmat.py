
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

"""Helpers for studying NCrystal total cross section curves (for filter tables).

All cross sections are in barn per atom, wavelengths in Angstrom, energies in eV.
"""
import NCrystalDev as NC
import numpy as np

WL_MIN_DEFAULT = 0.01
WL_MAX_DEFAULT = 100.0

EXTRA_CFGS = [
    'gasmix::air',
    'gasmix::He3/10bar',
    'gasmix::BF3/2atm',
    'phases<0.3*stdlib::Al_sg225.ncmat&0.7*stdlib::Cu_sg225.ncmat>',
    'phases<0.9*stdlib::Al_sg225.ncmat&0.1*stdlib::Polyethylene_CH2.ncmat>',
    'stdlib::Be_sg194.ncmat;temp=80K',
    'stdlib::Be_sg194.ncmat;temp=20K',
    'stdlib::Al_sg225.ncmat;temp=600K',
    'stdlib::Al_sg225.ncmat;density=0.5x',
    'stdlib::Al_sg225.ncmat;dcutoff=0.1',
    'stdlib::Si_sg227.ncmat;dcutoff=0.2',
    'stdlib::Al_sg225.ncmat;vdoslux=1',
    'stdlib::Al_sg225.ncmat;inelas=0',
    'stdlib::Al_sg225.ncmat;coh_elas=0;incoh_elas=0;inelas=0',
]


def all_cfgs():
    return [f.fullKey for f in NC.browseFiles(factory='stdlib')] + EXTRA_CFGS


class Mat:
    """Total cross section (scattering+absorption) of an isotropic material."""

    def __init__(self, cfg):
        self.cfg = cfg
        self.sc = NC.createScatter(cfg)
        self.ab = NC.createAbsorption(cfg)
        self.info = NC.createInfo(cfg)
        if not self.sc.isNonOriented():
            raise ValueError('oriented material')

    def sigma(self, wl):
        wl = np.asarray(wl, dtype=float)
        return self.sc.xsect(wl=wl) + self.ab.xsect(wl=wl)

    def sigma_scat(self, wl):
        return self.sc.xsect(wl=np.asarray(wl, dtype=float))

    def sigma_abs(self, wl):
        return self.ab.xsect(wl=np.asarray(wl, dtype=float))

    @property
    def numberdensity(self):
        return self.info.numberdensity  # atoms/Aa^3

    def edges(self, wlmin=0.0, wlmax=np.inf):
        """Sorted unique Bragg edge wavelengths (2d) of all crystalline phases."""
        res = []

        def collect(info):
            if info.isMultiPhase():
                for _, ph in info.phases:
                    collect(ph)
            elif info.hasHKLInfo():
                for h in info.hklObjects():
                    res.append(2.0 * h.d)
        collect(self.info)
        e = np.unique(np.asarray(res, dtype=float))
        return e[(e >= wlmin) & (e <= wlmax)]

    def edge_jumps(self, edges, releps=1e-9):
        """Absolute jumps sigma(e-)-sigma(e+) at the given edges."""
        edges = np.asarray(edges, dtype=float)
        return self.sigma(edges * (1 - releps)) - self.sigma(edges * (1 + releps))


def dense_reference(mat, wlmin=WL_MIN_DEFAULT, wlmax=WL_MAX_DEFAULT,
                    n=200000, edge_pad=1e-9):
    """Dense reference sampling, including points just below and above each
    Bragg edge. Returns (wl, sigma), sorted in wl."""
    wl = np.geomspace(wlmin, wlmax, n)
    e = mat.edges(wlmin, wlmax)
    wl = np.unique(np.concatenate([wl, e * (1 - edge_pad), e * (1 + edge_pad)]))
    return wl, mat.sigma(wl)
