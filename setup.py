
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2023 NCrystal developers                                   ##
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

from skbuild import setup  # This line replaces 'from setuptools import setup'

def _extract_version():
    import pathlib
    p = pathlib.Path(__file__).parent / 'NCrystal' / '__init__.py'
    if not p.exists():
        raise RuntimeError("Unable to find NCrystal/__init__.py (for version string extraction).")
    version_str = None
    for l in p.read_text().splitlines():
        if l.startswith('__version__'):
            delim = '"' if '"' in l else "'"
            version_str = l.split(delim)[1]
            break
    if not version_str:
        raise RuntimeError(f'Unable to find version string in {p}.')
    v = tuple( int(i) for i in version_str.split('.') )
    if not len(v)==3:
        raise RuntimeError(f'Invalid version string extracted from {p}.')
    return v

def get_version():
    x,y,z = _extract_version()
    return f'{x}.{y}.dev{z}' if z>=80 else f'{x}.{y}.{z}'

setup(
    name="ncrystal",
    version=get_version(),
    author='NCrystal developers (Thomas Kittelmann, Xiao Xiao Cai)',
    license = "Apache-2.0",
    description='Library for thermal neutron transport in crystals and other materials.',
    url='https://github.com/mctools/ncrystal/wiki',
    keywords='neutron,montecarlo,science',
    packages=['NCrystal'],
    python_requires='>=3.6, <4',
    install_requires=['numpy'],
    long_description='NCrystal is a library and associated tools which enables calculations for Monte Carlo simulations of thermal neutrons in crystals and other materials, supporting a range of physics including both coherent, incoherent, elastic and inelastic scatterings in a wide range of materials, including crystal powders, mosaic single crystals, layered single crystals, amorphous solids, and liquids. Multiphase materials or isotopically enriched materials are supported as well, and the framework furthermore supports phase-contrast (SANS) physics.',
    cmake_args=['-DNCRYSTAL_NOTOUCH_CMAKE_BUILD_TYPE=ON',
                '-DNCRYSTAL_MODIFY_RPATH=OFF',
                '-DNCRYSTAL_ENABLE_SETUPSH=OFF',
                '-DNCRYSTAL_ENABLE_DATA=EMBED',#the c++ library can not currently locate the files if not embedding
                '-DNCRYSTAL_ENABLE_MCSTAS=OFF',#Explicitly disable for now (was already moved into upstream McStas)
                '-DNCRYSTAL_ENABLE_GEANT4=OFF',#Explicitly disable for now (planning to move out of core NCrystal repo)
                '-DNCRYSTAL_SKIP_PYMODINST=ON',
                '-DCMAKE_INSTALL_LIBDIR=NCrystal/ncrystal_pyinst_data/lib',
                '-DCMAKE_INSTALL_INCLUDEDIR=NCrystal/ncrystal_pyinst_data/include',
                '-DCMAKE_INSTALL_DATADIR=NCrystal/ncrystal_pyinst_data/data',
                '-DNCrystal_DATAFILESDIR=NCrystal/ncrystal_pyinst_data/stdlib_data',
                ]
)
