
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
import codecs
import os.path

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

setup(
    name="ncrystal",
    version=get_version("NCrystal/__init__.py"),
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
                '-DNCRYSTAL_ENABLE_DATA=EMBED',
                '-DNCRYSTAL_ENABLE_MCSTAS=OFF',#Explicitly disable for now (was already moved into upstream McStas)
                '-DNCRYSTAL_ENABLE_GEANT4=OFF',#Explicitly disable for now (planning to move out of core NCrystal repo)
                '-DNCRYSTAL_SKIP_PYMODINST=ON',
                '-DCMAKE_INSTALL_LIBDIR=NCrystal/ncrystal_pyinst_data/lib',
                '-DCMAKE_INSTALL_INCLUDEDIR=NCrystal/ncrystal_pyinst_data/include',
                '-DCMAKE_INSTALL_DATADIR=NCrystal/ncrystal_pyinst_data/data',
                '-DNCrystal_DATAFILESDIR=NCrystal/ncrystal_pyinst_data/stdlib_data',
    ]
)
