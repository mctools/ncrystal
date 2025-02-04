
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

class Cfg:

    @property
    def ncrystal_namespace(self):
        return 'dev'

    def sbpkgname_comp(self, compname):
        orig = """NCAbsFact NCBkgdExtCurve NCCore NCElIncScatter NCFactories
                NCFreeGas NCInterfaces NCPowderBragg NCScatFact NCDump NCAbsOOV
                NCCfgUtils NCData NCExperimental NCFactory_Laz NCGasMix
                NCLCBragg NCPubUtils NCSAB NCSCBragg NCUtils NCAtomDB
                NCCInterface NCDynInfoUtils NCExtdUtils NCFactory_NCMAT
                NCInfoBld NCMiniMC NCQuickFact NCSABScatter NCThreads NCVDOS"""
        guess = dict( (e[2:].lower(),e) for e in orig.split() )
        if compname not in guess:
            e = {'extd_utils':'NCExtdUtils',
                 'misc':'NCPubUtils',
                 'ncmat':'NCFactory_NCMAT',
                 'lazlau':'NCFactory_Laz',
                 'sanshardsphere':'NCExperimental',
                 'stdscatfactory':'NCScatFact',
                 'elasincoh':'NCElIncScatter'
                 }.get(compname)
            if e is None:
                e = 'NC' + ''.join(e.capitalize() for e in compname.split('_'))
            return e
        return guess[compname]

    def sbpkgname_testlib(self, testlibname):
        return f'TestLib_{testlibname}'#NB: no NC prefix

    @property
    def sbpkgname_lib(self):
        return self.sbpkgname_comp('cinterface')

    @property
    def sbpkgname_data(self):
        return 'NCData'

    @property
    def sbpkgname_pymods(self):
        return 'NCrystalDev'

    @property
    def sbpkgname_ncrystalhh(self):
        return 'NCrystalDev'

    @property
    def sbpkgname_cli(self):
        return 'NCCmd'

    @property
    def sbpkgname_examples(self):
        return 'NCExamples'

    @property
    def sbpkgname_geant4(self):
        return 'NCG4'

    @property
    def sbld_instdir(self):
        return self.__instdir

    @property
    def sbld_mode(self):
        return self.__bldmode

    @property
    def ncrystal_version_str(self):
        return self.__version_str

    @property
    def ncrystal_version_int(self):
        return self.__version_int

    @property
    def ncrystal_version_major(self):
        return self.__version_tuple[0]

    @property
    def ncrystal_version_minor(self):
        return self.__version_tuple[1]

    @property
    def ncrystal_version_patch(self):
        return self.__version_tuple[2]

    def __init__(self):
        import sys
        if len(sys.argv)<3:
            raise SystemExit('Error: expects cmdline arguments'
                             ': <sbldinstdir> <sbldmode>')
        import pathlib
        self.__instdir = pathlib.Path(sys.argv[1])
        self.__bldmode = sys.argv[2]

        from .dirs import reporoot
        self.__version_str = (reporoot/'VERSION').read_text().strip()
        self.__version_tuple = tuple( int(i)
                                      for i in self.__version_str.split('.') )
        self.__version_int = sum(int(i)*j
                                 for i,j
                                 in zip(self.__version_tuple,(1000000,1000,1)))


cfg = Cfg()
