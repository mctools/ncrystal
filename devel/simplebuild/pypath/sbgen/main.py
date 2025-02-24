
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

from .file import add_file, create_files
from . import dirs
from .cfg import cfg
import os
import itertools

_pydep2sblddep = { 'ase':'ASE',
                   'gemmi':'Gemmi',
                   'numpy':'Numpy',
                   'pandas':'Pandas',
                   'toml':'PyToml',
                   'scipy':'Scipy',
                   'spglib':'Spglib',
                   'matplotlib':'matplotlib',
                   'mpmath':'mpmath' }

_pydeps2pkg_suffix = [ ( set(['numpy']), 'np' ),
                       ( set(['toml']), 'toml' ),
                       ( set(['numpy','matplotlib']), 'mpl' ),
                       ( set(['numpy','mpmath']), 'mpmath' ),
                       ( set(['numpy','ase','spglib','gemmi']), 'asg' )
                      ]

def determine_testpkg_by_pydeps( pydeps ):
    if not pydeps:
        return 'NCTestPy', set()
    #find shortest _pydeps2pkg_suffix which has all deps:
    bestset = None
    bestname = None
    for kset, kname in _pydeps2pkg_suffix:
        if pydeps <= kset and ( bestset is None or len(bestset)>len(kset)):
            bestset = kset
            bestname = kname
    if bestset is None:
        return None, pydeps#not possible, needs custom pkg later
    return 'NCTestPy%s'%bestname, bestset

def create_pkg_extratestscripts(pkgname):
    import pathlib
    create_pkginfo( pkgname, pkg_deps = ['NCrystalDev','NCCmd'] )
    for sf in sorted(pathlib.Path(__file__).parent.glob('_testscript_*.py')):
        n = sf.stem[len('_testscript_'):]
        add_file( f'pkgs/{pkgname}/scripts/test{n}', link_target = sf )

def create_pkginfo_pytestpkg( pkgname,
                              pydeps ):
    create_pkginfo( pkgname,
                    extdeps = [ get_dependency_from_pydep(e)
                                for e in (pydeps or [])],
                    pkg_deps = [cfg.sbpkgname_pymods] )

_created_custom_pydeps = set()
def get_dependency_from_pydep( pd):
    s = _pydep2sblddep.get(pd)
    if s is not None:
        return s
    #Need custom pydep
    name = 'Py%s'%pd.capitalize()
    if name in _created_custom_pydeps:
        return name#already created
    _created_custom_pydeps.add(name)
    version_extractor = 'm.__version__'

    content = """
execute_process(
    COMMAND "${Python_EXECUTABLE}" "-c"
    "import @pd@ as m; print(@ve@)"
    OUTPUT_VARIABLE tmp RESULT_VARIABLE tmp_ec
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
if ("x${tmp_ec}" STREQUAL "x0")
  set(HAS_@name@ 1)
  string(STRIP "${tmp}" ExtDep_@name@_VERSION)
  set(ExtDep_@name@_PREPEND_COMPILE_FLAGS "")
  set(ExtDep_@name@_PREPEND_LINK_FLAGS "")
else()
  set(HAS_@name@ 0)
endif()
""".replace('@name@',name).replace('@pd@',pd).replace('@ve@',version_extractor )
    add_file( dirs.genroot.joinpath('extdep_definitions',
                                    f'ExtDep_{name}.cmake'),
              content = content )
    return name

def extract_deps_from_needs(sf):
    with sf.open('rt') as fh:
        for line in fh:
            if line.startswith('# NEEDS: '):
                return set( line[9:].split() )

def create_pkginfo( pkgname,
                    extdeps = None,
                    extra_cflags = None,
                    pkg_deps = None ):
    c = 'USEPKG ' if pkg_deps else ''
    if pkg_deps:
        c+= ' '.join(sorted(pkg_deps))
    if extdeps:
        c += ' USEEXT '
        c+= ' '.join(sorted(extdeps))
    if extra_cflags:
        c += ' EXTRA_COMPILE_FLAGS '
        c += ' '.join(extra_cflags)
    add_file( f'pkgs/{pkgname}/pkg.info', content = f'package({c})\n' )

def ncapi_contents():
    c = dirs.srcroot.joinpath('include/NCrystal/ncapi.h.in').read_text()
    #NOTE: NCRYSTAL_NAMESPACED_ENVVARS since we also add
    #_do_namespace_envvars.py below!!
    c = c.replace('/* @NCRYSTAL_HOOK_FOR_ADDING_DEFINES@ */',
f"""
#define NCRYSTAL_NAMESPACE_PROTECTION {cfg.ncrystal_namespace}
#define NCRYSTAL_DISABLE_CMDLINEPLUGINMGR
#define NCRYSTAL_SIMPLEBUILD_DEVEL_MODE
#define NCRYSTAL_NO_CMATH_CONSTANTS
#define NCRYSTAL_NAMESPACED_ENVVARS
#define NCRYSTAL_VERSION_MAJOR {cfg.ncrystal_version_major}
#define NCRYSTAL_VERSION_MINOR {cfg.ncrystal_version_minor}
#define NCRYSTAL_VERSION_PATCH {cfg.ncrystal_version_patch}
#define NCRYSTAL_VERSION_STR "{cfg.ncrystal_version_str}"
#define NCRYSTAL_VERSION {cfg.ncrystal_version_int}
//Uncomment to debug deadlocks: #define NCRYSTAL_DEBUG_LOCKS
#ifndef NDEBUG
// Specifically test aligned allocs in simplebuild debug mode:
#  define NCRYSTAL_TRACKALIGNEDALLOC
#endif

""")
    return c

def create_custom_extdep( name, cflags = '', ldflags = '',
                          is_not_syshdrs = True,
                          includemap = None ):
    outdir = dirs.genroot.joinpath('extdep_definitions')

    content=f"""
set( HAS_{name} "1" )
set( ExtDep_{name}_VERSION "{cfg.ncrystal_version_str}" )
set( ExtDep_{name}_COMPILE_FLAGS "{cflags}")
set( ExtDep_{name}_LINK_FLAGS "{ldflags}")

"""
    #Make sure we intercept headers from system-installed NCrystal (default
    #priority is 1000):
    content += (f'set( ExtDep_{name}_FLAGPRIORITY "100" )\n')

    #Map includes like "NCrystal/core/NCBla.hh" to "NCCore/NCBla.hh"
    if includemap:
        incmapfilebn = f'extdep_incmap_{name}.txt'
        content += (f'set( ExtDep_{name}_INCLUDEMAPFILE '
                    '"${CMAKE_CURRENT_LIST_DIR}/%s" )\n'%incmapfilebn )
        incmap_content = ''
        for k,v in includemap:
            incmap_content += f'{k} {v}\n'
        add_file( outdir.joinpath(incmapfilebn), content = incmap_content )

    if is_not_syshdrs:
        #Do not add -isystem to our include paths:
        content += 'set( ExtDep_{name}_IsNotSystemHeaders "ON" )\n'

    add_file( outdir / f'ExtDep_{name}.cmake', content = content )

def platform_so_ending():
    import platform
    return '.dylib' if platform.system() == 'Darwin' else '.so'

def pkg_libname( pkgname ):
    libname = 'libPKG__%s'%pkgname
    return libname + platform_so_ending()

def ncconfig_h_contents():

    bin2incdir = os.path.relpath( dirs.srcroot.joinpath('include'),
                                  cfg.sbld_instdir.joinpath('bin') )

    cmakebuildtype = dict( debug = 'Debug',
                           reldbg = 'RelWithDebInfo',
                           release = 'Release' )[ cfg.sbld_mode ]


    libname = pkg_libname(cfg.sbpkgname_lib)
    expandvars = dict( NCLIBNAME = libname,
                       NCDATAPKGNAME = cfg.sbpkgname_data,
                       NCVERSION = cfg.ncrystal_version_str,
                       NCINTVERSION = str(cfg.ncrystal_version_int),
                       NCNAMESPACE = cfg.ncrystal_namespace,
                       NCBIN2INCDIR = bin2incdir,
                       CMAKEBUILDTYPE = cmakebuildtype,
                      )
    c = """
const char * nccfg_const_bin2libdir(void) { return "../lib"; }
const char * nccfg_const_bin2shlibdir(void) { return "../lib"; }
const char * nccfg_const_libname(void) { return "@NCLIBNAME@"; }
const char * nccfg_const_shlibname(void) { return "@NCLIBNAME@"; }
const char * nccfg_const_bin2libpath(void)
    { return "../lib/@NCLIBNAME@"; }
const char * nccfg_const_bin2shlibpath(void)
    { return "../lib/@NCLIBNAME@"; }
const char * nccfg_const_bin2datadir(void) { return "../data/@NCDATAPKGNAME@"; }
const char * nccfg_const_bin2incdir(void) { return "@NCBIN2INCDIR@"; }
const char * nccfg_const_bin2cmakedir(void)
    { return "cmakedir/not/available/in/simplebuild/devel/mode"; }
const char * nccfg_const_version(void) { return "@NCVERSION@"; }
const char * nccfg_const_intversion(void) { return "@NCINTVERSION@"; }
const char * nccfg_const_namespace(void) { return "@NCNAMESPACE@"; }
const char * nccfg_const_cmakebuildtype(void) { return "@CMAKEBUILDTYPE@"; }
int nccfg_boolopt_data(void) { return 1; }
int nccfg_boolopt_dynamic_plugins(void) { return 1; }
int nccfg_boolopt_embed_data(void) { return 0; }
int nccfg_boolopt_examples(void) { return 1; }
int nccfg_boolopt_modify_rpath(void) { return 1; }
int nccfg_boolopt_threads(void) { return 1; }
int nccfg_boolopt_expects_shlibdir_override(void) { return 0; }
"""
    for k,v in expandvars.items():
        c = c.replace('@%s@'%k,v)
    return c

def create_testplugin_pkg(pkgname,pkg_deps):
    create_pkginfo( pkgname,
                    pkg_deps = pkg_deps,
                    extra_cflags = ['-DNCPLUGIN_NAME=DummyPlugin'] )

    p = dirs.exsrcroot.joinpath('plugin')
    for sf in itertools.chain( p.joinpath('src').glob('*.hh'),
                               p.joinpath('src').glob('*.cc') ):
        add_file( f'pkgs/{pkgname}/libsrc/{sf.name}', link_target = sf )
    for sf in p.joinpath('data').glob('*.ncmat'):
        add_file( f'pkgs/{pkgname}/data/{sf.name}', link_target = sf )


def define_files():

    #TODO: Add an ncdevenv wrapper, which make sure that all ncrystal cmdline
    #thng have the normal names (and intercepts anything installed in the
    #system). Then we can do "ncdevenv bash" or "ncdevenv mcrun" and things like
    #that, in order to test in an environment with proper names. We should in
    #that case also make sure that python imports work with the module name
    #NCrystal and not just NCrystalDev.

    #TODO: Also autogenerate NCrystal.hh, adding all public includes? Problem
    #is, we need both python and cmake code for it. Might be easier with a
    #static one + a standalone test that verifies we did not forget anything.

    from . import dirs

    #NCrystal headers:
    ncapi_loc_in_genroot = 'incpath/NCrystal/ncapi.h'
    add_file( ncapi_loc_in_genroot, content=ncapi_contents() )

    #NCrystal source file compilation:
    from ncrystal_repo_tools.core_components import load_components
    name2comp = load_components()
    all_ncsbpkgs = []
    includemap = []
    for name,comp in sorted(name2comp.items()):
        sbpkgname = cfg.sbpkgname_comp(name)
        all_ncsbpkgs.append(sbpkgname)
        extdeps = ['NCDevHeaders']
        if name=='factories':
            extdeps.append('DL')
        if comp.name=='minimc':
            extdeps.append('Threads')
        tmp_extra_flags = [f'-I{dirs.srcroot}/src']
        if comp.name=='factories':
            tmp_extra_flags.append( f'-DNCRYSTAL_DATADIR={dirs.datadir}' )
        create_pkginfo( sbpkgname,
                        extdeps = extdeps,
                        pkg_deps = [ cfg.sbpkgname_comp(n)
                                     for n in comp.direct_depnames ],
                        extra_cflags = tmp_extra_flags )
        for sf in (comp.local_hdrs+comp.srcfiles):
            add_file( f'pkgs/{sbpkgname}/libsrc/{sf.name}', link_target = sf )

        #For traditional simplebuild includes:
        for sf in comp.hdrfiles or []:
            includemap.append( ( str(sf.relative_to(dirs.srcincroot)), f'{sbpkgname}/{sf.name}' ) )
            add_file( f'pkgs/{sbpkgname}/libinc/{sf.name}', link_target = sf )

    #Also symlink ncapi.h and NCrystal.hh:
    pkgname_core = cfg.sbpkgname_comp('core')
    add_file( f'pkgs/{pkgname_core}/libinc/ncapi.h',
              link_target = dirs.genroot.joinpath(ncapi_loc_in_genroot) )
    includemap.append( ( 'NCrystal/ncapi.h', f'{pkgname_core}/ncapi.h' ) )

    for fn in ['NCrystal.hh','NCPluginBoilerplate.hh']:
        add_file( f'pkgs/{cfg.sbpkgname_ncrystalhh}/libinc/{fn}',
                  link_target = dirs.srcroot.joinpath(f'include/NCrystal/{fn}') )

    #NCrystalDev package (python module):
    assert cfg.sbpkgname_lib in all_ncsbpkgs
    create_pkginfo( 'NCrystalDev', pkg_deps=all_ncsbpkgs + ['NCData'] )
    for sf in (dirs.pysrcroot/'src/NCrystal').glob('*.py'):
        add_file( f'pkgs/NCrystalDev/python/{sf.name}', link_target = sf )
    #Special marker used by _locatelib.py:
    add_file( 'pkgs/NCrystalDev/python/_is_sblddevel.py', content='' )

    #We added '#define NCRYSTAL_NAMESPACED_ENVVARS' in ncapi.h, so we also let
    #the python layer know:
    add_file( 'pkgs/NCrystalDev/python/_do_namespace_envvars.py', content='' )

    #Commandline scripts:
    create_pkginfo( cfg.sbpkgname_cli, pkg_deps=['NCrystalDev'])
    for sf in (dirs.pysrcroot/'src/NCrystal').glob('_cli_*.py'):
        cliname = sf.name[len('_cli_'):-len('.py')]
        sbscriptname = 'tool' if cliname == 'nctool' else cliname
        content = f"""#!/usr/bin/env python3
import NCrystalDev.{sf.stem} as mod
mod.main()
"""
        add_file( f'pkgs/{cfg.sbpkgname_cli}/scripts/{sbscriptname}',
                  content=content,
                  make_executable = True)
    for subpath,fn in [('include/NCrystal/internal/utils/NCCFileUtils.hh',
                        'NCCFileUtils.h'),
                       ('src/utils/NCCFileUtils.cc','NCCFileUtils.c'),
                       ('app_config/main.c','main.c')]:
        sf = dirs.srcroot / subpath
        add_file( f'pkgs/{cfg.sbpkgname_cli}/app_config/{fn}',
                  link_target = sf )
    add_file( f'pkgs/{cfg.sbpkgname_cli}/app_config/ncconfig_autogen.h',
              content = ncconfig_h_contents() )

    #Data files (even though NCLib already knows the path via a define):
    create_pkginfo( cfg.sbpkgname_data )
    for sf in dirs.datadir.glob('*.ncmat'):
        add_file( f'pkgs/{cfg.sbpkgname_data}/data/{sf.name}',
                  link_target = sf )

    #TODO: One custom extdep per pkg, and in an extra incpath, so we enforce the
    #deps to be correct??.
    create_custom_extdep( 'NCDevHeaders',
                          includemap = includemap,
                          cflags = (f'-I{dirs.genroot}/incpath '
                                    f'-I{dirs.srcroot}/include '
                                    '-DNCrystal_EXPORTS '#for NCCFileUtils.hh
                                    '-fno-math-errno') )
    #Examples:
    example_src =sorted(dirs.exsrcroot.glob('*.c*'))

    def add_compiled_examples(filelist,pkgname,*,override_name = None):
        if override_name:
            assert len(filelist)==1
        for f in filelist:
            assert f.name.startswith('ncrystal_example_')
            bn, ext = f.name[len('ncrystal_example_'):].split('.')
            bn = override_name or bn
            add_file( f'pkgs/{pkgname}/app_{bn}/main.{ext}',
                      link_target = f )
            app_name = f'sb_{pkgname}_{bn}'.lower()
            #also add test (no ref log):
            add_file( f'pkgs/{pkgname}/scripts/test{bn}',
                      content=f'#!/usr/bin/env bash\n{app_name}\n',
                      make_executable=True )

    create_pkginfo( cfg.sbpkgname_examples,
                    pkg_deps = [cfg.sbpkgname_lib] )
    add_compiled_examples(example_src,cfg.sbpkgname_examples)

    #Test libs:
    all_test_pkgs = []
    #all_test_pkgs_needing_ncrystal = []
    testlib_pkgnames = []
    for testlibdir in dirs.testroot.joinpath('libs').glob('lib_*'):
        tlname = testlibdir.name[4:]
        pkgname = cfg.sbpkgname_testlib(tlname)
        testlib_pkgnames.append(pkgname)
        incdir = testlibdir.joinpath('include',f'TestLib_{tlname}')
        all_test_pkgs.append(pkgname)
        create_pkginfo(pkgname,pkg_deps=[cfg.sbpkgname_ncrystalhh])
        sfs_cpp = list(testlibdir.glob('*.cc'))
        sfs_c = list(testlibdir.glob('*.c'))
        sfs_priv_h = list(testlibdir.glob('*.hh' if sfs_cpp else '*.h'))
        if sfs_cpp and sfs_c:
            raise RuntimeError("can not mix languages in sbld pkgs"
                               f" (problems in {testlibdir})")
        for sf in itertools.chain((sfs_cpp or sfs_c),sfs_priv_h):
            add_file( f'pkgs/{pkgname}/libsrc/{sf.name}',
                      link_target = sf )

        for sf in itertools.chain(incdir.glob('*.h'),incdir.glob('*.hh')):
            add_file( f'pkgs/{pkgname}/libinc/{sf.name}',
                      link_target = sf )

    #Test apps:
    all_test_pkgs.append('NCTestApps')
    create_pkginfo('NCTestApps',pkg_deps=testlib_pkgnames)
    for testappdir in dirs.testroot.joinpath('src').glob('app_*'):
        appname = testappdir.name[4:]
        for sf in itertools.chain( testappdir.glob('*.h'),
                                   testappdir.glob('*.hh'),
                                   testappdir.glob('*.c'),
                                   testappdir.glob('*.cc'),
                                   testappdir.glob('test.log') ):
            add_file( f'pkgs/NCTestApps/app_test{appname}/{sf.name}',
                      link_target = sf )

    #Test modules (shared libs for loading via python ctypes):
    testmods_pkgnames = []
    for moddir in dirs.testroot.joinpath('modules').glob('lib_*'):
        tlname = moddir.name[4:]
        pkgname = f'NCTestMod_{tlname}'
        testmods_pkgnames.append(pkgname)
        all_test_pkgs.append(pkgname)
        create_pkginfo(pkgname,
                       pkg_deps=[cfg.sbpkgname_ncrystalhh,
                                 'NCTestUtils'])
        for sf in itertools.chain( moddir.glob('*.h'),
                                   moddir.glob('*.hh'),
                                   moddir.glob('*.c'),
                                   moddir.glob('*.cc') ):
            add_file( f'pkgs/{pkgname}/libsrc/{sf.name}',
                      link_target = sf )

    #Test plugin module:
    all_test_pkgs.append( 'NCTestPlugin' )
    create_testplugin_pkg( 'NCTestPlugin',
                           pkg_deps = [ cfg.sbpkgname_ncrystalhh ] )

    ##Exercise the CMake code versus an NCrystal installed into the test env:
    #all_test_pkgs_needing_ncrystal.append( 'NCTestPlugin2' )
    #create_testplugin_pkg( 'NCTestPlugin2',
    #                       pkg_deps = [ cfg.sbpkgname_ncrystalhh ] )

    #Python modules for test scripts, common headers, and data:
    pkgname = 'NCTestUtils'
    all_test_pkgs.append(pkgname)
    for sf in dirs.testroot.joinpath('pypath/NCTestUtils').glob('*.py'):
        add_file( f'pkgs/{pkgname}/python/{sf.name}',
                  link_target = sf )
    for sf in itertools.chain( dirs.testroot.joinpath('include/NCTestUtils').glob('*.h'),
                               dirs.testroot.joinpath('include/NCTestUtils').glob('*.hh') ):
        add_file( f'pkgs/{pkgname}/libinc/{sf.name}',
                  link_target = sf )
    #Add special marker so it knows it is in sbld mode:
    add_file( f'pkgs/{pkgname}/python/_is_simplebuild.py', content='')
    extra_data_pkgs = []
    for sf in dirs.testroot.joinpath('data').iterdir():
        if sf.is_file() and '#' not in sf.name and '~' not in sf.name:
            add_file( f'pkgs/{pkgname}/data/{sf.name}', link_target = sf )
        if sf.is_dir() and not dirs.is_empty_dir( sf ):
            pkgname_edp = f'NCTestData_{sf.name}'
            extra_data_pkgs.append( (pkgname_edp, sf)  )

    create_pkginfo( pkgname,
                    pkg_deps = [pn for pn,_ in extra_data_pkgs] )

    for pkgname, srcdir in extra_data_pkgs:
        create_pkginfo( pkgname )
        for sf in srcdir.iterdir():
            if sf.is_file() and '#' not in sf.name and '~' not in sf.name:
                add_file( f'pkgs/{pkgname}/data/{sf.name}', link_target = sf )

    #Finally the test scripts. Here we must parse the files and look for NEEDS
    #lines to figure out which ones to include:
    #depset_2_testpyscripts
    pytest_pkgname2pydeps = {}
    extrapkg_pydeps = set()
    extrapkg_sflist = set()
    for sf in dirs.testroot.joinpath('scripts').glob('*.py'):
        if '#' in sf.name or '~' in sf.name:
            continue
        pydeps = extract_deps_from_needs(sf)
        pkgname, pkgpydeps = determine_testpkg_by_pydeps( pydeps )
        if pkgname is None:
            extrapkg_pydeps.update(pkgpydeps)
            extrapkg_sflist.add(sf)
            continue
        if pkgname not in pytest_pkgname2pydeps:
            pytest_pkgname2pydeps[pkgname] = pkgpydeps
        else:
            assert pytest_pkgname2pydeps[pkgname] == pkgpydeps
        outfn = f'pkgs/{pkgname}/scripts/test{sf.stem}'
        add_file( outfn, link_target = sf, make_executable = True )
        reflog = sf.parent.joinpath(sf.stem+'.log')
        if reflog.is_file():
            add_file( outfn + '.log', link_target = reflog )

    for pkgname, pydeps in pytest_pkgname2pydeps.items():
        all_test_pkgs.append(pkgname)
        create_pkginfo_pytestpkg( pkgname, pydeps )

    #spillover:
    if extrapkg_sflist:
        assert extrapkg_pydeps
        all_test_pkgs.append('NCTestPyExtra')
        create_pkginfo_pytestpkg( 'NCTestPyExtra', extrapkg_pydeps )
        for sf in extrapkg_sflist:
            outfn = f'pkgs/NCTestPyExtra/scripts/test{sf.stem}'
            add_file( outfn, link_target = sf, make_executable = True )
            reflog = sf.parent.joinpath(sf.stem+'.log')
            if reflog.is_file():
                add_file( outfn + '.log', link_target = reflog )

    #Link local _testscripts_*.py as test scripts:
    all_test_pkgs.append('NCExtraTestScripts')
    create_pkg_extratestscripts('NCExtraTestScripts')

    #For --require:
    create_pkginfo( 'NCTestAll', pkg_deps = all_test_pkgs )

def main():
    define_files()
    create_files()

#TODO: Make globbing utils for multiple extensions + guard against names with '#' or '~'
