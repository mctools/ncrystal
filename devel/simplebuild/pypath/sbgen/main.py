
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

from .file import add_file, create_files
from . import dirs
from .cfg import cfg
import os

def create_pkginfo( pkgname,
                    extdeps = None,
                    extra_cflags = None,
                    pkg_deps = None ):
    c = 'USEPKG ' if pkg_deps else ''
    if pkg_deps:
        c+= ' '.join(pkg_deps)
    if extdeps:
        c += ' USEEXT '
        c+= ' '.join(extdeps)
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
#define NCRYSTAL_SIMPLEBUILD_DEVEL_MODE
#define NCRYSTAL_NO_CMATH_CONSTANTS
#define NCRYSTAL_NAMESPACED_ENVVARS
#define NCRYSTAL_VERSION_MAJOR {cfg.ncrystal_version_major}
#define NCRYSTAL_VERSION_MINOR {cfg.ncrystal_version_minor}
#define NCRYSTAL_VERSION_PATCH {cfg.ncrystal_version_patch}
#define NCRYSTAL_VERSION_STR "{cfg.ncrystal_version_str}"
#define NCRYSTAL_VERSION {cfg.ncrystal_version_int}
""")
    return c

def create_custom_extdep( name, cflags = '', ldflags = '' ):
    content=f"""
set( HAS_{name} "1" )
set( ExtDep_{name}_VERSION "{cfg.ncrystal_version_str}" )
set( ExtDep_{name}_COMPILE_FLAGS "{cflags}")
set( ExtDep_{name}_LINK_FLAGS "{ldflags}")
"""
    add_file( dirs.genroot.joinpath('extdep_definitions',
                                    f'ExtDep_{name}.cmake'),
              content = content )

def ncconfig_h_contents():

    bin2incdir = os.path.relpath( dirs.srcroot.joinpath('include'),
                                  cfg.sbld_instdir.joinpath('bin') )

    cmakebuildtype = dict( debug = 'Debug',
                           reldbg = 'RelWithDebInfo',
                           release = 'Release' )[ cfg.sbld_mode ]

    expandvars = dict( NCLIBPKGNAME = cfg.sbpkgname_ncrystal_lib,
                       NCDATAPKGNAME = cfg.sbpkgname_ncrystal_data,
                       NCVERSION = cfg.ncrystal_version_str,
                       NCINTVERSION = str(cfg.ncrystal_version_int),
                       NCNAMESPACE = cfg.ncrystal_namespace,
                       NCBIN2INCDIR = bin2incdir,
                       CMAKEBUILDTYPE = cmakebuildtype,
                      )
    c = """
const char * nccfg_const_bin2libdir() { return "../lib"; }
const char * nccfg_const_bin2shlibdir() { return "../lib"; }
const char * nccfg_const_libname() { return "libPKG__@NCLIBPKGNAME@.so"; }
const char * nccfg_const_shlibname() { return "libPKG__@NCLIBPKGNAME@.so"; }
const char * nccfg_const_bin2libpath()
    { return "../lib/libPKG__@NCLIBPKGNAME@.so"; }
const char * nccfg_const_bin2shlibpath()
    { return "../lib/libPKG__@NCLIBPKGNAME@.so"; }
const char * nccfg_const_bin2datadir() { return "../data/@NCDATAPKGNAME@"; }
const char * nccfg_const_bin2incdir() { return "@NCBIN2INCDIR@"; }
const char * nccfg_const_bin2cmakedir()
    { return "cmakedir/not/available/in/simplebuild/devel/mode"; }
const char * nccfg_const_version() { return "@NCVERSION@"; }
const char * nccfg_const_intversion() { return "@NCINTVERSION@"; }
const char * nccfg_const_builtinplugins() { return ""; }
const char * nccfg_const_namespace() { return "@NCNAMESPACE@"; }
const char * nccfg_const_cmakebuildtype() { return "@CMAKEBUILDTYPE@"; }
int nccfg_boolopt_data() { return 1; }
int nccfg_boolopt_dynamic_plugins() { return 1; }
int nccfg_boolopt_embed_data() { return 0; }
int nccfg_boolopt_examples() { return 1; }
int nccfg_boolopt_modify_rpath() { return 1; }
int nccfg_boolopt_threads() { return 1; }
int nccfg_boolopt_expects_shlibdir_override() { return 0; }
"""
    for k,v in expandvars.items():
        c = c.replace('@%s@'%k,v)
    return c

def define_files():

    #TODO: Add an ncdevenv wrapper, which make sure that all ncrystal cmdline
    #thng have the normal names (and intercepts anything installed in the
    #system). Then we can do "ncdevenv bash" or "ncdevenv mcrun" and things like
    #that, in order to test in an environment with proper names. We should in
    #that case also make sure that python imports work with the module name
    #NCrystal and not just NCrystalDev.

    #TODO: keep structure from original NCrystalDev repo : NCCore, NCUtils, etc.

    #TODO: Wrap <reporoot>/tests as well.

    #TODO: Absorb NCVersion into ncapi.h and NCDefs.hh (that way we can also
    #avoid versions hardcoded in both NCVersion.hh and ncrystal.h, and the
    #CMakeLists.txt files, and have just the version in VERSION +
    #ncrystal_python/NCrystal/__init__.py. Idea: we could even avoid the latter
    #if we make the sdists for ncrystal_python dynamically, (but also generate
    #the pyproject.toml for them).

    #TODO: Also autogenerate NCrystal.hh, adding all public includes? Problem
    #is, we need both python and cmake code for it. Might be easier with a
    #static one + a standalone test that verifies we did not forget anything.

    from . import dirs

    #NCrystal headers:
    ncapi_loc_in_genroot = 'extraincpath/NCrystal/ncapi.h'
    fgen_ncapi = add_file( ncapi_loc_in_genroot,
                           content=ncapi_contents() )
    #FIXME: One custom extdep per pkg, and in an extraincpath, so we enforce the
    #deps to be correct.
    create_custom_extdep( 'NCDevHeaders',
                          cflags = (f'-I{dirs.genroot}/extraincpath '
                                    f'-I{dirs.srcroot}/include '
                                    '-fno-math-errno') )

    #NCrystal source file compilation:
    from .ncrystalsrc import load_components
    name2comp = load_components()
    all_ncsbpkgs = []
    for name,comp in sorted(name2comp.items()):
        #TODO: comp.srcdir_hdrs, comp.srcfiles
        sbpkgname = cfg.sbpkgname_ncrystal_comp(name)
        all_ncsbpkgs.append(sbpkgname)
        extdeps = ['NCDevHeaders']
        if name=='factories':
            extdeps.append('DL')
        create_pkginfo( sbpkgname,
                        extdeps = extdeps,
                        pkg_deps = [ cfg.sbpkgname_ncrystal_comp(n)
                                     for n in comp.direct_depnames ],
                        extra_cflags = [f'-I{dirs.srcroot}/src',
                                        f'-DNCRYSTAL_DATADIR={dirs.datadir}'] )#Fixme NCRYSTAL_DATADIR also only for factories
        for sf in (comp.srcdir_hdrs+comp.srcfiles):
            add_file( f'pkgs/{sbpkgname}/libsrc/{sf.name}', link_target = sf )

        #For traditional simplebuild includes:
        for sf in comp.hdrfiles or []:
            add_file( f'pkgs/{sbpkgname}/libinc/{sf.name}', link_target = sf )

    #Also symlink ncapi.h and NCrystal.hh:
    add_file( 'pkgs/%s/libinc/ncapi.h'%cfg.sbpkgname_ncrystal_comp('core'),
              link_target = dirs.genroot.joinpath(ncapi_loc_in_genroot) )
    add_file( 'pkgs/%s/libinc/NCrystal.hh'%cfg.sbpkgname_ncrystal_ncrystalhh,
              link_target = dirs.srcroot.joinpath('include/NCrystal/NCrystal.hh') )#fixme: autogen?

    #create_pkginfo( cfg.sbpkgname_ncrystal_lib,
    #                extdeps = ['DL','NCDevHeaders'],
    #                extra_cflags = [f'-I{dirs.srcroot}/src',
    #                                f'-DNCRYSTAL_DATADIR={dirs.datadir}'] )
    #for sf in (dirs.srcroot/'src').glob('*c'):
    #    add_file( f'pkgs/{cfg.sbpkgname_ncrystal_lib}/libsrc/{sf.name}',
    #              link_target = sf )

    #NCrystalDev package (python module):
    assert cfg.sbpkgname_ncrystal_lib in all_ncsbpkgs
    create_pkginfo( 'NCrystalDev', pkg_deps=all_ncsbpkgs )
    for sf in (dirs.pysrcroot/'NCrystal').glob('*.py'):
        add_file( f'pkgs/NCrystalDev/python/{sf.name}', link_target = sf )
    #Special marker used by _locatelib.py:
    add_file( 'pkgs/NCrystalDev/python/_is_sblddevel.py', content='' )

    #We added '#define NCRYSTAL_NAMESPACED_ENVVARS' in ncapi.h, so we also let
    #the python layer know:
    add_file( 'pkgs/NCrystalDev/python/_do_namespace_envvars.py', content='' )

    #Commandline scripts:
    create_pkginfo( cfg.sbpkgname_ncrystal_cli, pkg_deps=['NCrystalDev'])
    for sf in (dirs.pysrcroot/'NCrystal').glob('_cli_*.py'):
        cliname = sf.name[len('_cli_'):-len('.py')]
        sbscriptname = 'tool' if cliname == 'nctool' else cliname
        content = f"""#!/usr/bin/env python3
import NCrystalDev.{sf.stem} as mod
mod.main()
"""
        add_file( f'pkgs/{cfg.sbpkgname_ncrystal_cli}/scripts/{sbscriptname}',
                  content=content,
                  make_executable = True)
    for subpath,fn in [('include/NCrystal/internal/utils/NCCFileUtils.hh',
                        'NCCFileUtils.h'),
                       ('src/utils/NCCFileUtils.cc','NCCFileUtils.c'),
                       ('app_config/main.c','main.c')]:
        sf = dirs.srcroot / subpath
        add_file( f'pkgs/{cfg.sbpkgname_ncrystal_cli}/app_config/{fn}',
                  link_target = sf )
    add_file( f'pkgs/{cfg.sbpkgname_ncrystal_cli}/app_config/ncconfig_autogen.h',
              content = ncconfig_h_contents() )

    #Data files (even though NCLib already knows the path via a define):
    create_pkginfo( cfg.sbpkgname_ncrystal_data )
    for sf in dirs.datadir.glob('*.ncmat'):
        add_file( f'pkgs/{cfg.sbpkgname_ncrystal_data}/data/{sf.name}',
                  link_target = sf )

    #Examples:
    example_src =[]
    example_src_g4 =[]
    for f in dirs.exsrcroot.glob('*.c*'):
        if 'g4sim' in f.name:
            example_src_g4.append(f)
        else:
            example_src.append(f)

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

    create_pkginfo( cfg.sbpkgname_ncrystal_examples,
                    pkg_deps=[cfg.sbpkgname_ncrystal_lib] )
    add_compiled_examples(example_src,cfg.sbpkgname_ncrystal_examples)

    #Geant4 (until it moves elsewhere):
    create_custom_extdep( 'NCG4DevHeaders',
                          cflags = f'-I{dirs.ncg4srcroot}/include ' )
    create_pkginfo( cfg.sbpkgname_ncrystal_geant4,
                    pkg_deps=[cfg.sbpkgname_ncrystal_lib],
                    extdeps = ['Geant4','NCG4DevHeaders'],
                    extra_cflags = [f'-I{dirs.ncg4srcroot}/src'],
                   )
    for sf in (dirs.ncg4srcroot/'src').glob('*c'):
        add_file( f'pkgs/{cfg.sbpkgname_ncrystal_geant4}/libsrc/{sf.name}',
                  link_target = sf )

    for sf in dirs.ncg4srcroot.joinpath('include/G4NCrystal').glob('*.hh'):
        add_file( f'pkgs/{cfg.sbpkgname_ncrystal_geant4}/libinc/{sf.name}',
                  link_target = sf )

    add_compiled_examples( example_src_g4,
                           cfg.sbpkgname_ncrystal_geant4,
                           override_name = 'example' )

def main():
    define_files()
    create_files()
