#!/usr/bin/env python3

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

# NEEDS: numpy

import NCrystalDev as NC
from NCTestUtils.common import fix_ncrystal_version_printouts
import NCrystalDev.cli as nc_cli
import pathlib
import contextlib
import shlex
from NCTestUtils.common import ( print_text_file_with_snipping,
                                 ensure_error,
                                 work_in_tmpdir,
                                 fmt_args_as_str )

#default is to encode numbers in produced .laz/.lau files with 14 digits of
#precision (%.14g), but for robustness of the unit test we reduce this in any
#printed content to avoid spurious false positives due to FP instabilities:
test_precision = 10

#for re-loading those files and comparing the precision with the original
#NCrystal physics. Note, this could be a lot more precise if running with the
#default FP precision of ncmat2hkl (test_precision=14), but we can't since then
#we would get spurious test failures in the reflog diffing.:
max_diff_lvl = 1e-7

def test_pyapi( cfgstr, fmt, nstart = 30, nend = 20 ):
    from NCrystalDev.mcstasutils import cfgstr_2_hkl
    kwargs = dict(cfgstr=cfgstr,
                  tgtformat=fmt,
                  fp_format=f'%.{test_precision}g')
    args_str = fmt_args_as_str( **kwargs )
    hr=f"============= PyAPI >>{args_str}<< ===================="
    print(hr)
    res = '\n'.join(cfgstr_2_hkl(**kwargs))+'\n'
    print_text_file_with_snipping( res,
                                   nstart=nstart,
                                   nend=nend,
                                   prefix='res>')
    print('='*len(hr))
    return res

def test_cli( args, *,
              nstart = 30, nend = 20,
              outfile = None,
              in_tmp_dir = True ):
    if isinstance(args,str):
        args = shlex.split(args)
    args = args[:] + [f'--override-prec={test_precision}']
    is_stdout = outfile is None or outfile=='stdout'
    if is_stdout:
        import io
        iobuf = io.StringIO()
        cm_output = contextlib.redirect_stdout(iobuf)
    else:
        iobuf = None
        cm_output = contextlib.nullcontext()

    hr=f"============= CLI >>{shlex.join(args)}<< ===================="
    print(hr)
    created_content = None
    ctx = work_in_tmpdir if in_tmp_dir else contextlib.nullcontext
    with ctx(), cm_output:
        nc_cli.run('ncmat2hkl',*args)
        if outfile not in ('stdout',None):
            created_content = pathlib.Path(outfile).read_text()
            print_prefix = 'OUTFILE>'

    if iobuf:
        created_content = iobuf.getvalue()
        print_prefix = 'STDOUT>'
    print_text_file_with_snipping( created_content,
                                   nstart=nstart,
                                   nend=nend,
                                   prefix=print_prefix)
    print('='*len(hr))
    return created_content

def main():
    fix_ncrystal_version_printouts()
    import sys
    do_plot = '--plot' in sys.argv[1:]
    if do_plot:
        #from PyAna import *
        pass
    else:
        import NCTestUtils.enable_fpe # noqa F401

    test_cli(['-h'])

    def testcfg(cfgstr,*,outfile,use_pyapi):
        test_cfgstr(cfgstr,do_plot=do_plot,
                    outfile=outfile,use_pyapi=use_pyapi)

    outfile_opts = [None,'stdout','dummy.FMT']
    cfgstrs = [
        'stdlib::Al_sg225.ncmat;dcutoff=0.6',
        'stdlib::CaSiO3_sg2_Wollastonite.ncmat;dcutoff=3.2',
        'stdlib::Y2SiO5_sg15_YSO.ncmat;dcutoff=2.0',
        'stdlib::Mg2SiO4_sg62_MagnesiumSilicate.ncmat;dcutoff=1.2',
        'stdlib::PbO-alpha_sg129_Litharge.ncmat;dcutoff=0.9',
        'stdlib::MgCO3_sg167_MagnesiumCarbonate.ncmat;dcutoff=1.0',
        'stdlib::Al2O3_sg167_Corundum.ncmat;dcutoff=1.5',
        'stdlib::ZnO_sg186_ZincOxide.ncmat;dcutoff=0.5',
        'stdlib::C_sg194_pyrolytic_graphite.ncmat'#has embedded cfg
    ]


    for i,cfgstr in enumerate(cfgstrs):
        testcfg( cfgstr,
                 use_pyapi = bool(i%2),
                 outfile=outfile_opts[i%len(outfile_opts)] )

    with ensure_error(NC.NCFileNotFound,
                      ('Requested factory "stdlib" can'
                       ' not provide data: "unknown_file.ncmat".')):
        test_cli(['--format','laz','stdlib::unknown_file.ncmat'])

    with ensure_error(RuntimeError,
                      ('this converter does not'
                       ' handle non-crystalline materials')):
        test_cli(['--format','laz','solid::B4C/2.52gcm3/B_is_0.95_B10_0.05_B11'])

    data = """NCMAT v1
@CELL
    lengths 4.04958 4.04958 4.04958
    angles 90. 90. 90.
@ATOMPOSITIONS
    Al 0. 0.5 0.5
    Al 0. 0. 0.
    Al 0.5 0.5 0.
    Al 0.5 0. 0.5
@DEBYETEMPERATURE
    Al   410.3542
"""
    NC.registerInMemoryFileData('Al_nosg.ncmat',data)
    NC.registerInMemoryFileData('Al.ncmat',data+'@SPACEGROUP\n    225\n')

    testcfg( 'virtual::Al.ncmat;dcutoff=0.8',use_pyapi = True,outfile=None)
    with ensure_error(RuntimeError,
                      ('this converter does not handle crystalline'
                       ' materials without spacegroup number')):
        test_cli(['--format','laz','virtual::Al_nosg.ncmat;dcutoff=0.8'])

    with ensure_error(RuntimeError,
                      ('this converter does not handle configurations'
                       ' with density overrides (since resulting files might'
                       ' not load correctly by all clients)')):
        test_cli(['--format','laz','stdlib::Al_sg225.ncmat;dcutoff=0.6;density=0.8x'])


    from argparse import ArgumentError
    with ensure_error(ArgumentError,
                      ('Please specify output type with --format, either'
                       ' --format=laz or --format=lau (run'
                       ' with --help for details)')):
        test_cli(['stdlib::Al_sg225.ncmat'])
    with ensure_error(ArgumentError,
                      ("argument --format/-f: invalid choice: 'LAU'"
                       " (choose from laz, lau)")):
        test_cli(['stdlib::Al_sg225.ncmat','--format','LAU'])

    with work_in_tmpdir():
        cmdargs=['stdlib::Al_sg225.ncmat;dcutoff=0.6',
             '--format','lau',
             '-o','dummy.lau']
        test_cli(cmdargs, in_tmp_dir=False,outfile='dummy.lau')
        with ensure_error(ArgumentError,
                          'Output destination already exists: dummy.lau'):
            test_cli(cmdargs, in_tmp_dir=False,outfile='dummy.lau')

def cfgstr2hkl(cfgstr,fmt,outfile,use_pyapi):
    assert fmt in ('lau','laz')
    common = dict(nstart=100, nend=4)
    if use_pyapi:
        return test_pyapi(cfgstr,fmt,**common)
    if outfile is not None and 'FMT' in outfile:
        outfile = outfile.replace('FMT',fmt)
    args=[f'--format={fmt}',cfgstr]
    if outfile is not None and outfile!='stdout':
        args+=['-o',outfile]
    res = test_cli( args,outfile = outfile,**common )
    return res

def test_cfgstr(cfgstr, do_plot,outfile,use_pyapi):

    import numpy as np

    #NC.registerInMemoryFileData('foo.laz'c,fgstr2hkl(cfgstr,'laz'))
    def bragg_only(cfgstr):
        return cfgstr + ';inelas=0;elas=0;bragg=1'
    def make_oriented(cfgstr):
        return cfgstr+(';dir1=@crys_hkl:0,0,1@lab:0,0,1'
                       ';dir2=@crys_hkl:0,1,0@lab:0,1,0'
                       ';mos=10deg;dirtol=180deg')
    sc_orig = NC.createScatter(bragg_only(cfgstr))
    sc_orig_oriented = NC.createScatter(make_oriented(bragg_only(cfgstr)))
    hr="="*100+'\n'
    def prtitle( title ):
        print('\n'+hr*3 + f"====> {title} <=====\n" +hr*3+'\n')
    for fmt in ('laz','lau'):
        fn=f'foo.{fmt}'


        def smart_dump(infoobj):
            #like infoobj.dump(verbose=2), but snipping lines:
            dump_full = infoobj.dump_str(verbose=2)
            df_lines = dump_full.splitlines()
            for nstart,line in enumerate(df_lines):
                if 'HKL' in line and line.lstrip().lower().startswith('hkl planes'):
                    break
            print_text_file_with_snipping( dump_full,
                                           nstart=nstart+15,
                                           nend=5,
                                           prefix='')

        info=NC.createInfo(cfgstr)
        if do_plot:
            prtitle(f'NCrystal dump of cfgstr={cfgstr}')
            info.dump()
        wls = np.linspace(0.001,info.braggthreshold*1.2,20000)
        wls_sparse = np.linspace(wls[0],wls[-1],2000)

        prtitle(f'Conversion of {cfgstr} to {fmt}')
        _created_content = cfgstr2hkl(cfgstr,fmt,outfile,use_pyapi=use_pyapi)
        NC.registerInMemoryFileData(fn,_created_content)
        rd=NC.createTextData(fn).rawData
        assert rd == _created_content
        assert rd.endswith('\n')

        prtitle(f'NCrystal dump of converted {fmt}')
        smart_dump(NC.createInfo(fn))

        prtitle(f'Testing bragg component of converted {fmt} vs. original')
        sc=NC.createScatter(bragg_only(fn))
        sc.dump()
        xs=sc.xsect(wl=wls)
        xs_orig=sc_orig.xsect(wl=wls)
        diff=abs(xs-xs_orig).max()
        if do_plot:
            print (f'XS max-abs-diff converted {fmt} vs. original:',diff)
        else:
            print (f'XS max-abs-diff converted {fmt} vs. original < {max_diff_lvl}?: %s'%('yes' if diff<max_diff_lvl else 'no'))
        if do_plot:
            import matplotlib.pyplot as plt
            plt.plot(wls,xs,label='via laz')
            plt.plot(wls,xs_orig,label='orig',ls='--',lw=3,alpha=0.5)
            plt.legend()
            plt.show()
        assert diff < max_diff_lvl

        #Test that can be used for single crystal, with huge mosaicity so we can
        #get non-trivial cross sections without careful aiming:
        sc=NC.createScatter(make_oriented(bragg_only(fn)))
        sc.dump()
        xs = sc.xsect(wl=wls_sparse, direction=(0,0,1))
        xs_orig = sc_orig_oriented.xsect(wl=wls_sparse, direction=(0,0,1))
        diff=abs(xs-xs_orig).max()
        if do_plot:
            print (f'SC XS max-abs-diff converted {fmt} vs. original:',diff)
        else:
            print (f'SC XS max-abs-diff converted {fmt} vs. original < {max_diff_lvl}?: %s'%('yes' if diff<max_diff_lvl else 'no'))
        if do_plot:
            import matplotlib.pyplot as plt
            plt.plot(wls_sparse,xs,label='via laz')
            plt.plot(wls_sparse,xs_orig,label='orig',ls='--',lw=3,alpha=0.5)
            plt.legend()
            plt.show()

        assert diff < max_diff_lvl

if __name__ == '__main__':
    main()
