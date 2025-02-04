
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

from ._cliimpl import ( create_ArgumentParser,
                        cli_entry_point,
                        print )
from ._common import write_text
from .constants import ( constant_planck, constant_boltzmann )
from .vdos import vdos_units_2_eV
import pathlib
units_2_fact = vdos_units_2_eV
units_opts = ', '.join(sorted(units_2_fact.keys()))

_cache_np=[None]
def get_numpy():
    if _cache_np[0] is None:
        from ._numpy import _ensure_numpy, _np
        _ensure_numpy()
        _cache_np[0] = _np
    return _cache_np[0]

def parseArgs( progname, arglist, return_parser=False ):
    import textwrap
    helpw = 60
    descrw = helpw + 22

    descr="""WARNING: This script is somewhat outdated. Many users will instead
wish to use the PhononDOSAnalyser helper class from the Python API in the
NCrystal.vdos python module. For more information, see the jupyter notebook
about "Adding phonon information" on:
https://github.com/mctools/ncrystal-notebooks

Script which can read vdos curves from either .ncmat files or simple two-column
ascii text files (the two columns being energy grid and density) and help users
prepare output suitable for inclusion in .ncmat files. When the input is not an
.ncmat file, the user must specify the energy grid units by adding a line in the file
like for instance "#unit:THz" or "#unit:meV".

In case of NCMAT files with more than one VDOS, just post-fix the filename with
the element to investigate, separated by '@@', e.g. "Al2O3.ncmat@@O". In case of
.txt files with more than 2 columns, select the column representing the VDOS
density in the same manner, e.g. "Al2O3.txt@2" (the first column is always the
energy or frequency).

Thus it is possible to plot the curve, compare against vdos curves from other
.ncmat files, apply low-E truncation (cutoff), regularise binning, and overlay
with an ideal Debye spectrum for a given Debye energy or Debye
temperature. Finally, when running without --plot, it will output the resulting
spectrum into format which is ready for inclusion in .ncmat files.

"""
    descr = '\n\n'.join(textwrap.fill(e,descrw) for e in descr.split('\n\n'))

    from argparse import RawTextHelpFormatter
    parser = create_ArgumentParser( prog = progname,
                                    description=descr,
                                    formatter_class=RawTextHelpFormatter )
    def wrap(t):
        return textwrap.fill(' '.join(t.split()),width=helpw)

    parser.add_argument("FILE",help=wrap("Either .ncmat file with VDOS or a"
                                         " text file with two colums of data:"
                                         " egrid and density."))
    parser.add_argument('--plot','-p', action='store_true',help='Plot extracted spectrum')
    parser.add_argument("--cutoff",'-c',nargs='+',default=None,type=float,
                        help=wrap("""Emin cutoff points in eV (more than one can
                        be provided for simultaneous inspection with --plot)."""))
    parser.add_argument("--unit",'-u',default='meV',type=str,
                        help=wrap('Visualisation unit (ignored unless'
                                  ' --plot is supplied). Defaults to meV. '
                                  f'Possible options are {units_opts}.'))
    parser.add_argument("--ref",'-r',nargs='+',
                        action='append',type=str,
                        help=wrap("""Optionally provide list of .ncmat files
                        with vdos data, to superimpose on plots."""))
    parser.add_argument("--forceregular",'-f',type=int,nargs='?',default=0,
                        help=wrap("""Optionally provide this argument to
                        reparameterise with that amount of linearly spaced
                        points in [emin,emax+epsilon], where epsilon is chosen
                        so the grid can be extended to 0 with a whole number of
                        bins. This format will be directly used by NCrystal
                        without on-the-fly reparameterisation upon loading."""))
    parser.add_argument("--debye",'-d',nargs='?',default='',type=str,
                        help=wrap("""Set to debye temperature (unit K) or egrid
                        point (units like meV, eV, THz, ...) in order to plot
                        Debye spectrum with that parameter on top."""))
    parser.add_argument('--g1',default=0.0,type=float,metavar='TEMP',
                        help=wrap("""Use with --plot to show Sjolander's G1
                        function at the indicated temperature value (in kelvin)
                        instead of the DOS directly. This is the Symmetric G1
                        fct without a detailed balance factor, and it will be
                        plotted assuming gamma0=1.0."""))
    parser.add_argument('--stdout',action='store_true',help=wrap("""Produce no
                        output file but print vdos_egrid and vdos_density lines
                        to stdout (for scripting)"""))

    #We could expand the env var name to account for any namespace:
    #  from ._common import expand_envname
    #  dpienvvar = expand_envname('DPI')
    #But that will break the unit tests, so we just use:
    dpienvvar = 'NCRYSTAL_DPI'

    dpi_default=200
    parser.add_argument('--dpi', default=-1,type=int,
                        help=wrap(f"""Change plot resolution. Set to 0 to leave
                        matplotlib defaults alone.  (default value is
                        {dpi_default}, or whatever the {dpienvvar} env var is
                        set to)."""))

    if return_parser:
        return parser
    args = parser.parse_args(arglist)

    if args.forceregular is None:
        parser.error('Missing argument (number of points) to --forceregular.')

    if args.dpi==-1:
        from ._common import ncgetenv_int_nonneg
        args.dpi = ncgetenv_int_nonneg('DPI',dpi_default)
    if args.dpi > 3000:
        parser.error('Too high DPI value requested.')

    if args.ref:
        args.ref = [item for sublist in args.ref for item in sublist]
    if args.ref and not args.plot:
        parser.error('Option --ref requires --plot')
    if args.unit and args.unit not in units_2_fact:
        parser.error(f'Unknown unit {args.unit}. Valid options are {units_opts}')
    if args.debye and not args.plot:
        parser.error('Option --debye requires --plot')
    if args.stdout and args.plot:
        parser.error('Option --stdout can not be used with --plot')
    if args.cutoff and len(args.cutoff)>1 and not args.plot:
        parser.error('Option --cutoff can only have one argument when not using --plot')
    if args.cutoff and len(args.cutoff)>1 and args.forceregular:
        parser.error('Option --cutoff can only have one argument when using --forceregular')

    if args.debye:
        if args.debye.endswith('K'):
            args.debye = float(args.debye[0:-1])*constant_boltzmann
        else:
            #find (longest, so "meV" does not trigger "eV") fitting unit:
            ll=[ (len(u),u) for u in units_2_fact.keys() if args.debye.endswith(u) ]
            ll.sort()
            if not ll:
                parser.error("Option --debye requires unit (see --help)")
            unit = ll[-1][1]
            args.debye = units_2_fact[unit] * float(args.debye[0:-len(unit)])


    return args

def create_argparser_for_sphinx( progname ):
    return parseArgs(progname,[],return_parser=True)

def decodeFileName(filename):
    if '@@' in filename:
        path,select = filename.split('@@',1)
    else:
        path,select = filename,None
    import os
    bn=os.path.basename(path)
    return dict(path=path,select=select,basename=bn,
                title=bn if not select else '%s in %s'%((select if not select.isdigit(
                ) else f'column #{select}'),bn))

def getVDOSFromFile(fn):
    fnd = decodeFileName(fn)
    if fnd['basename'].endswith('.ncmat'):
        return getVDOSFromNCMAT(fn)
    return getVDOSFromTXT(fn)

def getVDOSFromTXT(fn):
    fn = decodeFileName(fn)
    select=None
    if fn['select']:
        assert fn['select'].isdigit(),"selection must be column number"
        select = int(fn['select'])
    #figure out unit:
    unit=None
    with open(fn['path']) as fh:
        for ll in fh:
            if ll.startswith('#') and 'unit' in ll:
                _=ll.split('#',1)[1].split(':',1)
                if not len(_)==2:
                    continue
                unit=_[1].strip()
                if unit not in units_2_fact.keys():
                    raise RuntimeError(f'Unknown unit "{unit}" specified in {fn["path"]}. Valid choices are: {units_opts}')
                break
    if not unit:
        raise RuntimeError(f'Missing energy/frequency unit in {fn["path"]}. Please put initial line with content like "#unit:THz". Valid units are: {units_opts}')
    usecols=None
    if select is not None:
        assert select>0
        usecols=(0,select)
    _ = get_numpy().genfromtxt(fn['path'],dtype=[('egrid','f8'),('density','f8')],usecols=usecols)
    egrid=_['egrid'].copy()
    density=_['density'].copy()
    density /= density.max()
    return egrid.copy()*units_2_fact[unit],density.copy()

def getVDOSFromNCMAT(fn):
    fnd = decodeFileName(fn)
    from . import core as NC
    info = NC.createInfo(fnd['path'])
    select = fnd['select']
    di_vdos = [di for di in info.dyninfos if isinstance(di,NC.Info.DI_VDOS)]
    if select is not None:
        ll=[]
        for di in di_vdos:
            dl=di.atomData.displayLabel()
            if dl!=select:
                print(f"NB: Ignoring (due to selection) VDOS for element {dl} in file {fn}.")
            else:
                ll+=[di]
        if not ll:
            raise RuntimeError(f'ERROR: Could not find VDOS in file {fn} for selected element: {select}')
        if not len(ll)==1:
            raise RuntimeError(f'ERROR: Multiple VDOS entries in file {fn} for selected element: {select} (which is rather odd!)')
        di_vdos = ll
    if len(di_vdos)>1:
        s=' '.join(di.atomData.displayLabel() for di in di_vdos)
        raise RuntimeError(f"ERROR: Multiple VDOS entries found in file {fn}. Please select one of them (by putting \"@@<element>\" after the file-name): {s}")
    elif not di_vdos:
        raise RuntimeError(f"ERROR: No vdos found in file {fn}")
    eg,ds = di_vdos[0].vdosOrigEgrid(),di_vdos[0].vdosOrigDensity()
    ds /= ds.max()
    if len(eg)==2:
        get_numpy()
        from ._numpy import _np_linspace
        eg = _np_linspace(eg[0],eg[1],len(ds))
    return eg.copy(),ds.copy()




@cli_entry_point
def main( progname, arglist ):
    args = parseArgs( progname, arglist )
    from ._numpy import _ensure_numpy, _np
    _ensure_numpy()

    file_decoded = decodeFileName(args.FILE)
    args_file_basename = file_decoded['basename']
    egrid,density = getVDOSFromFile(args.FILE)

    np = _np

    assert len(egrid) == len(density)
    print (f"Loaded VDOS with {len(density)} grid points from {args_file_basename}")

    np = get_numpy()
    def numpy_is_sorted(a):
        return np.all(a[:-1] <= a[1:])
    def numpy_is_strongly_sorted(a):
        return np.all(a[:-1] < a[1:])

    if not numpy_is_strongly_sorted(egrid):
        for i in range(len(egrid)-1):
            if not egrid[i] < egrid[i+1]:
                print("Problems detected in egrid points with values ",egrid[i],"and",egrid[i+1])
        raise RuntimeError('ERROR: egrid values (first column) of input file are not in sorted'
                         +' (ascending) order, or there are identical elements.')

    cutoffs=[]
    if args.cutoff:
        for c in args.cutoff:
            if c >= egrid[-1]:
                raise RuntimeError(f'ERROR: Cutoff value {c} is higher than highest point in egrid')
            i=np.searchsorted(egrid,c)
            assert i==0 or egrid[i-1]<c
            assert egrid[i]>=c
            cutoffs+=[ (i, egrid[i] ) ]
            print(f" => Mapping cutoff value of {c} to grid point at {cutoffs[-1][1]}")


    if args.forceregular or (not args.plot):
        if applyCutoff(egrid,density,cutoffs)[0][0]<=1e-5:
            raise RuntimeError(f"""
            ERROR: The first value in the loaded egrid is {egrid[0]} which is less than 1e-5eV.
            This is not allowed when using --forceregular or when not using --plot.
            Please use the --cutoff parameter to remove lowest part of input spectrum (perhaps
            after investigating the cutoff value with --plot).
            """)


    if args.forceregular:
        regularised_egrid,regularised_density = regularise(*applyCutoff(egrid,density,cutoffs),args.forceregular)

    if args.plot:
        vis_unit=args.unit
        vis_unit_fact = 1.0 / units_2_fact[vis_unit]

        def plt_plot(egrid,dos,*aargs,**kwargs):
            if args.g1 and args.g1 > 0.0:
                #G1 = f(E)/(E*2*gamma0*sinh(E/2kT)) [symmetric G1 that is]
                #u = E/2kT
                #asymmetric means another factor of exp(-u).
                #And exp(-u)/2sinh(u) = exp(-u) / (exp(u)-exp(-u) )  = 1 / ( exp(2u)-1)
                #And exp(+u)/2sinh(u) = exp(+u) / (exp(u)-exp(-u) )  = 1 / ( 1-exp(-2u) )
                #
                #So with gamma0=0 we get:
                egrid_eV = egrid.copy() / vis_unit_fact
                T = args.g1
                #gamma0 = 1.0
                two_u = egrid_eV / ( constant_boltzmann * T )
                #g1asym_neg = dos / ( egrid_eV * -np.expm1(-two_u) )
                #g1asym_pos = dos / ( egrid_eV * np.expm1(two_u) )
                g1sym = dos / (egrid_eV*np.sinh(0.5*two_u))
                plt.plot(egrid,g1sym,*aargs,**kwargs)
            else:
                plt.plot(egrid,dos,*aargs,**kwargs)

        import matplotlib as mpl
        mpl.rcParams['figure.dpi']=args.dpi
        #ability to quit plot windows with Q:
        if 'keymap.quit' in mpl.rcParams and 'q' not in mpl.rcParams['keymap.quit']:
            mpl.rcParams['keymap.quit'] = tuple(list(mpl.rcParams['keymap.quit'])+['q','Q'])
        import matplotlib.pyplot as plt
        plt.xlabel(vis_unit)
        plt_plot(egrid*vis_unit_fact,density,'o-',label=file_decoded['title'])
        if args.forceregular:
            plt_plot(regularised_egrid*vis_unit_fact,regularised_density,'x-',label='regularised')
        from ._numpy import _np_linspace
        for c_idx, c_val in cutoffs:
            d=density[c_idx]
            # f(x)=k*x^2, f(c_val)=d<=> k*c_val^2 = d <=> k = d/c_val^2
            x=_np_linspace(0.0,c_val,3000)
            plt_plot(x*vis_unit_fact,(d/c_val**2)*(x**2),label=f'with cutoff {c_val} eV')
        if args.debye:
            x=_np_linspace(0.0,max(egrid.max(),args.debye),1000)
            y = np.where(  x<=args.debye, x**2 * ( density.max() / args.debye**2 ), 0.0 )
            plt_plot(x*vis_unit_fact,y,
                     label=f'Debye spectrum (E_Debye={1000*args.debye:.5}meV, T_Debye={args.debye/constant_boltzmann:.5}K)')
        for r in (args.ref or []):
            eg,ds = getVDOSFromFile(r)
            plt_plot(eg*vis_unit_fact,ds,label=decodeFileName(r)['title'])
        plt.legend()
        plt.title(file_decoded['title'])
        plt.grid(ls=':')
        plt.show()
        return

    if args.forceregular:
        egrid, density = regularised_egrid,regularised_density
    else:
        egrid, density = applyCutoff(egrid,density,cutoffs)

    #Check if egrid is linspace:
    binwidth = (egrid[-1]-egrid[0])/(len(egrid)-1)
    is_linspace=True
    if not args.forceregular:
        for i in range(len(egrid)-1):
            bw=egrid[i+1]-egrid[i]
            if abs(binwidth-bw)>0.01*binwidth:
                is_linspace=False
                break
        if is_linspace:
            print('NB: Detected linearly spaced input egrid')

    #normalise so unity at highest point (gives smaller file sizes):
    density /= density.max()

    #remove excess trailing zeros
    while len(density)>10 and density[-2]==0.0 and density[-1]==0.0:
        density = density[0:-1]
        egrid = egrid[0:-1]

    from .ncmat import formatVectorForNCMAT

    egrid_cnt =''
    if is_linspace:
        egrid_cnt += f'  vdos_egrid {egrid[0]:.14} {egrid[-1]:.14}'
    else:
        egrid_cnt += formatVectorForNCMAT('vdos_egrid',egrid)
    egrid_cnt += '\n'
    egrid_cnt += formatVectorForNCMAT('vdos_density',density)

    if args.stdout:
        print("<<<GENERATED-CONTENT-BEGIN>>>")
        print(egrid_cnt,end='')
    else:
        outfn=pathlib.Path('converted_output.ncmat')
        content = f"""NCMAT v5
#Autogenerated file from {args_file_basename}.
@DENSITY
  1.0 g_per_cm3 #FIX{'ME'}!! Please replace with proper value, or remove and optionally provide crystal structure!
@DYNINFO
  element  <UNKNOWN-PLEASE-EDIT>
  fraction 1
  type     vdos\n"""
        content += egrid_cnt
        content += '\n'
        write_text( outfn, content )
        print(f"Wrote {outfn}")

def applyCutoff(egrid,density,cutoffs):
    if cutoffs:
        assert len(cutoffs)==1
        c_idx,c_val = cutoffs[0]
        return egrid[c_idx:], density[c_idx:]
    return egrid,density

def regularise(egrid,density,n):

    #first step back from any zeroes at the upper end:
    i=1
    while i < len(density) and density[-i]==0.0:
        i += 1
    safepeel = i-2
    if safepeel>=1:
        print (f"Ignoring {safepeel} last points while regularising since last {safepeel+1} points are 0.")
        egrid,density = egrid[0:-(safepeel)],density[0:-(safepeel)]
    emin,emax=egrid[0],egrid[-1]
    print('old range',emin,emax)
    THZ = constant_planck*1e12
    print('old range [THZ]',emin/THZ,emax/THZ)

    np = get_numpy()
    for k in range(1,1000000000):
        #k is number of bins below emin, an integral number by definition in a regularised grid.
        binwidth = emin/k
        nbins=int(np.floor((emax-emin)/binwidth))+1
        eps = (emin+nbins*binwidth)-emax
        assert eps>=0.0
        if nbins+1 >= n:
            break
    n=nbins+1
    binwidth = emin/k
    new_emax = emin + (n-1) * binwidth
    if abs( (new_emax-binwidth) - emax ) < 1e-3*binwidth:
        nbins -= 1
        n -= 1
        new_emax -= binwidth
    print (f" ==> Choosing regular grid with n={n} pts from emin={emin} to emax={new_emax} ({new_emax-emax} beyond old emax)")
    assert new_emax >= emax-binwidth*1.001e-3
    from ._numpy import _np_linspace
    new_egrid = _np_linspace(emin,new_emax,n)
    test=new_egrid[0] / ( (new_egrid[-1]-new_egrid[0])/(len(new_egrid)-1) )
    assert abs(round(test)-test)<1e-6,f'{test}'
    np = get_numpy()
    new_density = np.interp(new_egrid,egrid,density, left=0.0, right=0.0)
    print('last density values in new grid:',new_density[-5:])
    return new_egrid,new_density
