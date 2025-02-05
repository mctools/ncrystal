
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
                        print,
                        warn )

from . import constants as nc_constants
from . import core as nccore

import warnings
import pathlib

__pynecache=[None]
def import_pyne():
    global __pynecache
    if __pynecache[0] is not None:
        return __pynecache[0]

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")#silence annoying pyne.endf QAWarning's.
            import pyne.endf
    except ImportError:
        raise RuntimeError('Could not "import pyne.endf". You probably need to install PyNE (see http://pyne.io/)')

    #Monkey-patch pyne.endf reader, which was ignoring optional sections with teff
    #curves of non-principal atoms (see
    #https://github.com/pyne/pyne/issues/1209). The fix was accepted upstream, so
    #eventually we can remove this monkey-patching once we think all users will have
    #recent enough PyNE installations:
    _orig_rti = pyne.endf.Evaluation._read_thermal_inelastic
    def _fix_rti(_self):
        _orig_rti(_self)
        inel=_self.thermal_inelastic
        B=inel['B']
        for ii in range(len(B)//6-1):
            if B[6*(ii+1)]==0.0:
                if 'teff_%i'%(ii+1) in inel:
                    #If upstream already includes fix, do nothing.
                    continue
                _, teff = _self._get_tab1_record()
                inel['teff_%i'%(ii+1)] = teff#Store as teff_1, teff_2, etc. for now.
    pyne.endf.Evaluation._read_thermal_inelastic=_fix_rti
    __pynecache[0] = pyne
    return pyne

def elementZToName(z):
    d = nccore.atomDB(z,throwOnErrors=False)
    return d.elementName() if d else None

def guessElementFromMass(mass_amu):
    elem_data=[ (abs(d.averageMassAMU()-mass_amu),d) for d in [nccore.atomDB(z,throwOnErrors=False) for z in range(1,120)] if d ]
    elem_data.sort()
    print (f'Guessing mass={mass_amu:.8}u is element {elem_data[0][1].elementName()}')
    return elem_data[0][1].elementName()

def guessElementName(ZA,mass_amu):
    Z=ZA//1000
    if ZA==155:
        return 'Y'#YH2, whatever it is
    if ZA==140:
        return 'H'#benzene...  C6H6. Guessing H is dominant. But tsl-benzene.endf is actually a compound SAB...
    if ZA==1002:
        return 'D'
    if ZA in (133,134):
        return 'H'#solid or liquid methane (CH4 happens to have mass 12+1*4=16 close to oxygen, so can't guess from mass)
    if ZA in (112,113):
        return 'D'#ortho-D or para-D
    guess_by_mass = guessElementFromMass(mass_amu)
    guess_by_z = elementZToName(Z) if (Z>0 and Z<120) else None
    if guess_by_z:
        if guess_by_z!=guess_by_mass:
            raise RuntimeError(f"Mass and Z guess gives different element name (Z={Z}=>'{guess_by_z}', mass={mass_amu}amu=>'{guess_by_mass}')")
        return guess_by_z
    if Z==0:
        return guess_by_mass
    raise ValueError(f'Could not determine element name from ZA={ZA}, mass={mass_amu}amu')

#Custom bare-bones reader, so we can extract the data header as-is in the file:
def parse_endf6(src,name=None,selectfct=None):
    """Parse TSL data from ENDF-6 file according to ENDF manual. Select function can
    used to determine which mat,mf,mt lines to read. E.g. "selectfct=lambda
    mat,mf,mt:mf==7". Remember to select (mf,mt)=(1,451) for the description at
    the start of each material.

    """
    if isinstance(src,str) or isinstance(src,bytes) or hasattr(src,'__fspath__'):
        p=pathlib.Path(src).expanduser()
        with p.open('rt') as srcstream:
            return parse_endf6(srcstream,name=p.name)
    #according to the manual section 0.6, each line is broken down as follows:
    #Columns 1-66 content, 67-70 MAT number, 71-72 MF, 73-75 MT, 76-80 [optional] NS, 81+ undefined
    #Content can either be one 66 char long text field, or there are 6 fields of 11 chars.

    errprefix = '' if not name else f'in {name}'
    def require(check,errmsg):
        if not check:
            raise ValueError(f'{errprefix}: {errmsg}')

    def parse_endf6_line(line):
        require(len(line)>=75,"line too short")
        content=line[0:66]
        mat=int(line[66:70])
        mf=int(line[70:72])
        mt=int(line[72:75])
        return (content,mat,mf,mt)

    content,mat,mf,mt = parse_endf6_line(next(src))
    require((mat,mf,mt)==(1,0,0),'Not an ENDF-6 file?')
    global_header=content
    d={}#material -> material info
    for line in src:
        #print(l)
        content,mat,mf,mt = parse_endf6_line(line)
        if mat<1:
            continue#Ignore dummy lines (section dividers, etc.)
        if selectfct and not selectfct(mat,mf,mt):
            continue
        if mat not in d:
            d[mat]=dict()
        dd=d[mat]
        if mf not in dd:
            dd[mf]={}
        ddd=dd[mf]
        if mt not in ddd:
            ddd[mt]=[content]
        else:
            ddd[mt]+=[content]
    return dict(global_header=global_header,name=name),d

def extract_and_format_endf_file_header(filename):
    custom_parse = parse_endf6(filename,selectfct=lambda mat,mf,mt : (mf,mt)==(1,451))
    mats=custom_parse[1]
    assert len(mats)==1,("File contains more than one material number which"
                         " is not supported by the present script (please"
                         " contact NCrystal developers if you really"
                         " need this).")
    material_number,material_info = next(iter(mats.items()))
    return dict( filename = custom_parse[0]['name'],
                 material_number = material_number,
                 global_header = custom_parse[0]['global_header'],
                 infosection=material_info[1][451] )

def parse_endf_file(filename):
    print(f"Attempting to load ENDF file {filename}...")
    result={}
    warnings_=[]
    pyne = import_pyne()
    data=pyne.endf.Evaluation(filename)
    data.read()
    result['input_file']=pathlib.Path(filename)
    result['pynedata']=data
    print('Performing a few sanity checks...')
    if not hasattr(data,'thermal_inelastic') or 'B' not in data.thermal_inelastic or 'scattering_law' not in data.thermal_inelastic:
        raise RuntimeError('File does not appear to contain any inelastic neutron scattering kernel.')

    ti=data.thermal_inelastic

    if ti['ln(S)']:
        raise RuntimeError('The ln(S) flag is enabled. For the reasons given in the ENDF manual just before'
                           'section 7.4.1, it is not safe to simple convert these to S-values. We need ncmat to support a logsab_scaled option first.'
                           ' Another (easier?) alternative would be to use mpmath module here and calculate unscaled S-values and produce an ncmat file with S=...')

    assert ti['temperature_used'] in ('actual','0.0253 eV')
    if ti['temperature_used'].strip()=='actual':
        alphabetagrid_T0=None
    else:
        alphabetagrid_T0 = float(ti['temperature_used'].split('eV',1)[0].strip())/nc_constants.constant_boltzmann#Must stretch betagrid values relative to this

    #print(data.target)
    element_name=guessElementName(data.target['ZA'],data.target['mass'])
    zsymam=data.target['zsymam'].strip()
    if zsymam not in ['s-ch4','l-ch4','ortho-d','para-d','BENZINE'] and 'graphit' not in zsymam.lower():
        assert element_name in zsymam,f"{element_name} not in {zsymam}"
    result['element_name_principal']=element_name

    lenB = len(ti['B'])
    B  = [None] + list(ti['B'])#B-array, with fortran indexing (for easier reference with ENDF manual):
    A0 = float(B[3]) #Needed to unscale
    if B[1]==0.0:
        raise RuntimeError('Material has no principal scatterers - thus no S(alpha,beta)')
    suggested_emax = B[4]#In principle... but in practice not sure if this holds! Also, data.info['energy_max'] gives a different number (5.0)?!?
    count_principal = B[6]
    num_non_principal = ti['num_non_principal']

    non_principal_data = []
    for i in range(num_non_principal):
        offset = 6*(i+1)
        assert lenB >= 6*(i+2)
        a1=B[offset+1]#0.0: SCT (short collision time approx), 1.0: free-gas, 2.0: diffusive motion.
        if a1!=1.0:
            warnings_+=['WARNING: A non-free-gas model is proposed for the first non-principal atom']
            warnings_+=['         in the original file. This is currently not handled by NCrystal.']
            print(warnings_[-2])
            print(warnings_[-1])
        effective_mass = B[offset+3]#mass in units of neutron mass
        count = B[offset+6]#count per molecule or unit cell
        theElementName = guessElementFromMass(effective_mass*nc_constants.const_neutron_mass_amu)
        non_principal_data+=[[theElementName,count,effective_mass]]

    count_total = count_principal + sum(count_ for elemname,count_,effmass in non_principal_data)
    def format_fraction(c,ctot):
        return f'{c:.14g}/{ctot:.14g}' if c!=ctot else '1'
    for e in non_principal_data:
        elementName,thecount,effective_mass=e
        e += [format_fraction(thecount,count_total)]

    header = extract_and_format_endf_file_header(filename)
    result['header']=header
    result['warnings']=warnings_
    result['non_principal_data']=non_principal_data
    result['count_total']=count_total
    result['count_principal']=count_principal
    result['fraction_str_principal']=format_fraction(count_principal,count_total)
    result['suggested_emax']=suggested_emax

    ndatablock = len(ti['teff'].x)
    result_datablocks=[]
    result['result_datablocks'] = result_datablocks
    for idatablock in range(ndatablock):
        result_datablock={}
        result_datablocks += [result_datablock]
        temperature=ti['teff'].x[idatablock]
        temperature_effective=ti['teff'].y[idatablock]
        result_datablock['T']=temperature
        result_datablock['Teff']=temperature_effective
        sabt=ti['scattering_law'][:,:,idatablock] #indexing t,a,b
        sab=sabt.T
        result_datablock['sab']=sab
        #We must .copy() alpha/beta grids due to '*=' operators below:
        alphagrid = ti['alpha'].copy()*A0#undo funny ENDF alpha definition which divides by A0
        betagrid = ti['beta'].copy()
        if alphabetagrid_T0 is not None:
            #Undo funny ENDF scaling of grids at each temperature
            betagrid *= (alphabetagrid_T0/temperature)
            alphagrid *= (alphabetagrid_T0/temperature)
        result_datablock['alphagrid']=alphagrid.copy()
        result_datablock['betagrid']=betagrid.copy()

        assert ti['symmetric'] == bool(betagrid[0]>=0.0),"inconsistencies detected"

    return result

def format_endf_block_as_ncmatdyninfo_for_principal_element(parsed_endf_data,temperature,fraction_str=None):
    elem_name = parsed_endf_data["element_name_principal"]

    #Find block by temperature:
    temperature_exact,block_idx =  list(sorted((abs(temperature-_["T"]),_["T"],idx) for idx,_ in enumerate(parsed_endf_data['result_datablocks'])))[0][1:]
    assert abs(temperature_exact-temperature)<1e-6
    temperature=temperature_exact

    res =   '@DYNINFO\n'
    res += f"  # Scattering kernel for {elem_name} @ {temperature}K extracted from {parsed_endf_data['input_file'].name}\n"
    for w in parsed_endf_data['warnings']:
        res += f'\n  #  ----> {w}'
    res += '  #\n'
    res += '  # For reference the description embedded in input ENDF file (section MF1,MT451) is repeated here:\n'
    res += '  #\n'
    header=parsed_endf_data['header']
    for line in [header['global_header']]+header['infosection']:
        res+=f'  #     --->{(" "+line).rstrip()}\n'

    datablock = parsed_endf_data['result_datablocks'][block_idx]
    res+=f'  element  {elem_name}\n'
    if fraction_str is None:
        fraction_str = parsed_endf_data["fraction_str_principal"]
    else:
        assert parsed_endf_data["fraction_str_principal"]=='1'
    res+=f'  fraction {fraction_str}\n'
    res +='  type     scatknl\n'
    res+=f'  temperature {temperature:.10} #NB: ENDF file specified "effective temperature" as {datablock["Teff"]:.10}K\n'
    #NB: Not suggesting an emax value, since it seems to be unreliable! Same for B[4] from above...
    #res+=f"  egrid {parsed_endf_data['pynedata'].info['energy_max']}#Value (in eV) as specified in source ENDF file.\n"#argh... e.g. D2O@300K shows this can't be trusted!!!
    res += nccore.formatVectorForNCMAT('alphagrid',datablock['alphagrid'])
    res += nccore.formatVectorForNCMAT('betagrid',datablock['betagrid'])
    res += nccore.formatVectorForNCMAT('sab_scaled',datablock['sab'])
    return res

def parseArgs( progname, arglist, return_parser=False ):
    descr="""
Converts neutron scattering kernels ("thermal scattering laws, tsl") found in
ENDF files to NCMAT format (specifically @DYNINFO sections of NCMAT files).

Note that this script exists for the benefit of experts only. Most users are
recommended to simply download and use one of the pre-converted files shipped
with NCrystal or found at:

    https://github.com/mctools/ncrystal-extra/tree/HEAD/data

Only information for NCMAT @DYNINFO sections are extracted from ENDF files, so
it is recommended to provide the initial parts of the target NCMAT file in a
separate file specified via the --filehdr parameter. That file should include
both initial comments for the file as well as either @DENSITY or
@CELL/@SPACEGROUP/@ATOMPOSITIONS/@DEBYETEMPERATUE sections, depending on whether
or not the material is crystalline. Note that comments describing the origin of
the @DYNINFO sections extracted from the ENDF files will be inserted in the
relevant @DYNINFO sections. The --filehdr should include the string
'<<STDNOTICE>>' on a separate line, which will be expanded to a notice about
availability of files for other temperatures.

If the --filehdr argument is not provided, a dummy header will be
automatically generated, with a fake @DENSITY section, using the density of
water (1g/cm3). This allows the file to be technically valid, but will of
course result in incorrect physics in case the file is used for simulations.

As each ENDF file contains a list of temperatures for which scattering kernels
are available, one NCMAT file will be generated for each such temperature. The
resulting filename will be <basename>_T<temperature>K.ncmat, where <basename>
will default to the input filename unless --outbn is provided.

The script internally relies on the Python Nuclear Engineering Toolkit
(http://pyne.io/) for ENDF parsing, so this must be installed on the system.

Input ENDF files can for instance be downloaded from:

   https://www.nndc.bnl.gov/endf-b8.0/download.html

(download and open zip-file from the "Thermal Neutron Scattering Sublibrary").

"""

    import textwrap
    from argparse import RawTextHelpFormatter

    def wrap(t):
        return textwrap.fill(t,width=60)

    parser = create_ArgumentParser( prog = progname,
                                    description=descr,
                                    formatter_class=RawTextHelpFormatter)
    parser.add_argument('ENDFFILE',type=str,
                        help="Primary ENDF file with tsl data.")
    parser.add_argument('--secondary',metavar='ENDFFILE2:fraction',
                        type=str,
                        help=wrap("ENDF file for secondary element, "
                                  "along with the fraction of that element"))
    parser.add_argument("--filehdr",type=str,
                        help=wrap("File containing initial part"
                                  " of generated files."))
    parser.add_argument("--outbn",type=str,
                        help=wrap("Basename of generated ncmat files."))
    parser.add_argument("--ignoretemp",type=float,nargs='+',metavar='T',
                        help=wrap("If some temperature blocks in input should"
                                  " be ignored, provide the temperature values"
                                  " here (kelvin)."))

    def to_path(parser,fn):
        _ = pathlib.Path(fn)
        if not _.exists():
            parser.error(f'File not found: {_}')
        return _


    if return_parser:
        return parser

    args = parser.parse_args(arglist)

    args.ENDFFILE = to_path(parser,args.ENDFFILE)
    args.fraction1='1'
    args.fraction2=None
    if not args.secondary:
        args.ENDFFILE2=None
    else:
        if not args.secondary.count(':')==1:
            parser.error('Argument to --secondary must contain exactly one semicolon (:)')
        sn,fr2 = args.secondary.split(':',1)
        args.secondary = True
        if fr2.count('/')==1:
            a,b = [int(e) for e in fr2.split('/')]
            if b<=0 or not ( 0 < a < b ):
                parser.error('Invalid fraction specified in --secondary')
            args.fraction1=f'{b-a}/{b}'
            args.fraction2=f'{a}/{b}'
        else:
            args.fraction2=fr2
            if not ( 0.0 < float(fr2) < 1.0):
                parser.error('Invalid fraction specified in --secondary')
            args.fraction1=f'{1.0-float(fr2):10}'
        args.ENDFFILE2 = to_path(parser,sn)
    if args.filehdr:
        _filehdrtxt=to_path(parser,args.filehdr).read_text()
    else:
        _filehdrtxt=f"""
NCMAT v2
#
# Material with scattering kernel extracted from ENDF file.
#
<<STDNOTICE>>
#
@DENSITY
  1 g_per_cm3 #FIX{""}ME. Dummy number!!! (update or add unit cell sections and then remove @DENSITY section)
"""
    args.filehdr = [line.rstrip() for line in _filehdrtxt.strip().splitlines()]

    if not any('<<STDNOTICE>>' in line for line in args.filehdr):
        parser.error('ERROR: File specified with --filehdr should contain the string "<<STDNOTICE>>" on a separate line (just before the first data section)')

    if not args.outbn:
        if args.ENDFFILE2:
            parser.error('--outbn is required when providing two input files')
        else:
            args.outbn=args.ENDFFILE.name
            if args.outbn.endswith('.endf'):
                args.outbn=args.outbn[0:-5]
        assert args.outbn

    return args

def create_argparser_for_sphinx( progname ):
    return parseArgs(progname,[],return_parser=True)

@cli_entry_point
def main( progname, arglist ):
    args = parseArgs( progname, arglist )

    p1=parse_endf_file(args.ENDFFILE)
    p2=parse_endf_file(args.ENDFFILE2) if args.ENDFFILE2 else None

    temperatures1 = [ _["T"] for _ in p1['result_datablocks'] ]
    if p2:
        temperatures2 = [ _["T"] for _ in p2['result_datablocks'] ]
        temperatures_combined = [ t1 for t1 in temperatures1 if min(abs(t1-t2) for t2 in temperatures2)<1e-5 ]
        if len(temperatures_combined)<len(temperatures1):
            warn(f"{len(temperatures1)-len(temperatures_combined)} temperature blocks in {args.ENDFFILE} could not be used due to missing data in other file!")
            if len(temperatures_combined)<len(temperatures2):
                warn(f"{len(temperatures2)-len(temperatures_combined)} temperature blocks in {args.ENDFFILE} could not be used due to missing data in other file!")
    else:
        temperatures_combined = temperatures1

    temperatures_combined=sorted(list(set(temperatures_combined)))
    for t in (args.ignoretemp or []):
        removed=False
        for _ in temperatures_combined:
            if abs(t-_)<1e-6:
                temperatures_combined.remove(_)
                print (f"Ignoring T={_}K as requested")
                removed=True
                break
        if not removed:
            raise RuntimeError(f'Asked to ignore T={t}K, but no such temperature was found in the input')


    if p2 and (p1['non_principal_data'] or p2['non_principal_data']):
        raise RuntimeError("When providing two ENDF files as input both must files must be without non-principal"
                           " elements! If you know how to handle this, You can try to convert separately and combine the"
                           " resulting .ncmat files manually.")

    for t in temperatures_combined:
        fn=pathlib.Path(f'{args.outbn}_T{t}K.ncmat')
        with fn.open('wt') as fh:
            print(f'   -> Writing {fn}')
            fh.write('NCMAT v2\n')
            stdnotice = f'#\n# Notice: This NCMAT file is valid at T={t}K only.'
            if len(temperatures_combined)>1:
                stdnotice+=' Other files alternatively provide\n'
                stdnotice+='# the same material at temperatures:\n#\n'
                nperline=7
                t_to_write=list(_ for _ in temperatures_combined if _!=t)
                for i in range(0,len(t_to_write),nperline):
                    stdnotice += ('#       '+' '.join((f'{_}K' for _ in t_to_write[i:i+nperline]))+'\n')
            for line in args.filehdr:
                if line.startswith('NCMAT '):
                    continue
                if stdnotice and '<<STDNOTICE>>' in line:
                    line=line.replace('<<STDNOTICE>>',stdnotice)
                    stdnotice=''
                fh.write(line if line.endswith('\n') else f'{line}\n')
            if stdnotice:
                fh.write(stdnotice)
            for elementName,count_,effmass,fraction_str in p1['non_principal_data']:
                assert not p2
                fh.write('@DYNINFO\n')
                fh.write(f'  element  {elementName}\n')
                fh.write(f'  fraction {fraction_str}\n')
                fh.write('  type     freegas\n')
            fh.write(format_endf_block_as_ncmatdyninfo_for_principal_element(p1,t,
                                                                             args.fraction1 if p2 else None))
            if p2:
                fh.write(format_endf_block_as_ncmatdyninfo_for_principal_element(p2,t,args.fraction2))
        print('   -> Testing that NCrystal can load this file')
        nccore.createScatter(f'{fn};dcutoff=0.8')
    print("All done.")
