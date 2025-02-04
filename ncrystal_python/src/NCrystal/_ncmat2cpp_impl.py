
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

#NB: Important that this script can work WITHOUT any imports (std modules
#allowed, but other NCrystal modules are not), since the file is used as a
#script from the NCrystal CMake code to embed the data library.

_default_regfctname = ( 'NCrystal::registerInMemoryStaticFileData'
                        '(const std::string&,const char*)' )

def _find_data_list( keys, run_standalone = False ):
    if not keys:
        return True, []
    if ( isinstance(keys,list)
         and isinstance(keys[0],dict)
         and set(keys[0].keys())==set(['name',
                                       'read_text_function',
                                       'read_bytes_function']) ):
        #keys is already a datalist:
        return True, keys

    datalist = []
    bns=set()
    for f in keys:
        data = _find_data(f,run_standalone)
        if not data:
            return False, 'File not found: %s'%f
        assert data['name'] is not None
        if data['name'] in bns:
            return False, 'Name not unique in list: %s'%f
        datalist.append(data)
        bns.add(data['name'])
    datalist.sort(key = lambda d : d['name'])
    return True, datalist

def _find_data( key, run_standalone ):
    #Returns None if data is not found, otherwise a dictionary with name and
    #read functions.
    if run_standalone:
        #Not allowed to do any imports apart from python stdlib imports. So we
        #rely exclusively on finding files on disk.
        import pathlib
        f = pathlib.Path(key)
        if not f.exists():
            return None
        return dict( name = pathlib.Path(f).name,
                     read_text_function = f.read_text,
                     read_bytes_function = f.read_bytes )
    from .misc import AnyTextData
    from .core import createTextData,TextData
    if isinstance(key,AnyTextData):
        td = key
        is_textdata = False
    elif isinstance(key,TextData):
        td = key
        is_textdata = True
    elif '\n' not in key:
        td = createTextData(key)
        is_textdata = True
    else:
        td = AnyTextData(key)
        is_textdata = False
    name = td.dataSourceName if is_textdata else td.name

    if name is None:
        from .exceptions import NCBadInput
        raise NCBadInput( 'Can not accept unnamed text data when'
                          ' converting to C++ code, since a string'
                          ' key will be needed in the code.' )
    def td_read_text():
        return td.rawData if is_textdata else td.content
    def td_read_bytes():
        return td_read_text().encode('utf-8')
    return dict( name = name,
                 read_text_function = td_read_text,
                 read_bytes_function = td_read_bytes )

def parseArgs( progname, arglist, *, return_parser = False ):

    #Hidden option used by CMake:
    run_as_standalone_script = False
    while arglist and '--runasstandalonescript' in arglist:
        run_as_standalone_script = True
        arglist.remove('--runasstandalonescript')

    if not run_as_standalone_script:
        from ._cliimpl import create_ArgumentParser
    else:
        def create_ArgumentParser( *a, **kw ):
            import argparse
            return argparse.ArgumentParser(*a,**kw)


    #NOTE: Keep this description in sync with the doc-string in the Python API
    #function in ncmat2cpp.py:
    descr="""

Script which can be used to embed the content of .ncmat files (or actually any
ASCII/UTF8 excoded text files) directly into a C++ library. It does so by
reading the .ncmat files and creating C++ code which keeps the contents of the
files in static strings, and registers those strings with NCrystal, using the
original filename as key. Naturally, those file must be compiled along with the
rest of the C++ library, and the enclosing function must be invoked.

"""
    parser = create_ArgumentParser(prog = progname,
                                   description=descr)
    parser.add_argument( 'FILE', type=str, nargs='+',
                         help=( "One or more NCMAT (or other text) files. They"
                                " will be registered with a key equal to their"
                                " filename (without preceding directory name),"
                                " which must therefore be unique in the list") )
    parser.add_argument('--compact','-c', action='store_true',
                        help=("Set this option to strip all comments and"
                              " excess whitespace from the file data before"
                              " embedding as C++ code."))
    parser.add_argument('--validate','-v', action='store_true',
                        help="""If specified input files will be validated by confirming that they can be
                        loaded with NCrystal. For this to work, the NCrystal python module must be
                        available.""")
    parser.add_argument("--name",'-n',default='registerNCMATData',type=str,
                        help="""Name of C++ function to create which must be called in order to register the
                        data with NCrystal. If desired, it can contain namespace(s), e.g. a value of
                        "MyNameSpace::myFunction" will create a function "void myFunction()" in the
                        namespace "MyNameSpace.""")
    parser.add_argument("--regfctname",default=_default_regfctname,type=str,
                        help="""Name of C++ function used to register string objects with NCrystal.""")
    parser.add_argument("--include",nargs='+',type=str,action='append',
                        help="""One or more extra include statements for the top of the file. The file
                                NCrystal/core/NCDefs.hh will always be included by default (prevent this by
                                adding a special entry 'no-ncrystal-includes' to the list).""")
    parser.add_argument("--width",'-w',type=int,default=80,
                        help=("Wrap C++ code at this column width. Ignored"
                              " for text files <~60KB unless running with"
                              " --compact."))
    parser.add_argument('--outfile','-o',type=str,default='stdout',
                        help="Name of output file (default: stdout)")

    if return_parser:
        return parser

    args=parser.parse_args( arglist )
    args.run_as_standalone_script = run_as_standalone_script

    if not args.name or ' ' in args.name:
        parser.error('Invalid C++ function name provided to --name')

    ok, datalist_or_errmsg = _find_data_list( args.FILE,
                                              run_as_standalone_script )
    if not ok:
        parser.error(datalist_or_errmsg)
    args.files = datalist_or_errmsg
    args.FILE=None

    wmin=30
    wmax=999999
    if args.width>wmax:
        args.width=wmax
    if args.width < wmin:
        parser.error('Out of range value of --width (must be at least %i)'%wmin)

    #flatten args.include, so we get single list with 3 elements from:
    #        --inc 'foobla.hh' --inc '<vector>' Bla/Bla.hh
    args.include = [item for sublist in (args.include or []) for item in sublist]

    if args.run_as_standalone_script and args.validate:
        parser.error('Do not use --validate with hidden'
                     ' --runasstandalonescript option')

    return args

_sys_print = print
def files2cppcode(infiles,
                  *,outfile,
                  cppfunctionname='registerData',
                  compact=False,
                  width=140,
                  validate=False,
                  extra_includes=None,
                  regfctname=None,
                  quiet = False,
                  run_standalone = False
                  ):
    #NOTE: This function is called both from the CLI script (this file) and the
    #Python API function in ncmat2cpp.py

    #We would like to support AnyTextData, but since we must try to not do any
    #imports when this is used from CLI (to use it in the NCrystal CMake cfg
    #step to embed the standard lib), we should be a bit more careful:
    orig_print_fct = _sys_print
    if quiet:
        def print(*a,**kw):
            pass
    elif run_standalone:
        print = _sys_print
    else:
        from ._cliimpl import print
        orig_print_fct = print

    if regfctname is None:
        regfctname = _default_regfctname
    if regfctname.startswith('NCrystal::'):
        regfctname = 'NCRYSTAL_NAMESPACE::' + regfctname[len('NCrystal::'):]

    out=['// Code automatically generated by NCrystal/ncmat2cpp','']
    if extra_includes is None:
        extra_includes = []
    if 'no-ncrystal-includes' in extra_includes:
        extra_includes = [ e for e in extra_includes
                           if e!='no-ncrystal-includes' ]
    else:
        extra_includes.insert(0,'NCrystal/core/NCDefs.hh')

    for inc in extra_includes:
        if inc.startswith('#include'):
            out.append(inc)
        elif inc.startswith('<'):
            out.append('#include %s'%inc)
        else:
            out.append('#include "%s"'%inc)
    out.append('')

    large_files = False

    def fwddeclare(out,fctname,args_str=''):
        if '(' not in fctname:
            fctname+='()'
        _ = fctname.split('(',1)[0].split('::')
        namespaces,justname = _[0:-1],_[-1]
        if namespaces and namespaces[0]=='NCrystal':
            namespaces[0]='NCRYSTAL_NAMESPACE'
        argssignature='('+fctname.split('(',1)[1]
        tmp=''
        for ns in namespaces:
            tmp += 'namespace %s { '%ns
        tmp += 'void %s%s;'%(justname,argssignature)
        tmp += ' }'*len(namespaces)
        out+=[tmp]
        out+=['']

    fwddeclare(out,regfctname)
    if '::' in cppfunctionname:
        fwddeclare(out,cppfunctionname)

    out += [ 'void %s()'%cppfunctionname, '{' ]
    prefix='  '
    seen = set()

    ok, datalist_or_errmsg = _find_data_list( infiles )
    if not ok:
        if run_standalone:
            raise RuntimeError(datalist_or_errmsg)
        else:
            from .exceptions import NCBadInput
            raise NCBadInput(datalist_or_errmsg)
    else:
        datalist = datalist_or_errmsg

    for data in datalist:
        #for p in [pathlib.Path(f) for f in infiles]:
        fn= data['name']
        print(f"ncmat2cpp : Processing {fn}")
        assert fn not in seen, "ERROR: Multiple files in input named: %s"%fn
        seen.add(fn)

        def fmtline(line):
            if compact:
                ncmatcfg=None
                if 'NCRYSTALMATCFG[' in line:
                    #special case: preserve NCRYSTALMATCFG
                    _=line.split('NCRYSTALMATCFG[',1)[1]
                    errmsg = None
                    if 'NCRYSTALMATCFG' in _:
                        errmsg = "multiple NCRYSTALMATCFG entries in a single line"
                    elif ']' not in _:
                        errmsg = "NCRYSTALMATCFG[ entry without closing ] bracket"
                    if errmsg:
                        if run_standalone:
                            raise RuntimeError(errmsg)
                        else:
                            from .exceptions import NCBadInput
                            raise NCBadInput(errmsg)
                    ncmatcfg=_.split(']',1)[0].strip()
                    ncmatcfg.encode('utf8')#Just a check
                line=' '.join(line.split('#',1)[0].split())
                line.encode('ascii')#Just a check
                if ncmatcfg:
                    line += '#NCRYSTALMATCFG[%s]'%ncmatcfg
                if not line:
                    return ''
            else:
                line.encode('utf8')#Just a check
            return line.replace('"',r'\"')+r'\n'

        def as_c_str(strdata):
            try:
                strdata.encode('ascii')
                return '"%s"'%strdata
            except UnicodeEncodeError:
                pass
            try:
                strdata.encode('utf8')
                return 'u8"%s"'%strdata
            except UnicodeEncodeError:
                raise SystemExit('Invalid encoding encountered in'
                                 ' input (must be ASCII or UTF8)')


        raw_text_data = data['read_text_function']()

        if validate:
            assert not run_standalone,"standalone mode prevents --validate"
            print("Trying to validate: %s"%fn)
            from .misc import MaterialSource
            MaterialSource( raw_text_data ).load( doInfo = True,
                                                  doScatter = False,
                                                  doAbsorption = False )
            print('  -> OK')

        lines = list(raw_text_data.splitlines())
        assert lines,"file was empty: %s"%fn

        #string literals have a limit of 65K in the standard. For such large
        #files we must embed contents in const std::array<std::uint8_t, 12>
        #arrays.
        is_large = ( sum(len(as_c_str(fmtline(line))) for line in lines) > 60000 )

        if is_large:
            large_files=True

        out+= [prefix+"{"]
        out+= [prefix+"  // File %s%s%s"%(fn,
                                          (' (compact form without comments)' if compact else ''),
                                          (' (too large for string literals)' if is_large else ''))]

        if is_large:
            #NB: This could be used for non-text-data as well!
            raw_data_bytes = data['read_bytes_function']()
            out += [ prefix+'  static const std::array<std::uint8_t,%i> rawdata {'%(len(raw_data_bytes)+1)]
            n = len(raw_data_bytes)
            n_on_current_line = 0
            delim,currentline='',''
            _prefstr = prefix+'    '
            ndatawidth = width-len(_prefstr)
            for c in raw_data_bytes:
                ++n_on_current_line
                currentline += delim + str(c)
                delim = ','
                if len(currentline)>=ndatawidth:
                    out += [_prefstr + currentline+',']
                    delim=''
                    currentline = ''
            eol=',0};'
            if currentline:
                out += [_prefstr + currentline+eol]
            else:
                out[-1] = out[-1][:-1]+eol
            out += ['    const char * textdata = (const char*)(&rawdata[0]);']
        else:
            #count_all_entries = [0]
            def linepattern(strdata):
                return '    %s'%as_c_str(strdata)
            #fix_nentries_iout = len(out)
            out += [ prefix+'  const char * textdata =']
            if not compact:
                for line in lines:
                    out+= [prefix+linepattern(fmtline(line))]
            else:
                alldata=''
                for line in lines:
                    alldata += fmtline(line)
                while alldata:
                    n = max(10,width-(len(linepattern(''))))
                    if alldata[n-1:n+1]==r'\n':
                        n+=1#dont break up '\n' entries
                    out += [ prefix+linepattern(alldata[0:n])]
                    alldata = alldata[n:]
            out[-1]+=';'
        out+= [prefix+"  ::%s(\"%s\",textdata);"%(regfctname.split('(')[0],fn)]
        out+= [prefix+"}"]
    out+= ['}','']

    if large_files:
        out.insert(1,'#include <cstdint>')
        out.insert(1,'#include <array>')

    out = '\n'.join(out)
    if outfile is None:
        pass
    elif outfile == 'stdout':
        for line in out.splitlines():
            orig_print_fct(line)
    else:
        import pathlib
        of=pathlib.Path(outfile)
        if run_standalone:
            of.write_text(out,encoding='utf8')
        else:
            from ._common import write_text
            write_text(of,out)
        print('Wrote: %s'%of)
    return out

#Sphinx doc function. Signature always the following:
def create_argparser_for_sphinx( progname ):
    return parseArgs([progname],return_parser=True)

def main( progname, arglist ):
    args = parseArgs( progname, arglist )
    return files2cppcode( args.files,
                          outfile = args.outfile,
                          cppfunctionname = args.name,
                          compact = args.compact,
                          width = args.width,
                          extra_includes = args.include,
                          validate = args.validate,
                          regfctname = args.regfctname,
                          quiet = False,
                          run_standalone = args.run_as_standalone_script )

if __name__ == '__main__':
    #Running from CMake code to embed the standard data library.
    import sys
    import os
    progname = os.path.basename(sys.argv[0])
    arglist = sys.argv[1:] + ['--runasstandalonescript']
    main( progname, arglist )
