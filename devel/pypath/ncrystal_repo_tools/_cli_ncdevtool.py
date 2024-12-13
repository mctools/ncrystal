
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

def get_mode_list():
    import pathlib
    return sorted( m.name[9:-3]
                   for m in
                   pathlib.Path(__file__).parent.glob('_climode_*.py') )

def progname():
    #import os
    #return os.path.basename(sys.argv[0])
    return 'ncdevtool'

def usage():
    theprogname = progname()
    ml = get_mode_list()
    example_mode = 'grep'
    assert example_mode in ml
    modelist_indent = '\n     '
    modelist_str = modelist_indent + f'{modelist_indent}'.join(ml)
    print(f"""Usage:

  $> {theprogname} MODE <mode options>

Available MODEs are:
{modelist_str}

Use the -h or --help flag to get information about the options
available for a given mode. For example:

  $> {theprogname} {example_mode} --help
""")

def import_sibling_module( module_name ):
    import importlib
    pkgarg = __name__
    if pkgarg == '__main__':
        #Make running as python -m <packagename>.<thismodule> work:
        pkgarg = '%s.foo'%__package__
    return importlib.import_module(f'..{module_name}',pkgarg)

def main():
    import sys
    if len(sys.argv)<2 or sys.argv[1] in ('-h','--h','--he','--hel','--help'):
        usage()
        return
    mode = sys.argv[1]
    if mode not in get_mode_list():
        raise SystemExit(f'ERROR: Invalid mode "{mode}". Run without'
                         ' arguments to list available modes.')

    parser = NCDevUtilsArgParser( cliname = progname(),
                                  modename = mode,
                                  args = sys.argv[2:] )

    import_sibling_module('_climode_%s'%mode).main( parser )

class NCDevUtilsArgParser():

    def __init__(self, *, cliname, modename, args ):
        self.__parser = None
        self.__helpw = 59
        self.__descrw = 79
        self.__progname = f'{cliname} {modename}'
        self.__needs_wrap = False
        self.__args = args

    def get_raw_args( self ):
        #For special modes, not using the usual parser
        return self.__args

    def init( self, descr, **kwargs ):
        assert 'prog' not in kwargs
        assert 'descr' not in kwargs
        import argparse
        if 'formatter_class' not in kwargs:
            self.__needs_wrap = True
            import textwrap
            kwargs['formatter_class'] = argparse.RawTextHelpFormatter
            newdescr = ''
            for i, p in enumerate(descr.strip().split('\n\n')):
                if i:
                    newdescr += '\n\n'
                newdescr += textwrap.fill(' '.join(p.strip().split()),
                                          self.__descrw)
            #'`N' and '`' are preserved newlines and spaces respectively.
            descr = newdescr.replace('`N','\n').replace('`',' ')
        self.__parser = argparse.ArgumentParser( prog = self.__progname,
                                                 description=descr,
                                                 **kwargs )

    def __fixhelp( self, msg ):
        assert self.__parser is not None
        if self.__needs_wrap:
            import textwrap
            return textwrap.fill( ' '.join(msg.split()), width=self.__helpw )
        return msg

    def add_argument( self, *args, **kwargs ):
        assert 'help' in kwargs
        kwargs['help'] = self.__fixhelp( kwargs['help'] )
        self.__parser.add_argument(*args, **kwargs)

    def parse_args( self ):
        return self.__parser.parse_args( self.__args )

    def error( self, msg ):
        self.__parser.error(msg)

if __name__ == '__main__':
    main()
