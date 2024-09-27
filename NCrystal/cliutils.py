"""

Module providing access to the ncrystal commandline tools via the Python API.

"""

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

def cli_tool_list( canonical_names = True  ):
    """
    Get list of available NCrystal command line tools. The canonical_names flag
    controls whether to return canonical (ncrystal_ncmat2cpp, nctool, ...)  or
    short (ncmat2cpp, nctool, ...) names.
    """
    import pathlib
    short_names = [ f.name[5:-3] for f in
                    pathlib.Path(__file__).parent.glob('_cli_*.py') ]
    short_names.sort()
    if short_names:
        return short_names
    else:
        return [ _map_shortname_2_canonical_name(sn) for sn in short_names ]

def cli_tool_lookup( name ):
    """
    Try to lookup command line tool name, accepting various aliases, to for
    instance recognise both "ncmat2cpp" and "ncrystal_ncmat2cpp".

    Returns a dictionary with both a 'canonical_name' and 'short_name', where
    the former is the name of the corresponding command-line script, and the
    latter is the shortest name which can be used to identify the tool

    Returns None in case the name could not be resolved to an available tool.
    """
    #Note: We basically have to treat only ncrystal-config and nctool as special
    #cases.
    if name == 'ncrystal-config':
        #Special case:
        #FIXME: We don't actually have such a script here yet!
        short_name = 'config'
    elif name.startswith('ncrystal_'):
        short_name = name[9:]
    else:
        short_name = name
    if short_name not in cli_tool_list( canonical_names=False ):
        return None
    return dict( short_name = short_name,
                 canonical_name = _map_shortname_2_canonical_name(short_name) )


def run( cmdname, *arguments, check = True ):

    """Can be used to invoke ncrystal command-line tools such as nctool,
    ncrystal_ncmat2cpp, ncrystal_hfg2ncmat, etc. directly via the Python API
    without the need for spawning separate subprocesses or actually using a
    shell to execute command scripts. It also makes it easier to invoke the
    tools in a cross platform way, by staying exclusively in the Python API.

    Example:

    run('ncrystal_ncmat2cpp','Al_sg225.ncmat','-o','myal.cpp')

    Is the same as running the command in the terminal (here typical unix shell
    syntax):

    $> ncrystal_ncmat2cpp Al_sg225.ncmat -o myal.cpp

 on the commandline with the provided
    argument list. Returns results object without returncode, stdout, stderr,
    exc_type, and exc_value attributes. In case of no errors, returncode will be
    0 while exc_type and exc_type will both be None.

    FIXME: What about ncrystal-config?

    """
    #FIXME: Test that it gives a good experience in Jupyter

    #FIXME: support + unit test various aliases
    #$> python -mNCrystal ncrystal_ncmat2cpp ...
    #$> python -mNCrystal ncmat2cpp ...
    #$> python -mNCrystal nctool ...
    #$> python -mNCrystal ncrystal_nctool ... (NOPE?)
    #$> python -mNCrystal config ... #Perhaps implement this by shutil.which('ncrystal-config') and then IMPORTING the file directly? And of course, also then checking the version, etc. (we should most likely walk the PATH until we find an ncrystal-config with the correct version). Or we could simply have the ncrystal-config script as a local submodule named _ncrystal_config_copy.py ?
    #$> python -mNCrystal ncrystal-config ...

    #FIXME: The embedded tests could verify that all of the command line script
    #are available. That way, we would not forget to update the conda-forge
    #recipe when adding a script.

    resolved_cmd = cli_tool_lookup( cmdname )
    if resolved_cmd is None:
        from .exceptions import NCBadInput
        raise NCBadInput(f'Command line tool name "{cmdname}" not recognised')
    argv = [resolved_cmd['canonical_name']] + [a for a in arguments]

    import importlib
    clipymodname = '_cli_%s'%resolved_cmd['short_name']
    climod = importlib.import_module(f'..{clipymodname}', __name__)
    assert hasattr(climod,'main')

    #Fixme: unit test the check=False usage!! The code below was not tested
    #much, nor was the return value..

    ec = 0
    res = dict( returncode = 0,
                stdout = '',
                stderr = '',
                exc_type = None,
                exc_value = None )
    from ._common import print, ctxmgr_redirect_argparse_output
    try:
        with ctxmgr_redirect_argparse_output():
            climod.main( argv )
    except SystemExit as e:
        if len(e.args)==1 and isinstance(e.args[0],int):
            ec = e.args[0]
        else:
            ec = 1
        if check and ec != 0:
            raise RuntimeError(f'Command ended with exit code {ec}')
        else:
            res['exc_type'] = SystemExit
            res['exc_value'] = e.args
        res['returncode'] = ec
    except Exception as e:
        if check:
            raise
        res['returncode'] = 1
        res['exc_type'] = e.__class__
        res['exc_value'] = e.args

    class CLIResult:
        pass
    o = CLIResult()
    for k,v in res.items():
        setattr(o,k,v)
    return o

def _map_shortname_2_canonical_name( short_name ):
    return { 'config' : 'ncrystal-config',
             'nctool' : 'nctool' }.get( short_name,
                                        f'ncrystal_{short_name}' )

def _single_integer_or_None( arglist ):
    if len(arglist)==1 and isinstance(arglist[0],int):
        return arglist[0]

def _run_from_main_init():
    import sys
    import textwrap
    from ._common import print
    usagestr = textwrap.dedent("""
    Usage: provide name and arguments of NCrystal commandline-tool to run.

    Available tools:

    <<TOOLLIST>>
    Specifying -h or --help displays this message.""").lstrip()


    is_help_req = len(sys.argv)>1 and sys.argv[0] in ('-h','--help','/?')
    if is_help_req or len(sys.argv)<2 or sys.argv[0].startswith('-'):
        toolliststr = ''.join( f'       {t}\n' for t
                               in cli_tool_list( canonical_names = False ))
        print(usagestr.replace('<<TOOLLIST>>',toolliststr))
        raise SystemExit(0 if is_help_req else 1)

    res = run(sys.argv[1],*sys.argv[2:])
    if res.returncode==0:
        #TODO: Print the output?
        return
    if res['exc_type'] is SystemExit:
        ec = _single_integer_or_None(res['exc_value'])
        if ec is not None:
            #No additional message, just a pure exit
            raise SystemExit(ec)
        errstr = 'Error'
    else:
        errstr = res['exc_type'].__name__
    for a in res['exc_value']:
        print(f"{errstr}: {a}")
    raise SystemExit(res.returncode)

    #FIXME: Make sure we can run via:
    #$> python -mNCrystal ncrystal_ncmat2cpp Al_sg225.ncmat -o myal.cpp
    #Or even just:
    #$> python -mNCrystal ncmat2cpp Al_sg225.ncmat -o myal.cpp

#FIXME: Stop using SystemExit exceptions, but use a custom CLISysExit exception instead? In the actual cmdline usage we can then rethrow those??
