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
    from ._cliimpl import cli_tool_list_impl
    return cli_tool_list_impl( canonical_names = canonical_names )

def cli_tool_lookup( name ):
    """
    Try to lookup command line tool name, accepting various aliases, to for
    instance recognise both "ncmat2cpp" and "ncrystal_ncmat2cpp".

    Returns a dictionary with both a 'canonical_name' and 'short_name', where
    the former is the name of the corresponding command-line script, and the
    latter is the shortest name which can be used to identify the tool

    Returns None in case the name could not be resolved to an available tool.
    """
    from ._cliimpl import cli_tool_lookup_impl
    return cli_tool_lookup_impl( name )

def run( toolname, *arguments ):
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

    #FIXME: support + unit test various methods of invocation:
    #$> python -mNCrystal ncrystal_ncmat2cpp ...
    #$> python -mNCrystal ncmat2cpp ...
    #$> python -mNCrystal nctool ...
    #$> python -mNCrystal ncrystal_nctool ... (NOPE?)

    #FIXME: The embedded tests could verify that all of the command line script
    #are available. That way, we would not forget to update the conda-forge
    #recipe when adding a script.

    from ._cliimpl import _resolve_cmd_and_import_climod
    climod, argv = _resolve_cmd_and_import_climod( toolname, arguments )

    from ._cliimpl import ( ctxmgr_modify_argparse_creation,
                            _cli_call_from_pyapi_ctx )
    try:
        #Fixme: merge the two following context managers?:
        with _cli_call_from_pyapi_ctx():
            with ctxmgr_modify_argparse_creation(exit_on_error = False,
                                                 redirect_stdout_to_stderr = True):
                climod.main( argv )
    except SystemExit as e:
        #Map SystemExit to either a clean return or a RuntimeError.
        if str(e) in ('','0'):
            return#ended OK
        if len(e.args)==1 and isinstance(e.args[0],int):
            ec = e.args[0]
        elif str(e).isdigit():
            ec = str(e)
        else:
            ec = 1
        msg = f'Command ended with exit code {ec}'
        if not str(e).isdigit():
            msg = str(e)
        raise RuntimeError(msg) from e
