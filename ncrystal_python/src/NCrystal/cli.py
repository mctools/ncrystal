
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

"""

Module providing access to the ncrystal commandline tools via the Python API.

"""

def cli_tool_list( canonical_names = True  ):
    """
    Get list of available NCrystal command line tools. The canonical_names flag
    controls whether to return canonical (ncrystal_ncmat2cpp, nctool, ...)  or
    short (ncmat2cpp, nctool, ...) names.
    """
    from ._cliimpl import cli_tool_list_impl
    return cli_tool_list_impl( canonical_names = canonical_names )

def cli_tool_lookup( name ):
    """Try to lookup command line tool name, accepting various aliases, to for
    instance recognise both short names like "ncmat2cpp" and canonical names
    like "ncrystal_ncmat2cpp".

    Returns a dictionary with both a 'canonical_name', 'short_name', and
    'shellcmd'. Here, the 'canonical_name' is the canonical name of the
    corresponding command-line script in standard installations, the
    'short_name' is the shortest name which can be used to identify the tool,
    and 'shellcmd' is the actual name of the shell command in the present
    installation.

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

    """

    from ._cliimpl import _resolve_cmd_and_import_climod
    climod, argv = _resolve_cmd_and_import_climod( toolname, arguments )

    from ._cliimpl import ( ctxmgr_modify_argparse_creation,
                            _cli_call_from_pyapi_ctx )
    try:
        with _cli_call_from_pyapi_ctx():
            with ctxmgr_modify_argparse_creation(exit_on_error = False,
                                                 redirect_stderr_to_stdout = True):
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
