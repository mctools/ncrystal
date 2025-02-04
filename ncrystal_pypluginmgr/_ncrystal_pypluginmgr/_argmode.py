
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

_usage_txt = """

The ncrystal-pluginmanager command is invoked by the NCrystal library in order
to discover any plugins available in the environment.

Additionally it can be used to load and test a given dynamic plugin, through the
--test flag. For instance, once might test the plugin named "SomePlugin" with
the command:

ncrystal-pluginmanager --test SomePlugin

This requires ncrystal-python to be installed.

This is merely a convenience, since a similar test can be enabled by setting the
environment variables:

   NCRYSTAL_PLUGIN_RUNTESTS=1
   NCRYSTAL_REQUIRED_PLUGINS=SomePlugin

and then loading the NCrystal library and calling ensurePluginsLoaded().
"""

def main():
    import sys
    import os
    args = sys.argv[1:]
    if not args:
        from .cli import main
        main()
        return
    helpmode = '-h' in args or '--help' in args
    testmode = '--test' in args
    testlist = set([a for a in args if a != '--test'])
    argerror = False
    if testmode:
        if not testlist or any(a.startswith('-') for a in testlist ):
            argerror = True

    if helpmode or argerror:
        assert os.path.basename( sys.argv[0] ) == 'ncrystal-pluginmanager'
        print(_usage_txt.strip()+'\n')
        if argerror:
            raise SystemExit(1)
        return

    os.environ['NCRYSTAL_PLUGIN_RUNTESTS'] = '1'
    os.environ['NCRYSTAL_REQUIRED_PLUGINS'] = ':'.join(sorted(testlist))
    os.environ['NCRYSTAL_SLIMPYINIT'] = '1'

    try:
        import NCrystal # noqa F401
    except ModuleNotFoundError:
        raise SystemExit('ERROR: Can not use --test without'
                         ' the ncrystal-python package installed')
    import NCrystal.plugins as ncp
    #Something that triggers ensurePluginsLoaded():
    pluglist = ncp.browsePlugins(dump=False)
    #Are we still here? In that case there must have been no issues. As a sanity
    #check we verify the plugins are also in the returned pluglist:
    missing = testlist - set(pluginname for pluginname,_,_ in pluglist)
    if missing:
        raise SystemExit('ERROR: plugin missing (weird that it was'
                         ' not detected earlier: %s'%missing)

    import NCrystal.datasrc as nds
    import NCrystal.core as nccore
    files_to_test = []
    for f in nds.browseFiles(factory='plugins'):
        if f.name.split('/')[0] in testlist:
            files_to_test.append( f.fullKey )
    for f in files_to_test:
        print(f'Testing load of "{f}"')
        for fctname in ['createTextData','createInfo',
                        'createScatter','createAbsorption']:
            print(f'  -> {fctname}',flush=True)
            getattr(nccore,fctname)( f )

    print("All ok",flush=True)
