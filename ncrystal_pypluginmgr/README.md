NCrystal python plugin manager
==============================

Package providing an ncrystal-pluginmanager command, which searches for dynamic
NCrystal plugins (modular shared libraries) installed into module directories of
python modules named "ncplugin_*" in plugins/ subdirs.

When an ncrystal-pluginmanager command is available, the NCrystal library itself
will use it to get a list of such plugin files to load. This makes it possible
to simply install an NCrystal plugin via "pip install
<url-to-repo-where-ncrystal-plugin-is-developed>", and then have NCrystal pick
it up automatically. This is meant to be easier than the alternative of
compiling the plugin code manually, and then appending the resulting shared
library to the `NCRYSTAL_PLUGIN_LIST` environment variable.
