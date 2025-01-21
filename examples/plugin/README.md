Example dynamic plugin for NCrystal
===================================

This directory serves as a template which can be used as a starting point for
creating a custom NCrystal plugin. So copy the files and folders in this
directory to a new repository, decide on a name of the plugin, and you can start
modifying as needed.

First of all, decide on a name. In the example files, the name of the plugin is
`DummyPlugin`, but this should of course be updated to something more
appropriate. Once a suitable name has been chosen, be sure to update this files
`ncplugin_name.txt` and `pyproject.toml`, and replace `DummyPlugin` with the
name you have chosen. Additionally, the `@CUSTOM_DUMMYPLUGIN` section in any
NCMAT files in the `data/` directory should also match the name of the plugin.

After that, you should update the files in the `src/` folder to reflect the
actual physics of your plugin. In the example, the plugin allows materials to
replace the usual incoherent-elastic models of NCrystal with a custom silly
model in which incoherent-elastic physics has a cross section which is a
step-function in energy (or wavelength), and all scatterings are isotropic and
elastic. The plugin is activated for NCMAT files which have a
@CUSTOM_DUMMYPLUGIN section with two numbers (the cross section below the
wavelength theshold in barn, and the wavelength threshold in angstrom). The
primary files to modify are:

  * `src/NCPhysicsModel.cc`: This is where the actual physics model is
    implemented.

  * `src/NCPluginFactory.cc`: This is where it is defined how the new physics
    model is injected along-side the other existing physics in NCrystal.

  * `src/NCTestPlugin.cc`: This contains a small test function which exercises
    the plugin. This test function will be automatically called by NCrystal if
    the environment variable `NCRYSTAL_PLUGIN_RUNTESTS` is set to `1`.

To install and activate the plugin, one can either build it manually via CMake,
and then afterwards add the resulting shared library to the
`NCRYSTAL_PLUGIN_LIST` variable. However, a more convenient option is to simply
use Python build tools to install the plugin (presumably in a Python or Conda
environment), in which case NCrystal will automatically discover the presence of
the plugin and load it. So installing the plugin might be as simple as: `pip
install <plugindir>` where `<plugindir>` is the directory in which the
`pyproject.toml` and `CMakeLists.txt` files of the plugin resides.

Finally, any .ncmat file placed inside the `data/` folder will be included in
the plugin and automatically be made available to NCrystal users as well,
