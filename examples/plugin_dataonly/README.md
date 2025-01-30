Example plugin for NCrystal providing data files
================================================

This directory shows how a simple static plugin can be used to provide data
files for NCrystal.

Like any plugin, one should first of all, decide on a name. In the example here,
the name of the plugin is `DummyDataPlugin`, but this should of course be
updated to something more appropriate. Once a suitable name has been chosen, be
sure to update this name in the `pyproject.toml` file, in the `project.name`
field (i.e. replacing `DummyDataPlugin` with the name of your
choice). Additionally, one should rename the directory
`src/ncrystal_plugin_DummyDataPlugin`, once again replacing `DummyDataPlugin`
with the name of your choice. The data files should go in a `data/` subdirectory
of that directory (i.e. remove `dummy.ncmat` and put the files you wish to
provide there).

To install and activate the plugin, one must simply use Python build tools to
install the plugin (presumably in a Python or Conda environment), in which case
NCrystal will automatically discover the presence of the plugin and load it. So
installing the plugin might be as simple as: `pip install <plugindir>` where
`<plugindir>` is the directory in which the `pyproject.toml` resides.
