Dynamic C++ API example
=======================

Small example showing how a C++ application can use NCrystal without actual
build time dependency on NCrystal. The only requirement is that the project
embeds a copy of the relevant dynamic API interface header from the NCrystal
sources (in this example that is in `src/ncrystal_load.hh`), and is willing to
perform runtime discovery of the NCrystal shared library via the
`ncrystal-config` command, dynamic loading of it via dlopen/LoadLibrary, and
finally calling the special `ncrystal_access_dynamic_api` function.

This might sound slightly complicated, but in this example, all the nitty gritty
is encapsulated inside a set of utility files (`ncrystal_load.hh` and
`ncrystal_load.cc`), and the actual usage in the `main.cc` programme is
extremely simple. The advantage is clear: the application can have optional
support for NCrystal without any need for complicated build-time flags and
dependency management, and existing packages of a particular application can
work with future releases of NCrystal as they come out.

In a way, this provides C++ with a feature akin to what one might use in Python
to support optional dependencies:

```
try:
    import NCrystal
except ModuleNotFoundError
    NCrystal = None
    print("Sorry, you must install NCrystal to use this feature")
```

For this to work for a particular application, there must be a "Dynamic API
header" available in NCrystal with exactly the functionality needed by that
application. With the initial release of this feature in NCrystal 4.1.0, there
is only one such header file available:

`NCrystal/internal/dynapi/NCDynAPI_Type1_v1.hh`

If it is ever discovered that a small tweak is needed in that API, the old file
will NOT be touched, but rather a new file (`NCDynAPI_Type1_v2.hh`) will be made
available as well. If a rather different set of functionality is needed, it
might instead be decided to have a different type of API available
(`NCDynAPI_Type2_v1.hh`). If you have a particular need for functionality in
such a dynamic API, please get in touch with the NCrystal developers.
