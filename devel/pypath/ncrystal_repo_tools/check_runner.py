def get_available_checks_list():
    import pathlib
    return sorted( [ f.name[7:-3] for f in
                     pathlib.Path(__file__).parent.glob('_check_*.py') ] )

def run_check( name ):
    print()
    print(f">>>>> Running check: {name}")
    print()
    import importlib
    pymodname = '_check_%s'%name
    mod = importlib.import_module(f'..{pymodname}', __name__)
    assert hasattr(mod,'main')
    mod.main()

def run_all_checks():
    for c in get_available_checks_list():
        run_check(c)

# Planned additional tests:
#    * license blurbs
#    * shebangs consistent
#    * no excessive trailing whitespace (end of file, end of line)
#    * no tabs or non-ascii chars in source files.
#    * no include statements that violates the
#    * Ruff checks of all code
#
# We can also then have _tool_bla.py (mostly for accessing via ncdevtool:
#
#    * plot dependency graphs (between components)
#    * build the docs
#    * launch the local cmake based builds into some dir
#    * launch the ctests
#
# CI missing:
#
#    * exercising embeded data with e.g. tests.
#    * cpp check
#    * python coverage testing.
#    * static tests
