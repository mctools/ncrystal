# Make the sbgen module (from ./pypath/sbgen) and the ncrystal_repo_tools module
# (from (../pypath/ncrystal_repo_tools) available:
import pathlib
import sys
d = pathlib.Path(__file__).resolve().absolute().parent
d1 = d / 'pypath'
d2 = d.parent / 'pypath'
assert d1.is_dir()
assert d2.is_dir()
sys.path.insert(0,str(d2))
sys.path.insert(0,str(d1))
#Now launch the main generator entry point:
import sbgen.main # noqa E402
sbgen.main.main()

