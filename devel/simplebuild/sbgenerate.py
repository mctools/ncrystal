import pathlib
import sys
d = pathlib.Path(__file__).resolve().absolute().parent / 'pypath'
assert d.is_dir()
sys.path.insert(0,str(d))
import sbgen.main
sbgen.main.main()

