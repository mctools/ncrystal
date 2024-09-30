
from NCTestUtils.loadlib import Lib
lib = Lib('testfileutils')
lib.dump()
assert hasattr(lib,'nctest_file_exists')
assert hasattr(lib,'nctest_ncgetcwd')

import pathlib
d = pathlib.Path('.')
(d / 'foo.txt').write_text('foo')

assert lib.nctest_file_exists('./foo.txt')
assert lib.nctest_file_exists('foo.txt')
assert not lib.nctest_file_exists('bar.txt')
assert not lib.nctest_file_exists('./bar.txt')
assert not lib.nctest_file_exists('/some/file/that/does/not/exist')
print("Testing that two files are identical:")
_cwd = lib.nctest_ncgetcwd()
d_abs = d.resolve().absolute()
assert d.samefile(d_abs)
print("  1) %s"%d_abs)
print("  2) %s"%_cwd)
assert d_abs.samefile ( pathlib.Path(_cwd).resolve().absolute() )
assert d_abs.samefile ( _cwd )
assert d.samefile ( _cwd )

#print('cwd',lib.nctest_ncgetcwd())

#FIXME: Much more, including globbing and whatever is likely to cause issues on
#Windows.
