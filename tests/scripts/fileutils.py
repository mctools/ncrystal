
from NCTestUtils.loadlib import Lib
lib = Lib('testfileutils')
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
assert d.samefile ( lib.nctest_ncgetcwd() )

#print('cwd',lib.nctest_ncgetcwd())

#FIXME: Much more, including globbing and whatever is likely to cause issues on
#Windows.
