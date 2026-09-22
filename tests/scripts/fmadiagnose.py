#!/usr/bin/env python3

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2026 NCrystal developers                                   ##
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

# Exercises the ["util","fmadiagnose"] JSON query (similar to
# tests/src/app_fmadispatch but will be available in ncrystal-verify).

import platform

from NCrystalDev.misc import evaluate_query as ncquery


def classify_platform():
    #Independent (Python-side, not reusing any C++/CMake logic) guess at
    #whether this platform should have the technique enabled or disabled --
    #mirrors the static_assert cross-check in tests/src/app_fmadispatch. The
    #"disabled" side (architecture/OS only) is reliable from Python; the
    #"enabled" side can not check the actual compiler version used to build
    #the loaded NCrystal library, so it is a slightly weaker guess (x86 +
    #Linux/macOS only), but should hold in practice for any real wheel/conda
    #build or CI runner:
    mach = platform.machine().lower()
    system = platform.system()
    is_x86 = mach in ('x86_64', 'amd64', 'i386', 'i686', 'x86')
    is_arm = mach in ('arm64', 'aarch64') or mach.startswith('arm')
    is_windows = ( system == 'Windows' )
    is_linux_or_macos = system in ('Linux', 'Darwin')
    expect_enabled = is_x86 and is_linux_or_macos and not is_windows
    expect_disabled = is_arm or is_windows
    assert not ( expect_enabled and expect_disabled )
    return expect_enabled, expect_disabled


def main():
    expect_enabled, expect_disabled = classify_platform()

    res = ncquery(['util', 'fmadiagnose'])
    print("fmadiagnose query result:", res)

    enabled = res['enabled']
    correct = res['correct']
    optimised = res['optimised']
    speedup = res['speedup']

    assert isinstance(enabled, int) and enabled in (0, 1)
    assert correct is True, "dispatched and reference results differed!"

    if expect_enabled:
        assert enabled == 1, (
            "expected NCRYSTAL_FMADISPATCH to be enabled on this platform"
            " (x86, Linux/macOS), but it was not"
        )
    elif expect_disabled:
        assert enabled == 0, (
            "expected NCRYSTAL_FMADISPATCH to be disabled on this platform"
            " (ARM or Windows), but it was enabled"
        )
    else:
        print("Platform not confidently classified by this test;"
              " skipping the enabled/disabled expectation check.")

    if enabled and not optimised:
        print(f"Speedup: {speedup:.2f}x (unoptimised/debug build:"
              " not enforcing a minimum)")
    elif enabled:
        #Well below the ~5-15x typically measured (examples/
        #fmadispatch_investigate), to stay robust on a loaded/virtualised CI
        #machine while still failing hard on "no dispatch happened at all"
        #(which would show up as ~1x):
        assert speedup > 1.3, (
            f"NCRYSTAL_FMADISPATCH is enabled, but the measured speedup"
            f" ({speedup:.2f}x) is suspiciously low -- is it actually"
            f" being used?"
        )
        print(f"Speedup: {speedup:.2f}x (required: >1.3x)")
    else:
        print("NCRYSTAL_FMADISPATCH is not enabled on this build:"
              " skipping the speedup check.")


if __name__ == '__main__':
    main()
