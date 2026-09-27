#!/usr/bin/env bash
set -euo pipefail

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

# Windows-only developer tool, deliberately NOT wired into simplebuild
# (simplebuild does not support Windows, and this repo's day-to-day
# development does not happen there): records a Windows Performance
# Recorder (WPR) CPU-sampling trace of tests/src/app_perfvdos (the
# VDOS Gn expansion/FastConvolve hot path benchmark) and (best-effort)
# summarises it via the devel/win/winprofile .NET analyzer tool.
#
# Assumes:
#  - Running under git-bash on Windows (as on GitHub's windows-* runners,
#    or a real Windows dev machine with Git for Windows installed).
#  - wpr.exe on PATH (ships with Windows itself, not a separate install).
#  - dotnet SDK on PATH (needed to build/run the analyzer tool).
#  - <build-dir> already configured+built (CMake, target "perfvdos") with
#    debug info available (e.g. -DCMAKE_BUILD_TYPE=RelWithDebInfo -- plain
#    Release does not produce PDBs by default, and without PDBs the
#    analyzer cannot resolve our own code's function names, only raw
#    addresses).
#
# Usage: profile_perfvdos.sh <build-dir> <output-dir>
# Writes to <output-dir>: perfvdos.log, trace.etl, analysis.txt (if the
# analyzer succeeds), and copies of the relevant .pdb files (for local
# WPA -- Windows Performance Analyzer, a GUI tool -- analysis, regardless
# of whether the automated analysis above worked).

if [ $# -ne 2 ]; then
  echo "usage: $0 <build-dir> <output-dir>" >&2
  exit 1
fi
build_dir="$1"
out_dir="$2"
mkdir -p "$out_dir"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for cand in "$build_dir/tests/Release/perfvdos.exe" \
            "$build_dir/tests/RelWithDebInfo/perfvdos.exe" \
            "$build_dir/tests/Debug/perfvdos.exe"; do
  if [ -f "$cand" ]; then
    EXE="$cand"
    break
  fi
done
if [ -z "${EXE:-}" ]; then
  echo "Could not locate the built perfvdos executable under $build_dir" >&2
  exit 1
fi
echo "Using executable: $EXE"

#Windows has no rpath equivalent: running this directly (rather than via
#ctest, which sets this up itself -- see tests/CMakeLists.txt) needs the
#directory holding the (namespace-protected, hence "NCrystal-test.dll" not
#"NCrystal.dll") shared library prepended to PATH ourselves, or the exe
#fails to even start (confirmed to matter in practice, twice, for this
#exact style of direct-invocation Windows workflow):
dll="$(find "$build_dir" -iname 'NCrystal-test.dll' -o -iname 'NCrystal.dll' 2>/dev/null | head -n1)"
if [ -n "$dll" ]; then
  echo "Prepending to PATH: $(dirname "$dll")"
  export PATH="$(dirname "$dll"):$PATH"
else
  echo "WARNING: could not locate the NCrystal shared library under $build_dir" >&2
fi

trace="$out_dir/trace.etl"
rm -f "$trace"

echo "Starting WPR CPU trace..."
wpr -start CPU -filemode

#Always try to stop WPR even if the run itself fails or this script is
#interrupted, so we never leave a dangling system-wide trace session
#behind on the (ephemeral, but still) runner:
cleanup() { wpr -cancel >/dev/null 2>&1 || true; }
trap cleanup EXIT

set +e
"$EXE" > "$out_dir/perfvdos.log" 2>&1
ec=$?
set -e
echo "perfvdos exit code: $ec" | tee -a "$out_dir/perfvdos.log"

echo "Stopping WPR trace -> $trace"
wpr -stop "$trace"
trap - EXIT

#Copy PDBs alongside the trace, so a full WPA-based analysis remains
#possible locally even if the analyzer below has issues:
find "$build_dir" -iname '*.pdb' -exec cp {} "$out_dir/" \; 2>/dev/null || true

#Best-effort automated analysis (the caller decides whether a failure
#here should fail the overall job -- this script itself always exits 0
#past this point, since the trace+PDBs are already safely in out_dir):
echo "Running analyzer..."
if command -v dotnet >/dev/null 2>&1; then
  if dotnet run --project "$script_dir/winprofile" -c Release -- \
       "$trace" "perfvdos" "$out_dir" > "$out_dir/analysis.txt" 2>&1; then
    echo "Analyzer succeeded, see analysis.txt"
  else
    echo "Analyzer failed (see analysis.txt for details) -- the raw" \
         "trace.etl and .pdb files are still available for local WPA analysis."
  fi
else
  echo "dotnet not found on PATH, skipping automated analysis" \
       "(trace.etl and .pdb files are still available for local WPA analysis)." \
       > "$out_dir/analysis.txt"
fi

exit "$ec"
