#!/usr/bin/env bash
# Run a python script from this directory with the NCrystal development build
# (simplebuild), e.g.:
#
#   ./run.sh plot_filtertable.py stdlib::Al_sg225.ncmat
#
# The build must be done first, with "ncdevtool sb" (see devel/README.md).
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
exec "$HERE/../bin/ncdevtool" sbenv python "$@"
