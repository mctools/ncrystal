#!/bin/bash
set -e
set -u
SRCDIR="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )/../.." && pwd )"
test -f "${SRCDIR}/CMakeLists.txt"
test -f "${SRCDIR}/ncrystal_core/include/NCrystal/ncapi.h"

TGT="/tmp/${USER}/ncrystal_bldlocal"
rm -rf "${TGT}/bld" "${TGT}/inst"
mkdir -p "${TGT}/bld" "${TGT}/inst"
cd "${TGT}/bld"

THE_BUILD_TYPE=Release
cmake \
    -S "${SRCDIR}" \
    -B "${TGT}/bld" \
    -DCMAKE_BUILD_TYPE="${THE_BUILD_TYPE}" \
    -DCMAKE_INSTALL_PREFIX="${TGT}/inst" \
    -DNCRYSTAL_ENABLE_EXAMPLES=ON \
    -DNCRYSTAL_ENABLE_GEANT4=OFF \
    -DNCRYSTAL_BUILD_STRICT=ON \
    -DNCRYSTAL_ENABLE_SETUPSH=ON \
    "$@"

#-DNCRYSTAL_ENABLE_THREADS=OFF

cmake --build "${TGT}/bld" --config "${THE_BUILD_TYPE}"
echo "Build dir was: ${TGT}/bld"
ctest
cmake --install "${TGT}/bld"
echo "Build dir was: ${TGT}/bld"
