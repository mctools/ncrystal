#!/bin/bash
set -e
set -u
set -x
SRCDIR="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )/../.." && pwd )"
test -f "${SRCDIR}/CMakeLists.txt"
test -f "${SRCDIR}/ncrystal_core/include/NCrystal/ncapi.h.in"

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
    -DNCRYSTAL_ENABLE_TESTING=ON \
    "$@"

#-DNCRYSTAL_ENABLE_THREADS=OFF

cmake --build "${TGT}/bld" --config "${THE_BUILD_TYPE}"
echo "Build dir was: ${TGT}/bld"
ctest
cmake --install "${TGT}/bld"
echo "Build dir was: ${TGT}/bld"

#Run examples:
"${TGT}/inst/bin/ncrystal_example_customphysics"
"${TGT}/inst/bin/ncrystal_example_cpp"
"${TGT}/inst/bin/ncrystal_example_c"

#Test downstream compilation:
BUILDFLAGS=$("${TGT}/inst/bin/ncrystal-config" --show buildflags)
c++ -std=c++11 ${SRCDIR}/examples/ncrystal_example_cpp.cc ${BUILDFLAGS} -o ${TGT}/app_cpp
${TGT}/app_cpp
c++ -std=c++11 ${SRCDIR}/examples/ncrystal_example_customphysics.cc ${BUILDFLAGS} -o ${TGT}/app_customphysics
${TGT}/app_customphysics
#Due to missing c++ stdlib symbols otherwise, we link with it here:
c++ -std=c++11 ${SRCDIR}/examples/ncrystal_example_c.c ${BUILDFLAGS} -o ${TGT}/app_c
${TGT}/app_c
