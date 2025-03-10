#!/usr/bin/env bash
set -e
set -u
set -x
SRCDIR="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )/../.." && pwd )"
test -f "${SRCDIR}/CMakeLists.txt"
test -f "${SRCDIR}/ncrystal_core/include/NCrystal/ncapi.h.in"

TGT="/tmp/${USER}/ncrystal_bldlocal_plugin"
rm -rf "${TGT}/bld" "${TGT}/inst" "${TGT}/bld_plugin" "${TGT}/inst_plugin"
mkdir -p "${TGT}/bld"
cd "${TGT}"

THE_BUILD_TYPE=Release
cmake \
    -S "${SRCDIR}" \
    -B "${TGT}/bld" \
    -DCMAKE_BUILD_TYPE="${THE_BUILD_TYPE}" \
    -DCMAKE_INSTALL_PREFIX="${TGT}/inst" \
    -DNCRYSTAL_ENABLE_EXAMPLES=ON \
    -DNCRYSTAL_BUILD_STRICT=ON \
    -DNCRYSTAL_ENABLE_TESTING=OFF \
    "$@"

#-DNCRYSTAL_ENABLE_THREADS=OFF

cmake --build "${TGT}/bld" --config "${THE_BUILD_TYPE}"
echo "Build dir was: ${TGT}/bld"
ctest --build-config  "${THE_BUILD_TYPE}"  --output-on-failure --test-output-size-failed 10000 --test-output-truncation middle
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

#Test plugin:
mkdir -p "${TGT}/bld_plugin"
#cd "${TGT}/bld_plugin"
cmake \
    -S "${SRCDIR}/examples/plugin" \
    -B "${TGT}/bld_plugin" \
    -DCMAKE_BUILD_TYPE="${THE_BUILD_TYPE}" \
    -DCMAKE_INSTALL_PREFIX="${TGT}/inst_plugin" \
    -DNCrystal_DIR="$("${TGT}/inst/bin/ncrystal-config" --show cmakedir)"

cmake --build "${TGT}/bld_plugin" --config "${THE_BUILD_TYPE}"
echo "Plugin build dir was: ${TGT}/bld_plugin"
ctest --build-config  "${THE_BUILD_TYPE}"  --output-on-failure --test-output-size-failed 10000 --test-output-truncation middle
cmake --install "${TGT}/bld_plugin"
echo "Plugin build dir was: ${TGT}/bld_plugin"
echo "Plugin was installed in: ${TGT}/inst_plugin"

#export LD_LIBRARY_PATH="$("${TGT}/inst/bin/ncrystal-config" --show libdir):${LD_LIBRARY_PATH}"
export NCRYSTAL_PLUGIN_LIST=$(echo "${TGT}/inst_plugin/lib/libNCPlugin_DummyPlugin".*)
export NCRYSTAL_REQUIRED_PLUGINS="DummyPlugin"
export NCRYSTAL_PLUGIN_RUNTESTS=1
#Now run anything which triggers plugin-loading, to have the plugin-availability
#test carried out:
${TGT}/app_c
echo "All ok"
