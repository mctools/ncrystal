cmake_minimum_required(VERSION 3.20...3.30)

#Extract plugin name from ncplugin_name.txt:
file( STRINGS "${CMAKE_CURRENT_LIST_DIR}/ncplugin_name.txt" NCPlugin_NAME LIMIT_COUNT 1)
string(STRIP "${NCPlugin_NAME}" NCPlugin_NAME)

project( "NCPlugin_${NCPlugin_NAME}" VERSION 0.0.1 LANGUAGES CXX)

set( ncplugin_data_file_pattern "data/*.ncmat" )

if ( DEFINED SKBUILD_PROJECT_NAME )
  if ( NOT "${SKBUILD_PROJECT_NAME}" STREQUAL "ncrystal_plugin_${NCPlugin_NAME}" )
    message(
      FATAL_ERROR "Mismatch in plugin name hardcoded in "
      " ncplugin_name.txt and pyproject.toml's project.name field."
    )
  endif()
endif()

#Must always build against NCrystal:
if( NOT DEFINED "NCrystal_DIR" )
  #Need to invoke "ncrystal-config --show cmakedir" if we want to be able to
  #work with ncrystal-core installed via python wheels:
  execute_process(
    COMMAND ncrystal-config --show cmakedir
    OUTPUT_VARIABLE "NCrystal_DIR" OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()
find_package( NCrystal 4.0.0 REQUIRED )

function( ncrystal_srcfileglob varname pattern )
  #Glob while ignoring temporary files.
  file(
    GLOB tmpall LIST_DIRECTORIES false
    CONFIGURE_DEPENDS "${pattern}"
  )
  set(tmp "")
  foreach(fn ${tmpall})
    get_filename_component( bn "${fn}" NAME)
    if (bn MATCHES "(#|~| )+")
      message( WARNING "Ignoring file with invalid name: ${bn}")
    else()
      list(APPEND tmp "${fn}")
    endif()
  endforeach()
  set( ${varname} ${tmp} PARENT_SCOPE )
endfunction()

ncrystal_srcfileglob( plugin_srcfiles "${PROJECT_SOURCE_DIR}/src/*.cc" )
ncrystal_srcfileglob( plugin_hdrfiles  "${PROJECT_SOURCE_DIR}/src/*.hh" )
set( pluglib "NCPlugin_${NCPlugin_NAME}" )
add_library( ${pluglib} MODULE ${plugin_srcfiles} )
set_source_files_properties( ${plugin_srcfiles} PROPERTIES OBJECT_DEPENDS "${plugin_hdrfiles}" )
target_compile_definitions( ${pluglib} PRIVATE "NCPLUGIN_NAME=${NCPlugin_NAME}" "NCRYSTAL_NO_CMATH_CONSTANTS" )
target_link_libraries( ${pluglib} PRIVATE NCrystal::NCrystal )
target_include_directories( ${pluglib} PRIVATE "${PROJECT_SOURCE_DIR}/src" )

if ( ncplugin_data_file_pattern )
  file(GLOB plugin_datafiles LIST_DIRECTORIES false CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/${ncplugin_data_file_pattern}" )
else()
  set( plugin_datafiles "" )
endif()

if ( DEFINED SKBUILD_PROJECT_NAME )
  #Install in wheel platlib dir:
  set( pymoddir "${SKBUILD_PLATLIB_DIR}/ncrystal_plugin_${NCPlugin_NAME}")
  install( TARGETS ${pluglib} LIBRARY DESTINATION "${pymoddir}/plugins" )
  #Also add an empty __init__.py:
  file( TOUCH "${PROJECT_BINARY_DIR}/__init__.py" )
  INSTALL( FILES "${PROJECT_BINARY_DIR}/__init__.py" DESTINATION "${pymoddir}" )
  #Special location for data files:
  install( FILES ${plugin_datafiles} DESTINATION "${pymoddir}/data" )
else()
  install( TARGETS ${pluglib} LIBRARY DESTINATION lib )
  install( FILES ${plugin_datafiles} DESTINATION data )
endif()
