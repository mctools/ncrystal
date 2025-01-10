cmake_minimum_required(VERSION 3.28...3.30)

#Extract plugin name from ncplugin_name.txt:
file( STRINGS "${CMAKE_CURRENT_LIST_DIR}/ncplugin_name.txt" NCPlugin_NAME LIMIT_COUNT 1)
string(STRIP "${NCPlugin_NAME}" NCPlugin_NAME)

project( "NCPlugin_${NCPlugin_NAME}" VERSION 0.0.1 LANGUAGES CXX)

#Must always build against NCrystal:
if( NOT DEFINED "NCrystal_DIR" )
  execute_process(
    COMMAND ncrystal-config --show cmakedir
    OUTPUT_VARIABLE "NCrystal_DIR" OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()
find_package( NCrystal 3.9.80 REQUIRED )

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

function (
    ncrystal_add_plugin_library
    plugin_name
    source_directory
  )
  ncrystal_srcfileglob( plugin_srcfiles "${source_directory}/*.cc" )
  ncrystal_srcfileglob( plugin_hdrfiles  "${source_directory}/*.hh" )
  set( pluglib "NCPlugin_${plugin_name}" )
  add_library( ${pluglib} MODULE ${plugin_srcfiles} )
  set_source_files_properties( ${plugin_srcfiles} PROPERTIES OBJECT_DEPENDS "${plugin_hdrfiles}" )
  target_compile_definitions( ${pluglib} PRIVATE "NCPLUGIN_NAME=${plugin_name}" "NCRYSTAL_NO_CMATH_CONSTANTS" )
  target_link_libraries( ${pluglib} PRIVATE NCrystal::NCrystal )
  target_include_directories( ${pluglib} PRIVATE "${source_directory}" )
  install( TARGETS ${pluglib} LIBRARY DESTINATION lib )
endfunction()

ncrystal_add_plugin_library( "${NCPlugin_NAME}" "${PROJECT_SOURCE_DIR}/src" )
