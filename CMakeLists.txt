cmake_minimum_required(VERSION 2.6)
project(Physics)
message(STATUS "-------------------")
message(STATUS "PROJECT: Physics --")
message(STATUS "-------------------")
get_filename_component(DEVROOT  ${CMAKE_CURRENT_SOURCE_DIR}/../.. ABSOLUTE)
list(APPEND CMAKE_MODULE_PATH  ${DEVROOT}/Util/cmakeUtil/  )
include(common_functions)
set(BUILD_SHARED_LIBS "SHARED")
define_archos_suffixe(ARCHOS)

file( GLOB sources  Physics/src/* Physics/materials/*)

find_package(Trellis)
find_package(Xfem)
find_package(Xext)
find_package(SolverInterfaces)

list(APPEND EXTERNAL_INCLUDES  
  ${TRELLIS_INCLUDE_DIR} 
  ${XFEM_INCLUDE_DIR} 
  ${XEXT_INCLUDE_DIR} 
  ${SolverBase_INCLUDE_DIR}  
  ${SuperLu_INCLUDE_DIR}  
  )
 
add_library(Physics ${BUILD_SHARED_LIBS} ${sources} )
include_directories(Physics ${EXTERNAL_INCLUDES} )
find_devroot(DEVROOT)
set_target_properties(Physics PROPERTIES COMPILE_FLAGS " -w " )

install(TARGETS Physics DESTINATION ${DEVROOT}/lib/${ARCHOS})



set_tags(Physics-tags "${sources}")

add_doxy("${CMAKE_CURRENT_SOURCE_DIR}/Physics/src")
add_doxy("${CMAKE_CURRENT_SOURCE_DIR}/Physics/materials")
