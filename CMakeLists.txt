cmake_minimum_required(VERSION 2.8)
IF(POLICY CMP0020)
  CMAKE_POLICY(SET CMP0020 NEW)
ENDIF()
IF(POLICY CMP0028)
  CMAKE_POLICY(SET CMP0028 NEW)
ENDIF()
IF(POLICY CMP0037)
  CMAKE_POLICY(SET CMP0037 NEW)
ENDIF()
IF(POLICY CMP0047)
  CMAKE_POLICY(SET CMP0047 NEW)
ENDIF()
IF(POLICY CMP0042)
  CMAKE_POLICY(SET CMP0042 NEW)
ENDIF()

project(Gas-Kinetic-Scheme)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

if(MSVC)
   set(PLATFORM MSVC2015)
endif()

MESSAGE(STATUS "Build Project for ")
MESSAGE(STATUS ${PLATFORM})

#############################################################
###		SHARED OR STATIC LIBRARIES						  ###
#############################################################

include (${CMAKE_SOURCE_DIR}/cmake/general/LibUtilities.cmake)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries")
SharedLibs()


#############################################################
###		GENERAL INCLUDE DIRECTORIES						  ###
#############################################################

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

#############################################################
###		TEST SOURCES WITH/WITHOUT SUBFOLDERS			  ###
#############################################################

set(BUILD_TEST_SOURCE_GROUPS_WITH_SUBFOLDERS OFF CACHE BOOL "Test source groups with subfolders")

#############################################################
###		GENERAL FILE TYPES FOR LANGUAGE									  ###
#############################################################

include(${CMAKE_SOURCE_DIR}/cmake/general/LanguageFileTypes.cmake)
addCAndCPPFileTypes()

include(${CMAKE_SOURCE_DIR}/cmake/OpenMP/OpenMPCFlags.cmake)
addOpenMPCompilerFlags()

#############################################################
###		SUBDIRECTORIES									  ###
#############################################################

add_subdirectory(src/applications/libgks)
add_subdirectory(src/applications/main)
