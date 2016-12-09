include (GenerateExportHeader)
include(${CMAKE_SOURCE_DIR}/cmake/general/FileUtilities.cmake)

macro(sharedLibs)
	if(${BUILD_SHARED_LIBS})
		set(LIB_TYPE SHARED)
	else()
		set(LIB_TYPE STATIC)

		set(CompilerFlags
		        CMAKE_CXX_FLAGS
		        CMAKE_CXX_FLAGS_DEBUG
		        CMAKE_CXX_FLAGS_RELEASE
		        CMAKE_C_FLAGS
		        CMAKE_C_FLAGS_DEBUG
		        CMAKE_C_FLAGS_RELEASE
		        )
		foreach(CompilerFlag ${CompilerFlags})
		  string(REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
		endforeach()
	endif()
endmacro(sharedLibs)

macro(createLIB ${LIB_NAME} ${LIB_TYPE} ${MY_SRCS})
	add_library(${LIB_NAME} ${LIB_TYPE} ${MY_SRCS})
	if(${BUILD_SHARED_LIBS})
		GENERATE_EXPORT_HEADER	(	${LIB_NAME}
									BASE_NAME ${LIB_NAME}
									EXPORT_MACRO_NAME ${LIB_NAME}_EXPORT
									EXPORT_FILE_NAME ${CMAKE_SOURCE_DIR}/src/${LIB_NAME}/${LIB_NAME}_Export.h
									STATIC_DEFINE ${LIB_NAME}_BUILT_AS_STATIC
								)
	endif()
endmacro(createLIB)

#################################################################
###                     LIB OPTION                            ###
#################################################################
MACRO(setLibOptionOn)
	set(buildLib ON)
	set(buildTestSuite OFF)
	set (buildAllTests OFF)
ENDMACRO(setLibOptionOn)

MACRO(setTestSuiteOptionOn)
	set(buildLib OFF)
	set(buildTestSuite ON)
	set (buildAllTests OFF)
ENDMACRO(setTestSuiteOptionOn)

MACRO(setAllTestsOptionOn)
	set(buildLib OFF)
	set(buildTestSuite OFF)
	set (buildAllTests ON)
ENDMACRO(setAllTestsOptionOn)
