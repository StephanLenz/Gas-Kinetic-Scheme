macro(addFileEndingToCollect file_ending)
#input: file_ending --> appends it to the list of files to collect
	list (FIND files_to_collect ${file_ending} index)
	if (${index} EQUAL -1)
		set(files_to_collect ${files_to_collect} ${file_ending})
	endif()
#output: files_to_collect
endmacro(addFileEndingToCollect)

MACRO(OrganizeMainFiles)
#input: -
	unset(SPECIFIC_FILES)

	file(GLOB ALL_FILES ./*.*)
	set(SPECIFIC_FILES ${ALL_FILES})

	foreach(_FILE ${SPECIFIC_FILES})
		if(${_FILE} MATCHES "main.cpp")
			source_group(main FILES ${_FILE})
		else()
			source_group(general FILES ${_FILE})
		endif()
	endforeach()
#output: ${SPECIFIC_FILES}
ENDMACRO(OrganizeMainFiles)

MACRO(OrganizeSourceFiles)
#input: Lib, Path, global files_to_collect
	setPathAndLibName()

	collectAllFilesFrom(${path} "${files_to_collect}")

	SetSourceGroup(${libName})

#output:
ENDMACRO(OrganizeSourceFiles)

MACRO(setPathAndLibName)
#input: -

	set(libName ${LIB_NAME})
 	set(path ${CMAKE_CURRENT_LIST_DIR})

#output: libName, path
ENDMACRO(setPathAndLibName)

MACRO(collectAllFilesFrom path file_endings)
#input: path from files to collect
unset(allFilesInPath)

	foreach(_ending ${file_endings})
		FILE(GLOB Files_With_Ending ${path}/${_ending})
		set(allFilesInPath ${allFilesInPath} ${Files_With_Ending})
	endforeach()

	addFilesToAllSources("${allFilesInPath}")

#output: -
ENDMACRO(collectAllFilesFrom)

MACRO(addFilesToAllSources filesToAdd)
#input: MY_SRCS, TEST_SRCS

	separateFiles("${filesToAdd}")

	SET(MY_SRCS ${MY_SRCS} ${PACKAGE_SRCS})
	SET(TEST_SRCS ${TEST_SRCS} ${PACKAGE_TEST_SRCS})

#output: MY_SRCS, TEST_SRCS
ENDMACRO(addFilesToAllSources)

MACRO(separateFiles filesToSeperate)
#input: allFilesInPath (all source and test files)
unset(PACKAGE_SRCS)
unset(PACKAGE_TEST_SRCS)

foreach (_FILE ${filesToSeperate})
	if(${_FILE} MATCHES "Test" OR ${_FILE} MATCHES "Mock")
		SET(PACKAGE_TEST_SRCS ${PACKAGE_TEST_SRCS} ${_FILE})
	else()
		SET(PACKAGE_SRCS ${PACKAGE_SRCS} ${_FILE})
	endif()
endforeach()

#output: PACKAGE_SRCS, PACKAGE_TEST_SRCS
ENDMACRO(separateFiles)

MACRO(SetSourceGroup lib_name)
#input: LIB_NAME PACKAGE_SRCS
unset(folderAfterName)
setFolderAfterName(${lib_name})

if(buildLib)
	SetSourceGroupForSourceFiles(\${folderAfterName} \${PACKAGE_SRCS})
else()
	SetSourceGroupForTestFiles(${lib_name} \${folderAfterName} \${PACKAGE_TEST_SRCS})
endif()

#output: -
ENDMACRO(SetSourceGroup)

MACRO(setFolderAfterName targetName)
#input: targetName (e.g. lib name, exe name)
unset(folderAfterName)

	string(REPLACE "/" ";" currentDirList ${CMAKE_CURRENT_LIST_DIR})

	SET(findTargetName 0)

	FOREACH(FOLDER ${currentDirList})
		if(findTargetName)
			SET(folderAfterName ${folderAfterName}\\${FOLDER})
		endif()

		if(${FOLDER} STREQUAL ${targetName})
			SET(findTargetName 1)
		endif()
	ENDFOREACH()

#output: folderAfterName (sourcegroups)
ENDMACRO(setFolderAfterName)

MACRO(SetSourceGroupForSourceFiles sourceGroup sourceFiles)
#input: folderAfterName, PACKAGE_SRCS
	SOURCE_GROUP(${sourceGroup} FILES ${sourceFiles})
#output: -
ENDMACRO(SetSourceGroupForSourceFiles)

MACRO(SetSourceGroupForTestFiles libName path testFiles)
#input: libName path testFiles

if(${BUILD_TEST_SOURCE_GROUPS_WITH_SUBFOLDERS})
	TestSourceGroups_WithSubFolders()
else()
	TestSourceGroups_WithoutSubFolders(${libName} \${path} \${testFiles})
endif()

#output: -
ENDMACRO(SetSourceGroupForTestFiles)

MACRO(TestSourceGroups_WithoutSubFolders libName path testFiles)
#input: libName path testFiles

if(buildAllTests)
	set(_SOURCE_GROUP ${libName}\\${path})
elseif(buildTestSuite)
	set(_SOURCE_GROUP ${path})
endif()

source_group(${_SOURCE_GROUP} FILES ${testFiles})

#output: -
ENDMACRO(TestSourceGroups_WithoutSubFolders)

MACRO(TestSourceGroups_WithSubFolders)
string(LENGTH ${CMAKE_CURRENT_LIST_DIR} PATH_LENGTH)
MATH(EXPR PATH_LENGTH ${PATH_LENGTH}+1) #because of the last slash

foreach(testFile ${PACKAGE_TEST_SRCS})
	string(LENGTH ${testFile} FILE_NAME_LENGTH)

	if(${testFile} MATCHES Test.h)
		MATH(EXPR ENDING_LENGTH 6)
	elseif(${testFile} MATCHES Test.cpp)
		MATH(EXPR ENDING_LENGTH 8)
	elseif(${testFile} MATCHES Mocks.h)
		MATH(EXPR ENDING_LENGTH 7)
	elseif(${testFile} MATCHES Mocks.cpp)
		MATH(EXPR ENDING_LENGTH 9)
	elseif(${testFile} MATCHES .cpp)
		MATH(EXPR ENDING_LENGTH 4)
	elseif(${testFile} MATCHES .h)
		MATH(EXPR ENDING_LENGTH 2)
	endif()

	MATH(EXPR FILE_NAME_LENGTH ${FILE_NAME_LENGTH}-${PATH_LENGTH}-${ENDING_LENGTH})
	string(SUBSTRING ${testFile} ${PATH_LENGTH} ${FILE_NAME_LENGTH} FOLDER_NAME)

	if(buildTestSuite)
		set(ROOT ${folderAfterName})
	elseif(buildAllTests)
		set(ROOT ${LIB_NAME}\\${folderAfterName})
	endif()

	source_group(${ROOT}\\${FOLDER_NAME} FILES ${testFile})

endforeach()

ENDMACRO(TestSourceGroups_WithSubFolders)
