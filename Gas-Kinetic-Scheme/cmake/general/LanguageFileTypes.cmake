include(${CMAKE_SOURCE_DIR}/cmake/general/FileUtilities.cmake)

macro(addCAndCPPFileTypes)
	addFileEndingToCollect("*.h")
	addFileEndingToCollect("*.c")
	addFileEndingToCollect("*.cpp")
endmacro(addCAndCPPFileTypes)

macro(addObjCAndObjCPPFileTypesToCollect)
	addFileEndingToCollect("*.m")
	addFileEndingToCollect("*.mm")
	addFileEndingToCollect("*.h")
endmacro(addObjCAndObjCPPFileTypesToCollect)
