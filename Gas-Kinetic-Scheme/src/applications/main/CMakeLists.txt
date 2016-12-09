SET(EXE_NAME gks)

#clear variables
unset(MY_SRCS)

#################################################################
###		INCLUDE CMAKE MACROS    						      ###
#################################################################

include (${CMAKE_SOURCE_DIR}/cmake/general/LibUtilities.cmake)
include (${CMAKE_SOURCE_DIR}/cmake/general/FileUtilities.cmake)

#################################################################
###   	LINK DIRECTORIES		   			    	  	      ###
#################################################################


#################################################################
###   	EXCECUTABLE		   			    	  				  ###
#################################################################

OrganizeMainFiles()

SET(MY_SRCS ${MY_SRCS} ${SPECIFIC_FILES})
ADD_EXECUTABLE(${EXE_NAME} ${MY_SRCS})

#################################################################
###		LINK LIBRARY									      ###
#################################################################

target_link_libraries(${EXE_NAME} libgks)