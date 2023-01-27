########################################################################
##
##
##
## This scripts configure a GENERATOR with COMPILERS
##
##
##
########################################################################
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS    ON)

MESSAGE( "" )
MESSAGE( STATUS "Detecting configuration parameters" )
########################################################################
##
##
## Dectecting BUILD_TYPE
##
##
########################################################################
MESSAGE( STATUS "Detecting BUILD_TYPE...'${BUILD_TYPE}'" )
#MESSAGE( FATAL_ERROR "stop" )

SET(KNOW_BUILD_TYPE OFF)

#default
IF("" STREQUAL "${BUILD_TYPE}")
	SET(BUILD_TYPE "Release")
ENDIF()

IF("Debug" STREQUAL "${BUILD_TYPE}")
	SET(KNOWN_BUILD_TYPE ON)
ENDIF()

IF("Release" STREQUAL "${BUILD_TYPE}")
	SET(KNOWN_BUILD_TYPE ON)
ENDIF()

IF(NOT KNOWN_BUILD_TYPE)
	MESSAGE( FATAL_ERROR "Unknown BUILD_TYPE=${BUILD_TYPE}" )
ENDIF()

########################################################################
##
##
## Dectecting GENERATOR
##
##
########################################################################
MESSAGE( STATUS "Detecting GENERATOR..." )
SET(KNOWN_GENERATOR OFF) #check if we know it
SET(MULTI_GENERATOR OFF) #generator handles mutliple builds (Xcode, Visual Studio...): default is OFF
SET(TRACE_COMPILERS ON)  #if trace of compiler version may be kept (gcc,clang...)    : default is ON

#default
IF( "" STREQUAL "${GENERATOR}" )
    FIND_PROGRAM(Y_NINJA ninja)
    IF( "Y_NINJA-NOTFOUND" STREQUAL "${Y_NINJA}")
		SET(GENERATOR "Unix Makefiles")
	ELSE()
		SET(GENERATOR "Ninja")
	ENDIF()
ENDIF()

########################################################################
##
## Unix Makefiles: single output with compilers version
##
########################################################################
IF( "${GENERATOR}" STREQUAL "Unix Makefiles" )
	SET(KNOWN_GENERATOR ON)
ENDIF()

########################################################################
##
## Ninja: single output with compilers version
##
########################################################################
IF( "${GENERATOR}" STREQUAL "Ninja" )
	SET(KNOWN_GENERATOR ON)
ENDIF()

########################################################################
##
## Visual Studio: mutliple outputs, single compilers
##
########################################################################
IF( "${GENERATOR}" MATCHES "Visual Studio.*" )
	SET(KNOWN_GENERATOR ON)
	SET(MULTI_GENERATOR ON)
	SET(COMPILERS "microsoft") #cl.exe in MSVS
	SET(VERSION "")
	SET(TRACE_COMPILERS OFF)
ENDIF()

########################################################################
##
## Codeblocks: single output with compilers version
##
########################################################################
IF( "${GENERATOR}" MATCHES "CodeBlocks.*" )
	SET(KNOWN_GENERATOR ON)
ENDIF()

########################################################################
##
## Kate: single output with compilers version
##
########################################################################
IF( "${GENERATOR}" MATCHES "Kate - .*" )
	SET(KNOWN_GENERATOR ON)
ENDIF()

########################################################################
##
## Xcode: mutliple outputs, single compilers
##
########################################################################
IF( "${GENERATOR}" STREQUAL "Xcode" )
	SET(KNOWN_GENERATOR  ON)
	SET(MULTI_GENERATOR  ON)
	SET(COMPILERS   "clang") #AppleClang in Xcode
	SET(VERSION     "")
	SET(TRACE_COMPILERS OFF)
ENDIF()

########################################################################
##
## CodeLite: single output with compilers version
##
########################################################################
IF( "${GENERATOR}" MATCHES "CodeLite - .*" )
	SET(KNOWN_GENERATOR ON)
ENDIF()

IF(NOT KNOWN_GENERATOR)
	MESSAGE( FATAL_ERROR "Unknown GENERATOR=${GENERATOR}" )
ENDIF()


########################################################################
##
##
## Dectecting COMPILERS
##
##
########################################################################
MESSAGE( STATUS "Detecting COMPILERS..." )

SET(KNOWN_COMPILERS OFF)

#default if not set
IF( "" STREQUAL "${COMPILERS}")
	SET(COMPILERS "gnu")
	FIND_PROGRAM(UNAME uname)
	IF( "${UNAME}" STREQUAL "UNAME-NOTFOUND" )
		MESSAGE( STATUS "no uname => ${COMPILERS}" )
	ELSE()
 		EXECUTE_PROCESS(COMMAND         "${UNAME}" -s
						OUTPUT_VARIABLE CMAKE_SYSTEM_NAME
						OUTPUT_STRIP_TRAILING_WHITESPACE)
		MESSAGE( STATUS "CMAKE_SYSTEM_NAME is <${CMAKE_SYSTEM_NAME}>" )
		STRING(COMPARE EQUAL "Darwin"  ${CMAKE_SYSTEM_NAME} Y_DARWIN)
		STRING(COMPARE EQUAL "FreeBSD" ${CMAKE_SYSTEM_NAME} Y_FREEBSD)
 		STRING(COMPARE EQUAL "OpenBSD" ${CMAKE_SYSTEM_NAME} Y_OPENBSD)
 		IF(Y_DARWIN OR Y_FREEBSD OR Y_OPENBSD)
			SET(COMPILERS "clang")
 		ENDIF()
 	ENDIF()
ENDIF()

IF( "${COMPILERS}" STREQUAL "gnu" )
	SET(CC  "gcc${VERSION}")
	SET(CXX "g++${VERSION}")
	SET(KNOWN_COMPILERS ON)
ENDIF()

IF( "${COMPILERS}" STREQUAL "clang" )
	SET(CC  "clang${VERSION}")
	SET(CXX "clang++${VERSION}")
	SET(KNOWN_COMPILERS ON)
ENDIF()

IF( "${COMPILERS}" STREQUAL "intel" )
	SET(CC  "icc${VERSION}")
	SET(CXX "icpc${VERSION}")
	SET(KNOWN_COMPILERS ON)
	SET(VERSION "")
ENDIF()

IF( "${COMPILERS}" STREQUAL "microsoft" )
	SET(CC  "cl.exe")
	SET(CXX "cl.exe")
	SET(KNOWN_COMPILERS ON)
	SET(VERSION "")
ENDIF()

IF(NOT KNOWN_COMPILERS)
	MESSAGE( FATAL_ERROR "Unknown COMPILERS=${COMPILERS}")
ENDIF()

########################################################################
##
##
## Dectecting SOURCE
##
##
########################################################################
MESSAGE( STATUS "Detecting SOURCE...")
IF("" STREQUAL "${SOURCE}")
	SET(SOURCE "src")
ENDIF()


########################################################################
##
##
## prepare the build path
##
##
########################################################################
IF(TRACE_COMPILERS)
	SET(BUILD_ROOT "forge/${COMPILERS}${VERSION}")
ELSE()
	SET(BUILD_ROOT "forge")
ENDIF()

STRING(REGEX REPLACE " " "" __GENPATH ${GENERATOR})
IF(MULTI_GENERATOR)
	SET(BUILD_PATH "${BUILD_ROOT}/${__GENPATH}")
ELSE()
	SET(BUILD_PATH "${BUILD_ROOT}/${__GENPATH}/${BUILD_TYPE}")
ENDIF()

MESSAGE( "" )
MESSAGE( STATUS "Configuration summary:")
MESSAGE( STATUS "\t(*) BUILD_TYPE = [${BUILD_TYPE}]" )
MESSAGE( STATUS "\t(*) GENERATOR  = [${GENERATOR}]"  )
MESSAGE( STATUS "\t(*) COMPILERS  = [${COMPILERS}]"  )
MESSAGE( STATUS "\t(*) VERSION    = '${VERSION}'"    ) 
MESSAGE( STATUS "\t(*) SOURCE     = '${SOURCE}'"     )
MESSAGE( STATUS "\t(*) BUILD_PATH = '${BUILD_PATH}'" )

########################################################################
##
##
## and call cmake
##
##
########################################################################
MESSAGE( "" )
MESSAGE( STATUS "Creating build directory...")
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_PATH})
MESSAGE( STATUS "Calling CMake...")
IF(MULTI_GENERATOR)
	EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND}
							-G${GENERATOR}
							-DCMAKE_C_COMPILER=${CC}
							-DCMAKE_CXX_COMPILER=${CXX}
							${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE}
							WORKING_DIRECTORY ${BUILD_PATH})
ELSE()
	EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND}
							-G${GENERATOR}
							-DCMAKE_C_COMPILER=${CC}
							-DCMAKE_CXX_COMPILER=${CXX}
							-DCMAKE_BUILD_TYPE=${BUILD_TYPE}
							${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE}
							WORKING_DIRECTORY ${BUILD_PATH})
ENDIF()
MESSAGE("")
MESSAGE( STATUS "Ready to be compiled for [${COMPILERS}${VERSION}] with [${GENERATOR}]")
