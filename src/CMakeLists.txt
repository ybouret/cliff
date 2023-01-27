CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

PROJECT(cliff)

SET(YROOT "$ENV{HOME}/program/upsylon/sdk" )

SET(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/..)

IF( "${YROOT}" STREQUAL "" )
        MESSAGE( FATAL_ERROR "Undefined YROOT" )
ENDIF( "${YROOT}" STREQUAL "" )

MESSAGE(STATUS "YROOT=${YROOT}")
INCLUDE("${YROOT}/share/y/y.cmake")
INCLUDE_DIRECTORIES("${YROOT}/include")
LINK_DIRECTORIES("${YROOT}/lib")

FUNCTION(CREATE SOMETHING)
        IF( ${ARGC} LESS 1 )
                MESSAGE( FATAL_ERROR "Missing Argument for CREATE" )
        ENDIF()
        MESSAGE( STATUS "CREATE args: '${ARGV0}' / '${ARGN}' for '${SOMETHING}'" )
ADD_EXECUTABLE(${ARGV0} ${ARGN})
TARGET_LINK_Y(${ARGV0} lua++ mkl)
INSTALL( TARGETS ${ARGV0} RUNTIME DESTINATION bin )
ENDFUNCTION(CREATE)


CREATE(dimer
main.cpp
types.hpp
hcell.hpp
psi.hpp
pcell.hpp
xcell.hpp
ko.hpp
nhe.hpp

hcell.lua
pcell.lua
xcell.lua
)
