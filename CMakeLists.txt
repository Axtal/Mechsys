########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################

PROJECT (MECHSYS)

CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
ENABLE_TESTING         ()

OPTION(MAKE_VERBOSE       "Show additional messages during compilation/link?" OFF)
OPTION(MAKE_ALL_WARNINGS  "Make with all warnings (-Wall)"                    ON )
OPTION(MAKE_DEBUG_SYMBOLS "Make with debug symbols (-g3)"                     ON )
OPTION(MAKE_OPTIMIZED     "Make optimized (-O3)"                              OFF)

IF(MAKE_VERBOSE)
	SET(CMAKE_VERBOSE_MAKEFILE TRUE)
ENDIF(MAKE_VERBOSE)

IF(MAKE_ALL_WARNINGS)
	ADD_DEFINITIONS(-Wall)
ENDIF(MAKE_ALL_WARNINGS)

IF(MAKE_DEBUG_SYMBOLS)
	ADD_DEFINITIONS(-g3)
ENDIF(MAKE_DEBUG_SYMBOLS)

IF(MAKE_OPTIMIZED)
	ADD_DEFINITIONS(-O3)
ENDIF(MAKE_OPTIMIZED)

ADD_DEFINITIONS(-fmessage-length=0) # Each error message will appear on a single line; no line-wrapping will be done.

INCLUDE_DIRECTORIES(${MECHSYS_SOURCE_DIR}/lib/)

SUBDIRS (lib tst src)

# Dependencies
SET(Boost_USE_STATIC_LIBS ON)
INCLUDE(FindPythonLibs)
INCLUDE(FindBoost)
INCLUDE(FindThreads)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindUseBLITZ.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindGOTOBLAS.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindUMFPACK.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindSUPERLU.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindSUPERLUD.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindMPICXX.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindIGRAPH.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindCGAL.cmake)

SET(FLAGS "-DDO_DEBUG")
SET(LIBS  "jrs_triangle")

IF(PYTHONLIBS_FOUND)
	INCLUDE_DIRECTORIES   (${PYTHON_INCLUDE_PATH})
	EXECUTE_PROCESS (COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
	SET(LIBS ${LIBS} ${PYTHON_LIBRARIES})
ENDIF(PYTHONLIBS_FOUND)

IF(Boost_FOUND)
	INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
	SET(FLAGS "${FLAGS} -DUSE_BOOST")
ENDIF(Boost_FOUND)

SET(LINALG_FOUND FALSE)
IF(GOTOBLAS_FOUND)
	SET(LIBS  ${LIBS}  ${GOTOBLAS_LIBRARIES})
	SET(LINALG_FOUND TRUE)
ELSE(GOTOBLAS_FOUND)
	ENABLE_LANGUAGE(Fortran)
	INCLUDE(FindLAPACK)
	IF(LAPACK_FOUND)
		SET(LIBS  ${LIBS}  ${LAPACK_LIBRARIES})
		SET(LINALG_FOUND TRUE)
	ENDIF(LAPACK_FOUND)
ENDIF(GOTOBLAS_FOUND)

IF(UMFPACK_FOUND)
	INCLUDE_DIRECTORIES(${UMFPACK_INCLUDE_DIRS})
	SET(FLAGS "${FLAGS} -DHAVE_UMFPACK")
	SET(LIBS   ${LIBS} ${UMFPACK_LIBRARIES})
ENDIF(UMFPACK_FOUND)

IF(SUPERLU_FOUND)
	INCLUDE_DIRECTORIES(${SUPERLU_INCLUDE_DIRS})
	SET(FLAGS "${FLAGS} -DHAVE_SUPERLU")
	SET(LIBS   ${LIBS} ${SUPERLU_LIBRARIES})
ENDIF(SUPERLU_FOUND)

IF(SUPERLUD_FOUND)
	INCLUDE_DIRECTORIES(${SUPERLUD_INCLUDE_DIRS})
	SET(FLAGS "${FLAGS} -DHAVE_SUPERLUD")
	SET(LIBS   ${LIBS} ${SUPERLUD_LIBRARIES})
ENDIF(SUPERLUD_FOUND)

IF(IGRAPH_FOUND)
	SET(LIBS  ${LIBS}  ${IGRAPH_LIBRARIES})
ENDIF(IGRAPH_FOUND)

IF(CGAL_FOUND)
	SET(LIBS  ${LIBS}  ${CGAL_LIBRARIES})
ENDIF(CGAL_FOUND)
