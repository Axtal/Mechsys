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
OPTION(MAKE_ALL_WARNINGS  "Make with all warnings (-Wall)"                    OFF)
OPTION(MAKE_DEBUG_SYMBOLS "Make with debug symbols (-g)"                      OFF)
OPTION(MAKE_PROFILING     "Make with profiling information (-pg)"             OFF)
OPTION(MAKE_OPTIMIZED     "Make optimized (-O3)"                              ON)
OPTION(MAKE_PYTHON_MODULE "Create python module (mechsys)"                    OFF)
OPTION(MAKE_TST           "Compile tests"                                     OFF)
OPTION(MAKE_TFEM          "Compile FEM tests"                                 OFF)
OPTION(MAKE_TLBM          "Compile LBM tests"                                 OFF)
OPTION(MAKE_TDEM          "Compile DEM tests"                                 ON)

# needed when building Debian package
IF(NOT INSTALL_ROOT)
	SET(INSTALL_ROOT /)
ENDIF(NOT INSTALL_ROOT)

#SET(CMAKE_CXX_COMPILER /opt/intel/cc/10.1.018/bin/icpc)

IF(MAKE_VERBOSE)
	SET(CMAKE_VERBOSE_MAKEFILE TRUE)
ENDIF(MAKE_VERBOSE)

IF(MAKE_ALL_WARNINGS)
	ADD_DEFINITIONS(-Wall)
ENDIF(MAKE_ALL_WARNINGS)

IF(MAKE_DEBUG_SYMBOLS)
	ADD_DEFINITIONS(-g)
ENDIF(MAKE_DEBUG_SYMBOLS)

IF(MAKE_PROFILING)
	SET(CMAKE_CXX_FLAGS -pg)
ENDIF(MAKE_PROFILING)

IF(MAKE_OPTIMIZED)
	ADD_DEFINITIONS(-O3)
ENDIF(MAKE_OPTIMIZED)

IF(MAKE_STATIC)
	ADD_DEFINITIONS(-static)
ENDIF(MAKE_STATIC)

ADD_DEFINITIONS(-fPIC -fmessage-length=0) # Each error message will appear on a single line; no line-wrapping will be done.

INCLUDE_DIRECTORIES(${MECHSYS_SOURCE_DIR}/lib/)

# sub-directories
SET(ALLDIRS)

IF(MAKE_PYTHON_MODULE)
	SET(ALLDIRS ${ALLDIRS} src)
ENDIF(MAKE_PYTHON_MODULE)

IF(MAKE_TST)
	SET(ALLDIRS ${ALLDIRS} tst)
ENDIF(MAKE_TST)

IF(MAKE_TFEM)
	SET(ALLDIRS ${ALLDIRS} tfem)
ENDIF(MAKE_TFEM)

IF(MAKE_TLBM)
	SET(ALLDIRS ${ALLDIRS} tlbm)
ENDIF(MAKE_TLBM)

IF(MAKE_TDEM)
	SET(ALLDIRS ${ALLDIRS} tdem)
ENDIF(MAKE_TDEM)

SUBDIRS(${ALLDIRS})

########################################################################### Dependencies
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
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindCGAL.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindFLTK.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindGSL.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindMTL.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindVORO.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindTRIANGLE.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindTETGEN.cmake)

SET(FLAGS  "")
SET(LIBS   "")
SET(LFLAGS "")

IF(PYTHONLIBS_FOUND AND MAKE_PYTHON_MODULE)
	INCLUDE_DIRECTORIES (${PYTHON_INCLUDE_PATH})
	EXECUTE_PROCESS (COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
	FIND_LIBRARY (BOOST_PYTHON NAMES boost_python-mt PATHS /usr/lib)
	IF(NOT BOOST_PYTHON)
		FIND_LIBRARY (BOOST_PYTHON NAMES boost_python PATHS /usr/lib)
	ENDIF(NOT BOOST_PYTHON)
	SET(LIBS ${LIBS} ${PYTHON_LIBRARIES} ${BOOST_PYTHON})
ENDIF(PYTHONLIBS_FOUND AND MAKE_PYTHON_MODULE)

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
	ELSE(LAPACK_FOUND)
		INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindLAPACK.cmake)
		SET(LIBS ${LIBS} ${LAPACK_LIBRARIES} "gfortran")
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

IF(CGAL_FOUND)
	IF(GOTOBLAS_FOUND)
		SET(LIBS  ${LIBS}  ${CGAL_LIBRARIES})
	ELSE(GOTOBLAS_FOUND)
		ENABLE_LANGUAGE(Fortran)
		INCLUDE(FindBLAS)
		SET(LIBS  ${LIBS}  ${CGAL_LIBRARIES})
		SET(LIBS  ${LIBS}  "gfortran")
	ENDIF(GOTOBLAS_FOUND)
	SET(FLAGS "${FLAGS} -DUSE_CGAL")
ENDIF(CGAL_FOUND)

IF(FLTK_FOUND)
	SET(FLAGS "${FLAGS} ${FLTK_CFLAGS}")
	SET(LFLAGS "${LFLAGS} ${FLTK_LFLAGS}")
ENDIF(FLTK_FOUND)

IF(GSL_FOUND)
	SET(LIBS ${LIBS} ${GSL_LIBRARIES})
ENDIF(GSL_FOUND)

IF(VORO_FOUND)
	INCLUDE_DIRECTORIES(${VORO_INCLUDE_DIRS})
	SET(LIBS  ${LIBS}  ${VORO_LIBRARIES})
ENDIF(VORO_FOUND)

IF(TRIANGLE_FOUND)
	INCLUDE_DIRECTORIES(${TRIANGLE_INCLUDE_DIRS})
	SET(LIBS  ${LIBS}  ${TRIANGLE_LIBRARIES})
ENDIF(TRIANGLE_FOUND)

IF(TETGEN_FOUND)
	INCLUDE_DIRECTORIES(${TETGEN_INCLUDE_DIRS})
	SET(LIBS  ${LIBS}  ${TETGEN_LIBRARIES})
ENDIF(TETGEN_FOUND)

IF(MTL_FOUND)
	INCLUDE_DIRECTORIES(${MTL_INCLUDE_DIRS})
ENDIF(MTL_FOUND)

IF(GSL_FOUND)
	INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
ENDIF(GSL_FOUND)

IF(BLITZ_FOUND)
	INCLUDE_DIRECTORIES(${BLITZ_INCLUDE_DIRS})
ENDIF(BLITZ_FOUND)

FIND_PACKAGE(wxWidgets COMPONENTS base core aui)
IF(wxWidgets_FOUND)
    INCLUDE(${wxWidgets_USE_FILE})
    SET(LIBS ${LIBS} ${wxWidgets_LIBRARIES})
ENDIF(wxWidgets_FOUND)

########################################################################### FEM

SET(FEM_FOUND 0)
IF(BLITZ_FOUND AND LINALG_FOUND AND Boost_FOUND AND MTL_FOUND AND TRIANGLE_FOUND AND TETGEN_FOUND)
	SET(FEM_FOUND 1)
ENDIF(BLITZ_FOUND AND LINALG_FOUND AND Boost_FOUND AND MTL_FOUND AND TRIANGLE_FOUND AND TETGEN_FOUND)
