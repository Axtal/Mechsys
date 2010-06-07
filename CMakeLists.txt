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

PROJECT                (MECHSYS)
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
ENABLE_TESTING         ()

# Flags
OPTION(A_MAKE_VERBOSE       "Show additional messages during compilation/link?" OFF)
OPTION(A_MAKE_ALL_WARNINGS  "Make with all warnings (-Wall)"                    ON)
OPTION(A_MAKE_DEBUG_SYMBOLS "Make with debug symbols (-g)"                      ON)
OPTION(A_MAKE_PROFILING     "Make with profiling information (-pg)"             OFF)
OPTION(A_MAKE_OPTIMIZED     "Make optimized (-O3)"                              OFF)
OPTION(A_MAKE_WX_MONO       "Use wxWidgets monolithic libraries ?"              ON)

# Features
OPTION(A_WITH_PYTHON_MODULE "Create python module (mechsys)"                    OFF)
OPTION(A_WITH_TEST          "Compile tests"                                     OFF)
OPTION(A_WITH_TFEM          "Compile FEM tests"                                 OFF)
OPTION(A_WITH_TLBM          "Compile LBM tests"                                 OFF)
OPTION(A_WITH_TDEM          "Compile DEM tests"                                 ON)
OPTION(A_WITH_TSPH          "Compile SPH tests"                                 ON)
OPTION(A_WITH_OPTIONAL      "Make with optional dependencies"                   OFF)

# Options
OPTION(A_USE_MPI            "Compile parallel code with OpenMPI ?"              OFF)
OPTION(A_OPTIMIZE_DEM_MEM   "Optimize DEM code for memory instead of speed"     OFF)

# needed when building Debian package
IF(NOT INSTALL_ROOT)
	SET(INSTALL_ROOT /)
ENDIF(NOT INSTALL_ROOT)

ADD_DEFINITIONS     (-fmessage-length=0) # Each error message will appear on a single line; no line-wrapping will be done.
INCLUDE_DIRECTORIES (${MECHSYS_SOURCE_DIR})

### FLAGS ###############################################################################################

IF(A_MAKE_VERBOSE)
	SET(CMAKE_VERBOSE_MAKEFILE TRUE)
ENDIF(A_MAKE_VERBOSE)

IF(A_MAKE_ALL_WARNINGS)
	ADD_DEFINITIONS(-Wall)
ENDIF(A_MAKE_ALL_WARNINGS)

IF(A_MAKE_DEBUG_SYMBOLS)
	ADD_DEFINITIONS(-g)
ENDIF(A_MAKE_DEBUG_SYMBOLS)

IF(A_MAKE_PROFILING)
	SET(CMAKE_CXX_FLAGS -pg)
ENDIF(A_MAKE_PROFILING)

IF(A_MAKE_OPTIMIZED)
	ADD_DEFINITIONS(-O3)
ENDIF(A_MAKE_OPTIMIZED)

IF(A_MAKE_WX_MONO)
	SET(WXMONO TRUE)
ENDIF(A_MAKE_WX_MONO)

if(A_OPTIMIZE_DEM_MEM)
    ADD_DEFINITIONS (-DOPTIMIZE_DEM_MEM)
endif(A_OPTIMIZE_DEM_MEM)

### FEATURES ############################################################################################

SET(ALLDIRS) # sub-directories

IF(A_WITH_PYTHON_MODULE)
	SET(ALLDIRS ${ALLDIRS} src)
ENDIF(A_WITH_PYTHON_MODULE)

IF(A_WITH_TEST)
	SET(ALLDIRS ${ALLDIRS} test)
ENDIF(A_WITH_TEST)

IF(A_WITH_TFEM)
	SET(ALLDIRS ${ALLDIRS} tfem)
ENDIF(A_WITH_TFEM)

IF(A_WITH_TLBM)
	SET(ALLDIRS ${ALLDIRS} tlbm)
ENDIF(A_WITH_TLBM)

IF(A_WITH_TDEM)
	SET(ALLDIRS ${ALLDIRS} tdem)
ENDIF(A_WITH_TDEM)

IF(A_WITH_TSPH)
	SET(ALLDIRS ${ALLDIRS} tsph)
ENDIF(A_WITH_TSPH)

SUBDIRS(${ALLDIRS})

### DEPENDENCIES ########################################################################################

# essential dependencies
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindDEPS.cmake)
if(NOT DEPS_OK)
    MESSAGE("Not all dependencies were satisfied")
endif(NOT DEPS_OK)

# optionall dependencies
if(A_WITH_OPTIONAL)
    INCLUDE (FindMPI)
    INCLUDE (FindThreads)
    INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindOPTDEPS.cmake)
endif(A_WITH_OPTIONAL)

if(A_WITH_PYTHON_MODULE)
    INCLUDE (FindPythonLibs)
    if(PYTHONLIBS_FOUND)
        INCLUDE_DIRECTORIES (${PYTHON_INCLUDE_PATH})
        EXECUTE_PROCESS (COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
        FIND_LIBRARY (BOOST_PYTHON NAMES boost_python-mt PATHS /usr/lib)
        IF(NOT BOOST_PYTHON)
            FIND_LIBRARY (BOOST_PYTHON NAMES boost_python PATHS /usr/lib)
        ENDIF(NOT BOOST_PYTHON)
        SET(LIBS ${LIBS} ${PYTHON_LIBRARIES} ${BOOST_PYTHON})
    else(PYTHONLIBS_FOUND)
        MESSAGE("To compile the Python module, boost-python must be installed")
    endif(PYTHONLIBS_FOUND)
endif(A_WITH_PYTHON_MODULE)
