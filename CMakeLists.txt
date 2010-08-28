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
CMAKE_MINIMUM_REQUIRED (VERSION 2.8)
ENABLE_TESTING         ()

# Features                                                                         
OPTION(A_WITH_PYTHON_MODULE "Create python module (mechsys)" OFF)
OPTION(A_WITH_TEST          "Compile general tests"          OFF)
OPTION(A_WITH_TFEM          "Compile FEM tests"              OFF)
OPTION(A_WITH_TLBM          "Compile LBM tests"              OFF)
OPTION(A_WITH_TDEM          "Compile DEM tests"              OFF)
OPTION(A_WITH_TSPH          "Compile SPH tests"              OFF)

# needed when building Debian package
IF(NOT INSTALL_ROOT)
	SET(INSTALL_ROOT /)
ENDIF(NOT INSTALL_ROOT)

INCLUDE_DIRECTORIES (${MECHSYS_SOURCE_DIR})

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

IF(ALLDIRS)
    SUBDIRS(${ALLDIRS})
ENDIF(ALLDIRS)

### DEPENDENCIES ########################################################################################

INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindDEPS.cmake)
if(MISSING)
    MESSAGE("Missing dependencies =${MISSING}")
endif(MISSING)

### PYTHON MODULE #######################################################################################

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
