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

INCLUDE(FindThreads)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindUseBLITZ.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindGOTOBLAS.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindUMFPACK.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindSUPERLU.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindSUPERLUD.cmake)
INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindMPICXX.cmake)

ENABLE_LANGUAGE(Fortran)
INCLUDE(FindLAPACK)

SET(SERIAL_LIBS    "")
SET(SERIAL_FLAGS   "-DDO_DEBUG")
SET(PARALLEL_LIBS  "")
SET(PARALLEL_FLAGS "-DDO_DEBUG")

IF(UMFPACK_FOUND)
	INCLUDE_DIRECTORIES(${UMFPACK_INCLUDE_DIRS})
	SET(SERIAL_LIBS  ${SERIAL_LIBS} ${UMFPACK_LIBRARIES})
	SET(SERIAL_FLAGS "${SERIAL_FLAGS} -DHAVE_UMFPACK")
ENDIF(UMFPACK_FOUND)

IF(SUPERLU_FOUND)
	INCLUDE_DIRECTORIES(${SUPERLU_INCLUDE_DIRS})
	SET(SERIAL_LIBS  ${SERIAL_LIBS} ${SUPERLU_LIBRARIES})
	SET(SERIAL_FLAGS "${SERIAL_FLAGS} -DHAVE_SUPERLU")
ENDIF(SUPERLU_FOUND)

IF(SUPERLUD_FOUND)
	INCLUDE_DIRECTORIES(${SUPERLUD_INCLUDE_DIRS})
	SET(PARALLEL_LIBS  ${PARALLEL_LIBS} ${SUPERLUD_LIBRARIES})
	SET(PARALLEL_FLAGS "${PARALLEL_FLAGS} -DHAVE_SUPERLUD")
ENDIF(SUPERLUD_FOUND)

# SERIAL + LAPACK
IF(BLITZ_FOUND AND LAPACK_FOUND)
	ADD_LIBRARY           (msysfem  msysfem.cpp)
	TARGET_LINK_LIBRARIES (msysfem ${SERIAL_LIBS} ${LAPACK_LIBRARIES})
ENDIF(BLITZ_FOUND AND LAPACK_FOUND)
