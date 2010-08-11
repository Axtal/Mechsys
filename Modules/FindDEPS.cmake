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

SET (FLAGS   "${FLAGS}")
SET (LIBS     ${LIBS})
SET (LFLAGS  "${LFLAGS}")
SET (MISSING "")

SET (Boost_USE_STATIC_LIBS ON)
ENABLE_LANGUAGE (Fortran)

if(A_MAKE_WX_MONO)
    FIND_PACKAGE (wxWidgets COMPONENTS base core)
else(A_MAKE_WX_MONO)
    FIND_PACKAGE (wxWidgets COMPONENTS base core aui)       #  1
endif(A_MAKE_WX_MONO)

INCLUDE (FindMPI)                                           #  2
INCLUDE (FindVTK)                                           #  3
INCLUDE (FindHDF5)                                          #  4
INCLUDE (FindBoost)                                         #  5
INCLUDE (FindLAPACK)                                        #  6
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindBLITZ.cmake    ) #  7
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindGSL.cmake      ) #  8
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindMTL.cmake      ) #  9
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindVORO.cmake     ) # 10
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindTETGEN.cmake   ) # 11
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindTRIANGLE.cmake ) # 12
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindUMFPACK.cmake  ) # 13
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindMUMPS.cmake    ) # 14
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindSCALAPACK.cmake) # 15
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindPARMETIS.cmake ) # 16
#INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindSUPERLU.cmake  ) # 17
#INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindSUPERLUD.cmake ) # 18
#INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindFLTK.cmake     ) # 19
#INCLUDE(${MECHSYS_SOURCE_DIR}/Modules/FindCGAL.cmake     ) # 20

# 1
if(wxWidgets_FOUND)
    INCLUDE (${wxWidgets_USE_FILE})
    SET (LIBS ${LIBS} ${wxWidgets_LIBRARIES})
    ADD_DEFINITIONS (-DHAS_WXW)
else(wxWidgets_FOUND)
    if(A_USE_WXW)
        SET (MISSING "${MISSING} wxWidgets")
    endif(A_USE_WXW)
endif(wxWidgets_FOUND)

# 2
if(MPI_FOUND)
    INCLUDE_DIRECTORIES (${MPI_INCLUDE_PATH})
    SET (FLAGS  "${FLAGS}  ${MPI_COMPILE_FLAGS}")
    SET (LIBS    ${LIBS}   ${MPI_LIBRARIES}     )
    SET (LFLAGS "${LFLAGS} ${MPI_LINK_FLAGS}"   )
else(MPI_FOUND)
    SET (MISSING "${MISSING} OpenMPI")
endif(MPI_FOUND)

# 3
if(VTK_FOUND)
    INCLUDE_DIRECTORIES (${VTK_INCLUDE_DIRS})
    SET (LIBS  ${LIBS} vtkRendering)
    SET (FLAGS "${FLAGS} -DVTK_EXCLUDE_STRSTREAM_HEADERS")
else(VTK_FOUND)
    SET (MISSING "${MISSING} VTK")
endif(VTK_FOUND)

# 4
if(HDF5_FOUND)
    ADD_DEFINITIONS (-DH5_NO_DEPRECATED_SYMBOLS -DH5Gcreate_vers=2 -DH5Gopen_vers=2)
	INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
	SET (LIBS ${LIBS} ${HDF5_LIBRARIES})
else(HDF5_FOUND)
    SET (MISSING "${MISSING} HDF5")
endif(HDF5_FOUND)

# 5
if(Boost_FOUND)
	INCLUDE_DIRECTORIES (${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
    SET (MISSING "${MISSING} boost")
endif(Boost_FOUND)

# 6
if(LAPACK_FOUND)
    SET (LIBS ${LIBS} ${LAPACK_LIBRARIES})
else(LAPACK_FOUND)
    INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindLocLAPACK.cmake)
    if(LocLAPACK_FOUND)
        SET (LIBS ${LIBS} ${LocLAPACK_LIBRARIES} "gfortran")
    else(LocLAPACK_FOUND)
        SET (MISSING "${MISSING} LaPACK")
    endif(LocLAPACK_FOUND)
endif(LAPACK_FOUND)

# 7
if(BLITZ_FOUND)
	INCLUDE_DIRECTORIES (${BLITZ_INCLUDE_DIRS})
else(BLITZ_FOUND)
    SET (MISSING "${MISSING} blitz++")
endif(BLITZ_FOUND)

# 8
if(GSL_FOUND)
	INCLUDE_DIRECTORIES (${GSL_INCLUDE_DIRS})
	SET (LIBS ${LIBS} ${GSL_LIBRARIES})
else(GSL_FOUND)
    SET (MISSING "${MISSING} GSL")
endif(GSL_FOUND)

# 9
if(MTL_FOUND)
	INCLUDE_DIRECTORIES (${MTL_INCLUDE_DIRS})
else(MTL_FOUND)
    SET (MISSING "${MISSING} MTL4")
endif(MTL_FOUND)

# 10
if(VORO_FOUND)
	INCLUDE_DIRECTORIES (${VORO_INCLUDE_DIRS})
else(VORO_FOUND)
    SET (MISSING "${MISSING} Voro++")
endif(VORO_FOUND)

# 11
if(TETGEN_FOUND)
	INCLUDE_DIRECTORIES (${TETGEN_INCLUDE_DIRS})
	SET (LIBS ${LIBS} ${TETGEN_LIBRARIES})
else(TETGEN_FOUND)
    SET (MISSING "${MISSING} Tetgen")
endif(TETGEN_FOUND)

# 12
if(TRIANGLE_FOUND)
	INCLUDE_DIRECTORIES (${TRIANGLE_INCLUDE_DIRS})
	SET (LIBS ${LIBS} ${TRIANGLE_LIBRARIES})
else(TRIANGLE_FOUND)
    SET (MISSING "${MISSING} Triangle")
endif(TRIANGLE_FOUND)

# 13
if(UMFPACK_FOUND)
	INCLUDE_DIRECTORIES (${UMFPACK_INCLUDE_DIRS})
	SET(LIBS ${LIBS} ${UMFPACK_LIBRARIES})
    ADD_DEFINITIONS (-DHAS_UMFPACK)
else(UMFPACK_FOUND)
    SET (MISSING "${MISSING} UMFPACK")
endif(UMFPACK_FOUND)

# 14
if(MUMPS_FOUND)
	INCLUDE_DIRECTORIES (${MUMPS_INCLUDE_DIRS})
	SET(LIBS ${LIBS} ${MUMPS_LIBRARIES} -lmpi_f77)
    ADD_DEFINITIONS (-DHAS_MUMPS)
else(MUMPS_FOUND)
    SET (MISSING "${MISSING} MUMPS")
endif(MUMPS_FOUND)

# 15
if(SCALAPACK_FOUND)
	SET(LIBS ${LIBS} ${SCALAPACK_LIBRARIES}) # must be after MUMPS
else(SCALAPACK_FOUND)
    SET (MISSING "${MISSING} ScaLAPACK")
endif(SCALAPACK_FOUND)

# 16
if(PARMETIS_FOUND)
    INCLUDE_DIRECTORIES (${PARMETIS_INCLUDE_DIRS})
    SET (LIBS ${LIBS} ${PARMETIS_LIBRARIES})
    ADD_DEFINITIONS (-DHAS_PARMETIS)
else(PARMETIS_FOUND)
    SET (MISSING "${MISSING} ParMETIS")
endif(PARMETIS_FOUND)

# 17
if(SUPERLU_FOUND)
	INCLUDE_DIRECTORIES (${SUPERLU_INCLUDE_DIRS})
	SET (LIBS ${LIBS} ${SUPERLU_LIBRARIES})
    ADD_DEFINITIONS (-DHAS_SUPERLU)
else(SUPERLU_FOUND)
    #SET (MISSING "${MISSING} SuperLU")
endif(SUPERLU_FOUND)

# 18
if(SUPERLUD_FOUND)
	INCLUDE_DIRECTORIES (${SUPERLUD_INCLUDE_DIRS})
	SET (LIBS ${LIBS} ${SUPERLUD_LIBRARIES})
    ADD_DEFINITIONS (-DHAS_SUPERLUD)
else(SUPERLUD_FOUND)
    #SET (MISSING "${MISSING} SuperLUd")
endif(SUPERLUD_FOUND)

# 19
if(FLTK_FOUND)
	SET(FLAGS "${FLAGS} ${FLTK_CFLAGS}")
	SET(LFLAGS "${LFLAGS} ${FLTK_LFLAGS}")
    ADD_DEFINITIONS(-DHAS_FLTK)
else(FLTK_FOUND)
    #SET (MISSING "${MISSING} FLTK")
endif(FLTK_FOUND)

# 20
if(CGAL_FOUND)
    SET (LIBS ${LIBS} ${CGAL_LIBRARIES} "gfortran")
	ADD_DEFINITIONS(-DHAS_CGAL)
else(CGAL_FOUND)
    #SET (MISSING "${MISSING} CGAL")
endif(CGAL_FOUND)
