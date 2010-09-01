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

# Flags
OPTION(A_MAKE_VERBOSE       "Show additional messages during compilation/linking?" OFF)
OPTION(A_MAKE_ALL_WARNINGS  "Make with all warnings (-Wall)"                       ON)
OPTION(A_MAKE_DEBUG_SYMBOLS "Make with debug symbols (-g)"                         ON)
OPTION(A_MAKE_PROFILING     "Make with profiling information (-pg)"                OFF)
OPTION(A_MAKE_OPTIMIZED     "Make optimized (-O3)"                                 OFF)
OPTION(A_MAKE_WXW_MONO      "Use wxWidgets monolithic libraries ?"                 ON)
OPTION(A_MAKE_TERM_WHITEBG  "Select colors for terminal with white background ?"   OFF)
OPTION(A_MAKE_TERM_NOCOLORS "Don't use colors when printing to terminal ?"         OFF)
OPTION(A_MAKE_STDVECTOR     "Use std::vector instead of own implemenatation ?"     OFF)
                                                                                   
# Options                                                                          
OPTION(A_USE_MPI            "Use OpenMPI ? "                                       OFF)
OPTION(A_USE_MTL4           "Use MTL4 instead of included Vector/Matrix library ?" ON)
OPTION(A_USE_WXW            "Use wxWidgets ?"                                      OFF)
OPTION(A_USE_VTK            "Use VTK ?"                                            OFF)
OPTION(A_USE_HDF5           "Use HDF5 ?"                                           OFF)
OPTION(A_USE_SUPERLU        "Use SueperLU"                                         OFF)
OPTION(A_USE_SUPERLUD       "Use SuperLUd"                                         OFF)
OPTION(A_USE_FLTK           "Use FLTK"                                             OFF)
OPTION(A_USE_CGAL           "Use CGAL"                                             OFF)
OPTION(A_USE_PARMETIS       "Use ParMETIS"                                         OFF)
OPTION(A_USE_MUMPS          "Use MUMPS"                                            OFF)

ADD_DEFINITIONS(-fmessage-length=0) # Each error message will appear on a single line; no line-wrapping will be done.

### FLAGS ###############################################################################################

IF(A_MAKE_VERBOSE)
	SET (CMAKE_VERBOSE_MAKEFILE TRUE)
ENDIF(A_MAKE_VERBOSE)

IF(A_MAKE_ALL_WARNINGS)
	ADD_DEFINITIONS (-Wall)
ENDIF(A_MAKE_ALL_WARNINGS)

IF(A_MAKE_DEBUG_SYMBOLS)
	ADD_DEFINITIONS (-g)
ENDIF(A_MAKE_DEBUG_SYMBOLS)

IF(A_MAKE_PROFILING)
	ADD_DEFINITIONS (-pg)
    SET (LFLAGS "${LFLAGS} -pg")
ENDIF(A_MAKE_PROFILING)

IF(A_MAKE_OPTIMIZED)
	ADD_DEFINITIONS (-O3)
ENDIF(A_MAKE_OPTIMIZED)

IF(A_MAKE_WXW_MONO)
    SET (WXW_COMPONENTS base core)
ELSE(A_MAKE_WXW_MONO)
    SET (WXW_COMPONENTS base core aui)
ENDIF(A_MAKE_WXW_MONO)

IF(A_MAKE_TERM_WHITEBG AND NOT A_MAKE_TERM_NOCOLORS)
    ADD_DEFINITIONS (-DTERM_WHITEBG)
ENDIF(A_MAKE_TERM_WHITEBG AND NOT A_MAKE_TERM_NOCOLORS)

IF(A_MAKE_TERM_NOCOLORS)
    ADD_DEFINITIONS (-DTERM_NOCOLORS)
ENDIF(A_MAKE_TERM_NOCOLORS)

IF(A_MAKE_STDVECTOR)
    ADD_DEFINITIONS (-DUSE_STDVECTOR)
ENDIF(A_MAKE_STDVECTOR)

### FIND DEPENDENCIES AND SET FLAGS AND LIBRARIES #######################################################

SET (FLAGS   "${FLAGS}")
SET (LIBS     ${LIBS})
SET (LFLAGS  "${LFLAGS}")
SET (MISSING "")

SET (Boost_USE_STATIC_LIBS ON)
ENABLE_LANGUAGE (Fortran)

FIND_PACKAGE (wxWidgets COMPONENTS ${WXW_COMPONENTS})       #  1
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
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindSUPERLU.cmake  ) # 17
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindSUPERLUD.cmake ) # 18
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindFLTK.cmake     ) # 19
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindCGAL.cmake     ) # 20
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindPROC.cmake     ) # 21
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindSPHASH.cmake   ) # 22

# 1
if(wxWidgets_FOUND AND A_USE_WXW)
    INCLUDE (${wxWidgets_USE_FILE})
    SET (LIBS ${LIBS} ${wxWidgets_LIBRARIES})
    ADD_DEFINITIONS (-DHAS_WXW)
else(wxWidgets_FOUND AND A_USE_WXW)
    if(A_USE_WXW)
        SET (MISSING "${MISSING} wxWidgets")
    endif(A_USE_WXW)
endif(wxWidgets_FOUND AND A_USE_WXW)

# 2
if(MPI_FOUND AND A_USE_MPI OR A_USE_MUMPS)
    INCLUDE_DIRECTORIES (${MPI_INCLUDE_PATH})
    SET (FLAGS  "${FLAGS}  ${MPI_COMPILE_FLAGS}")
    SET (LIBS    ${LIBS}   ${MPI_LIBRARIES}     )
    SET (LFLAGS "${LFLAGS} ${MPI_LINK_FLAGS}"   )
    ADD_DEFINITIONS (-DHAS_MPI)
else(MPI_FOUND AND A_USE_MPI OR A_USE_MUMPS)
    if(A_USE_MPI OR A_USE_MUMPS)
        SET (MISSING "${MISSING} OpenMPI")
    endif(A_USE_MPI OR A_USE_MUMPS)
endif(MPI_FOUND AND A_USE_MPI OR A_USE_MUMPS)

# 3
if(VTK_FOUND AND A_USE_VTK)
    INCLUDE_DIRECTORIES (${VTK_INCLUDE_DIRS})
    SET (LIBS  ${LIBS} vtkRendering)
    SET (FLAGS "${FLAGS} -DVTK_EXCLUDE_STRSTREAM_HEADERS")
else(VTK_FOUND AND A_USE_VTK)
    if(A_USE_VTK)
        SET (MISSING "${MISSING} VTK")
    endif(A_USE_VTK)
endif(VTK_FOUND AND A_USE_VTK)

# 4
if(HDF5_FOUND AND A_USE_HDF5)
    ADD_DEFINITIONS (-DH5_NO_DEPRECATED_SYMBOLS -DH5Gcreate_vers=2 -DH5Gopen_vers=2)
	INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
	SET (LIBS ${LIBS} ${HDF5_LIBRARIES})
else(HDF5_FOUND AND A_USE_HDF5)
    if(A_USE_HDF5)
        SET (MISSING "${MISSING} HDF5")
    endif(A_USE_HDF5)
endif(HDF5_FOUND AND A_USE_HDF5)

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
if(MTL_FOUND AND A_USE_MTL4)
	INCLUDE_DIRECTORIES (${MTL_INCLUDE_DIRS})
    ADD_DEFINITIONS (-DUSE_MTL4)
else(MTL_FOUND AND A_USE_MTL4)
    if(A_USE_MTL4)
        SET (MISSING "${MISSING} MTL4")
    endif(A_USE_MTL4)
endif(MTL_FOUND AND A_USE_MTL4)

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
if(MUMPS_FOUND AND A_USE_MUMPS)
	INCLUDE_DIRECTORIES (${MUMPS_INCLUDE_DIRS})
	SET(LIBS ${LIBS} ${MUMPS_LIBRARIES} -lmpi_f77)
    ADD_DEFINITIONS (-DHAS_MUMPS)
else(MUMPS_FOUND AND A_USE_MUMPS)
    if(A_USE_MUMPS)
        SET (MISSING "${MISSING} MUMPS")
    endif(A_USE_MUMPS)
endif(MUMPS_FOUND AND A_USE_MUMPS)

# 15
if(SCALAPACK_FOUND AND A_USE_MUMPS)
	SET(LIBS ${LIBS} ${SCALAPACK_LIBRARIES}) # must be after MUMPS
else(SCALAPACK_FOUND AND A_USE_MUMPS)
    if(A_USE_MUMPS)
        SET (MISSING "${MISSING} ScaLAPACK")
    endif(A_USE_MUMPS)
endif(SCALAPACK_FOUND AND A_USE_MUMPS)

# 16
if(PARMETIS_FOUND AND A_USE_PARMETIS OR A_USE_MUMPS)
    INCLUDE_DIRECTORIES (${PARMETIS_INCLUDE_DIRS})
    SET (LIBS ${LIBS} ${PARMETIS_LIBRARIES})
    ADD_DEFINITIONS (-DHAS_PARMETIS)
else(PARMETIS_FOUND AND A_USE_PARMETIS OR A_USE_MUMPS)
    if(A_USE_PARMETIS OR A_USE_MUMPS)
        SET (MISSING "${MISSING} ParMETIS")
    endif(A_USE_PARMETIS OR A_USE_MUMPS)
endif(PARMETIS_FOUND AND A_USE_PARMETIS OR A_USE_MUMPS)

# 17
if(SUPERLU_FOUND AND A_USE_SUPERLU)
	INCLUDE_DIRECTORIES (${SUPERLU_INCLUDE_DIRS})
	SET (LIBS ${LIBS} ${SUPERLU_LIBRARIES})
    ADD_DEFINITIONS (-DHAS_SUPERLU)
else(SUPERLU_FOUND AND A_USE_SUPERLU)
    if(A_USE_SUPERLU)
        SET (MISSING "${MISSING} SuperLU")
    endif(A_USE_SUPERLU)
endif(SUPERLU_FOUND AND A_USE_SUPERLU)

# 18
if(SUPERLUD_FOUND AND A_USE_SUPERLUD)
	INCLUDE_DIRECTORIES (${SUPERLUD_INCLUDE_DIRS})
	SET (LIBS ${LIBS} ${SUPERLUD_LIBRARIES})
    ADD_DEFINITIONS (-DHAS_SUPERLUD)
else(SUPERLUD_FOUND AND A_USE_SUPERLUD)
    if(A_USE_SUPERLUD)
        SET (MISSING "${MISSING} SuperLUd")
    endif(A_USE_SUPERLUD)
endif(SUPERLUD_FOUND AND A_USE_SUPERLUD)

# 19
if(FLTK_FOUND AND A_USE_FLTK)
	SET(FLAGS "${FLAGS} ${FLTK_CFLAGS}")
	SET(LFLAGS "${LFLAGS} ${FLTK_LFLAGS}")
    ADD_DEFINITIONS(-DHAS_FLTK)
else(FLTK_FOUND AND A_USE_FLTK)
    if(A_USE_FLTK)
        SET (MISSING "${MISSING} FLTK")
    endif(A_USE_FLTK)
endif(FLTK_FOUND AND A_USE_FLTK)

# 20
if(CGAL_FOUND AND A_USE_CGAL)
    SET (LIBS ${LIBS} ${CGAL_LIBRARIES} "gfortran")
	ADD_DEFINITIONS(-DHAS_CGAL)
else(CGAL_FOUND AND A_USE_CGAL)
    if(A_USE_CGAL)
        SET (MISSING "${MISSING} CGAL")
    endif(A_USE_CGAL)
endif(CGAL_FOUND AND A_USE_CGAL)

# 21
if(PROC_FOUND)
    INCLUDE_DIRECTORIES (${PROC_INCLUDE_DIRS})
    SET (LIBS ${LIBS} ${PROC_LIBRARIES})
	ADD_DEFINITIONS(-DHAS_PROC)
else(PROC_FOUND)
    SET (MISSING "${MISSING} proc")
endif(PROC_FOUND)

# 22
if(SPHASH_FOUND)
    INCLUDE_DIRECTORIES (${SPHASH_INCLUDE_DIRS})
    SET (LIBS ${LIBS} ${SPHASH_LIBRARIES})
	ADD_DEFINITIONS(-DHAS_SPHASH)
else(SPHASH_FOUND)
    SET (MISSING "${MISSING} sparsehash")
endif(SPHASH_FOUND)
