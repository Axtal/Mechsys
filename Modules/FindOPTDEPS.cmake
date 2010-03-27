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

SET (FLAGS  "${FLAGS}")
SET (LIBS    ${LIBS})
SET (LFLAGS "${LFLAGS}")

FIND_PACKAGE (wxWidgets COMPONENTS base core aui)

INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindCGAL.cmake    )
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindFLTK.cmake    )
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindSUPERLU.cmake )
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindSUPERLUD.cmake)

SET (CGAL_OK FALSE)
if(CGAL_FOUND)
	if(GOTOBLAS_FOUND)
		SET (LIBS ${LIBS} ${CGAL_LIBRARIES})
        SET (CGAL_OK TRUE)
	else(GOTOBLAS_FOUND)
		ENABLE_LANGUAGE (Fortran)
		INCLUDE         (FindBLAS)
        if(BLAS_FOUND)
            SET (LIBS ${LIBS} ${CGAL_LIBRARIES} "gfortran")
            SET (CGAL_OK TRUE)
        endif(BLAS_FOUND)
	endif(GOTOBLAS_FOUND)
endif(CGAL_FOUND)
if(CGAL_OK)
	SET(FLAGS "${FLAGS} -DHAVE_CGAL")
endif(CGAL_OK)

if(FLTK_FOUND)
	SET(FLAGS "${FLAGS} ${FLTK_CFLAGS}")
	SET(LFLAGS "${LFLAGS} ${FLTK_LFLAGS}")
endif(FLTK_FOUND)

if(SUPERLU_FOUND)
	INCLUDE_DIRECTORIES (${SUPERLU_INCLUDE_DIRS})
	SET (FLAGS "${FLAGS} -DHAVE_SUPERLU")
	SET (LIBS   ${LIBS} ${SUPERLU_LIBRARIES})
endif(SUPERLU_FOUND)

if(SUPERLUD_FOUND)
	INCLUDE_DIRECTORIES (${SUPERLUD_INCLUDE_DIRS})
	SET (FLAGS "${FLAGS} -DHAVE_SUPERLUD")
	SET (LIBS   ${LIBS} ${SUPERLUD_LIBRARIES})
endif(SUPERLUD_FOUND)

IF(wxWidgets_FOUND)
    INCLUDE (${wxWidgets_USE_FILE})
    SET (LIBS ${LIBS} ${wxWidgets_LIBRARIES})
ENDIF(wxWidgets_FOUND)
