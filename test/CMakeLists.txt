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

SET(PROGS_A # tests requiring FLTK
  tplotxy
  )

SET(PROGS_B # tests requiring CGAL
  alphashape01)

SET(PROGS_C # tests requiring GSL
  quadrature
  monte_carlo
  test_numdiff1
  test_dunitvec
  test_ode1
  test_timesteps)

SET(PROGS_D # other tests 
 linalg01
 linalg02
 linalg03
 linalg04
 tlaexpr
 linsolver
 test_maps
 test_readfiles
 test_array
 test_colors
 serialize
 test_paragrid)

SET(PROGS_F # mesh
  mesh01
  mesh02
  mesh03
  mesh04
  delaunay)

SET(PROGS_G # tests requiring wxWidgets
  wxplotxy
  )


IF(FLTK_FOUND AND A_USE_FLTK)
	FOREACH(var ${PROGS_A})
		ADD_EXECUTABLE        (${var} "${var}.cpp")
		TARGET_LINK_LIBRARIES (${var} ${LIBS})
		SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
	ENDFOREACH(var)
ENDIF(FLTK_FOUND AND A_USE_FLTK)


IF(CGAL_FOUND AND A_USE_CGAL)
	FOREACH(var ${PROGS_B})
		ADD_EXECUTABLE        (${var} "${var}.cpp")
		TARGET_LINK_LIBRARIES (${var} ${LIBS})
		SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
	ENDFOREACH(var)
ENDIF(CGAL_FOUND AND A_USE_CGAL)


IF(GSL_FOUND)
	FOREACH(var ${PROGS_C})
		ADD_EXECUTABLE        (${var} "${var}.cpp")
		TARGET_LINK_LIBRARIES (${var} ${LIBS})
		SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
	ENDFOREACH(var)
ENDIF(GSL_FOUND)


FOREACH(var ${PROGS_D})
    ADD_EXECUTABLE        (${var} "${var}.cpp")
    ADD_TEST              (${var} ${var})
    TARGET_LINK_LIBRARIES (${var} ${LIBS})
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
ENDFOREACH(var)


IF(TETGEN_FOUND AND TRIANGLE_FOUND)
	FOREACH(var ${PROGS_F})
		ADD_EXECUTABLE        (${var} "${var}.cpp")
		TARGET_LINK_LIBRARIES (${var} ${LIBS})
		SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
	ENDFOREACH(var)
ENDIF(TETGEN_FOUND AND TRIANGLE_FOUND)


IF(wxWidgets_FOUND AND A_USE_WXW)
	FOREACH(var ${PROGS_G})
		ADD_EXECUTABLE        (${var} "${var}.cpp")
		TARGET_LINK_LIBRARIES (${var} ${LIBS})
		SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
	ENDFOREACH(var)
ENDIF(wxWidgets_FOUND AND A_USE_WXW)


if (A_USE_VTK)
    SUBDIRS(graph models vtk)
else (A_USE_VTK)
    SUBDIRS(graph models)
endif (A_USE_VTK)
