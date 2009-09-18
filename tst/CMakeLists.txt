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
  tunnel01)

SET(PROGS_B # tests requiring CGAL
  alphashape01)

SET(PROGS_C # tests requiring GSL
  quadrature
  monte_carlo)

SET(PROGS_D # other tests 
 linalg01
 linalg02
 linsolver
 serialize)

SET(PROGS_F # mesh
  mesh01)

IF(BLITZ_FOUND AND LINALG_FOUND AND Boost_FOUND AND FLTK_FOUND)
	FOREACH(var ${PROGS_A})
		ADD_EXECUTABLE        (${var} "${var}.cpp")
		TARGET_LINK_LIBRARIES (${var} ${LIBS})
		SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
	ENDFOREACH(var)
ENDIF(BLITZ_FOUND AND LINALG_FOUND AND Boost_FOUND AND FLTK_FOUND)


IF(BLITZ_FOUND AND LINALG_FOUND AND CGAL_FOUND AND Boost_FOUND)
	FOREACH(var ${PROGS_B})
		ADD_EXECUTABLE        (${var} "${var}.cpp")
		TARGET_LINK_LIBRARIES (${var} ${LIBS})
		SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
	ENDFOREACH(var)
ENDIF(BLITZ_FOUND AND LINALG_FOUND AND CGAL_FOUND AND Boost_FOUND)


IF(GSL_FOUND)
	FOREACH(var ${PROGS_C})
		ADD_EXECUTABLE        (${var} "${var}.cpp")
		TARGET_LINK_LIBRARIES (${var} ${LIBS})
		SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
	ENDFOREACH(var)
ENDIF(GSL_FOUND)


IF(BLITZ_FOUND AND LINALG_FOUND)
	FOREACH(var ${PROGS_D})
		ADD_EXECUTABLE        (${var} "${var}.cpp")
		ADD_TEST              (${var} ${var})
		TARGET_LINK_LIBRARIES (${var} ${LIBS})
		SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
	ENDFOREACH(var)
ENDIF(BLITZ_FOUND AND LINALG_FOUND)


IF(TETGEN_FOUND AND TRIANGLE_FOUND)
	FOREACH(var ${PROGS_F})
		ADD_EXECUTABLE        (${var} "${var}.cpp")
		TARGET_LINK_LIBRARIES (${var} ${LIBS})
		SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
	ENDFOREACH(var)
ENDIF(TETGEN_FOUND AND TRIANGLE_FOUND)
