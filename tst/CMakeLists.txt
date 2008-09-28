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

IF(BLITZ_FOUND AND LINALG_FOUND)

	ADD_EXECUTABLE        (tcamclay01   tcamclay01.cpp)
	ADD_TEST              (tcamclay01   tcamclay01)
	TARGET_LINK_LIBRARIES (tcamclay01   ${LIBS})
	SET_TARGET_PROPERTIES (tcamclay01   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (ttruss01   ttruss01.cpp)
	ADD_TEST              (ttruss01   ttruss01)
	TARGET_LINK_LIBRARIES (ttruss01   ${LIBS})
	SET_TARGET_PROPERTIES (ttruss01   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (ttri6      ttri6.cpp)
	ADD_TEST              (ttri6      ttri6)
	TARGET_LINK_LIBRARIES (ttri6      ${LIBS})
	SET_TARGET_PROPERTIES (ttri6      PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tquad4     tquad4.cpp)
	ADD_TEST              (tquad4     tquad4)
	TARGET_LINK_LIBRARIES (tquad4     ${LIBS})
	SET_TARGET_PROPERTIES (tquad4     PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tquad8     tquad8.cpp)
	ADD_TEST              (tquad8     tquad8)
	TARGET_LINK_LIBRARIES (tquad8     ${LIBS})
	SET_TARGET_PROPERTIES (tquad8     PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tpstrain01 tpstrain01.cpp)
	ADD_TEST              (tpstrain01 tpstrain01)
	TARGET_LINK_LIBRARIES (tpstrain01 ${LIBS})
	SET_TARGET_PROPERTIES (tpstrain01 PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tpstrain02 tpstrain02.cpp)
	ADD_TEST              (tpstrain02 tpstrain02)
	TARGET_LINK_LIBRARIES (tpstrain02 ${LIBS})
	SET_TARGET_PROPERTIES (tpstrain02 PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tsolid01   tsolid01.cpp)
	ADD_TEST              (tsolid01   tsolid01)
	TARGET_LINK_LIBRARIES (tsolid01   ${LIBS})
	SET_TARGET_PROPERTIES (tsolid01   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tsolid02   tsolid02.cpp)
	ADD_TEST              (tsolid02   tsolid02)
	TARGET_LINK_LIBRARIES (tsolid02   ${LIBS})
	SET_TARGET_PROPERTIES (tsolid02   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tmesh01    tmesh01.cpp)
	ADD_TEST              (tmesh01    tmesh01)
	TARGET_LINK_LIBRARIES (tmesh01    ${LIBS})
	SET_TARGET_PROPERTIES (tmesh01    PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tmesh02    tmesh02.cpp)
	ADD_TEST              (tmesh02    tmesh02)
	TARGET_LINK_LIBRARIES (tmesh02    ${LIBS})
	SET_TARGET_PROPERTIES (tmesh02    PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tmesh03    tmesh03.cpp)
	ADD_TEST              (tmesh03    tmesh03)
	TARGET_LINK_LIBRARIES (tmesh03    ${LIBS})
	SET_TARGET_PROPERTIES (tmesh03    PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tpoisson01 tpoisson01.cpp)
	ADD_TEST              (tpoisson01 tpoisson01)
	TARGET_LINK_LIBRARIES (tpoisson01 ${LIBS})
	SET_TARGET_PROPERTIES (tpoisson01 PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tlaplace01 tlaplace01.cpp)
	ADD_TEST              (tlaplace01 tlaplace01)
	TARGET_LINK_LIBRARIES (tlaplace01 ${LIBS})
	SET_TARGET_PROPERTIES (tlaplace01 PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tsegfaults tsegfaults.cpp)
	ADD_TEST              (tsegfaults tsegfaults)
	TARGET_LINK_LIBRARIES (tsegfaults ${LIBS})
	SET_TARGET_PROPERTIES (tsegfaults PROPERTIES COMPILE_FLAGS "${FLAGS}")

ENDIF(BLITZ_FOUND AND LINALG_FOUND)

IF(IGRAPH_FOUND)
	ADD_EXECUTABLE        (ttree01 ttree01.cpp)
	ADD_TEST              (ttree01 ttree01)
	TARGET_LINK_LIBRARIES (ttree01 ${LIBS})
ENDIF(IGRAPH_FOUND)

SUBDIRS(reddy)
