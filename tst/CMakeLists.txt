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
	IF(FLTK_FOUND)
		ADD_EXECUTABLE        (ttunnel   ttunnel.cpp)
		TARGET_LINK_LIBRARIES (ttunnel   ${LIBS})
		SET_TARGET_PROPERTIES (ttunnel   PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
	ENDIF(FLTK_FOUND)

	ADD_EXECUTABLE        (tserialize   tserialize.cpp)
	ADD_TEST              (tserialize   tserialize)
	TARGET_LINK_LIBRARIES (tserialize   ${LIBS} boost_python boost_serialization)
	SET_TARGET_PROPERTIES (tserialize   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tmesh01   tmesh01.cpp)
	ADD_TEST              (tmesh01   tmesh01)
	TARGET_LINK_LIBRARIES (tmesh01   ${LIBS} boost_python)
	SET_TARGET_PROPERTIES (tmesh01   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (ttriax   ttriax.cpp)
	ADD_TEST              (ttriax   ttriax)
	TARGET_LINK_LIBRARIES (ttriax   ${LIBS} boost_python)
	SET_TARGET_PROPERTIES (ttriax   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tbiot01   tbiot01.cpp)
	ADD_TEST              (tbiot01   tbiot01)
	TARGET_LINK_LIBRARIES (tbiot01   ${LIBS})
	SET_TARGET_PROPERTIES (tbiot01   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tterz01   tterz01.cpp)
	ADD_TEST              (tterz01   tterz01)
	TARGET_LINK_LIBRARIES (tterz01   ${LIBS})
	SET_TARGET_PROPERTIES (tterz01   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tembank01   tembank01.cpp)
	ADD_TEST              (tembank01   tembank01)
	TARGET_LINK_LIBRARIES (tembank01   ${LIBS})
	SET_TARGET_PROPERTIES (tembank01   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tembank02   tembank02.cpp)
	ADD_TEST              (tembank02   tembank02)
	TARGET_LINK_LIBRARIES (tembank02   ${LIBS})
	SET_TARGET_PROPERTIES (tembank02   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (texcav01   texcav01.cpp)
	ADD_TEST              (texcav01   texcav01)
	TARGET_LINK_LIBRARIES (texcav01   ${LIBS})
	SET_TARGET_PROPERTIES (texcav01   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (texcav02   texcav02.cpp)
	ADD_TEST              (texcav02   texcav02)
	TARGET_LINK_LIBRARIES (texcav02   ${LIBS})
	SET_TARGET_PROPERTIES (texcav02   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tlinalg01   tlinalg01.cpp)
	ADD_TEST              (tlinalg01   tlinalg01)
	TARGET_LINK_LIBRARIES (tlinalg01   ${LIBS})
	SET_TARGET_PROPERTIES (tlinalg01   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tcamclay01   tcamclay01.cpp)
	ADD_TEST              (tcamclay01   tcamclay01)
	TARGET_LINK_LIBRARIES (tcamclay01   ${LIBS})
	SET_TARGET_PROPERTIES (tcamclay01   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tbeam01   tbeam01.cpp)
	ADD_TEST              (tbeam01   tbeam01)
	TARGET_LINK_LIBRARIES (tbeam01   ${LIBS})
	SET_TARGET_PROPERTIES (tbeam01   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tbeam02   tbeam02.cpp)
	ADD_TEST              (tbeam02   tbeam02)
	TARGET_LINK_LIBRARIES (tbeam02   ${LIBS})
	SET_TARGET_PROPERTIES (tbeam02   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tspring01   tspring01.cpp)
	ADD_TEST              (tspring01   tspring01)
	TARGET_LINK_LIBRARIES (tspring01   ${LIBS})
	SET_TARGET_PROPERTIES (tspring01   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (ttruss01   ttruss01.cpp)
	ADD_TEST              (ttruss01   ttruss01)
	TARGET_LINK_LIBRARIES (ttruss01   ${LIBS})
	SET_TARGET_PROPERTIES (ttruss01   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (ttruss02   ttruss02.cpp)
	ADD_TEST              (ttruss02   ttruss02)
	TARGET_LINK_LIBRARIES (ttruss02   ${LIBS})
	SET_TARGET_PROPERTIES (ttruss02   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (ttri6      ttri6.cpp)
	ADD_TEST              (ttri6      ttri6)
	TARGET_LINK_LIBRARIES (ttri6      ${LIBS})
	SET_TARGET_PROPERTIES (ttri6      PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tquad4     tquad4.cpp)
	ADD_TEST              (tquad4     tquad4)
	TARGET_LINK_LIBRARIES (tquad4     ${LIBS})
	SET_TARGET_PROPERTIES (tquad4     PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tpstrain01 tpstrain01.cpp)
	ADD_TEST              (tpstrain01 tpstrain01)
	TARGET_LINK_LIBRARIES (tpstrain01 ${LIBS})
	SET_TARGET_PROPERTIES (tpstrain01 PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tpstrain02 tpstrain02.cpp)
	ADD_TEST              (tpstrain02 tpstrain02)
	TARGET_LINK_LIBRARIES (tpstrain02 ${LIBS})
	SET_TARGET_PROPERTIES (tpstrain02 PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tqnorm2d tqnorm2d.cpp)
	ADD_TEST              (tqnorm2d tqnorm2d)
	TARGET_LINK_LIBRARIES (tqnorm2d ${LIBS})
	SET_TARGET_PROPERTIES (tqnorm2d PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tqnorm3d tqnorm3d.cpp)
	ADD_TEST              (tqnorm3d tqnorm3d)
	TARGET_LINK_LIBRARIES (tqnorm3d ${LIBS})
	SET_TARGET_PROPERTIES (tqnorm3d PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tcylin01 tcylin01.cpp)
	ADD_TEST              (tcylin01 tcylin01)
	TARGET_LINK_LIBRARIES (tcylin01 ${LIBS})
	SET_TARGET_PROPERTIES (tcylin01 PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tsolid01   tsolid01.cpp)
	ADD_TEST              (tsolid01   tsolid01)
	TARGET_LINK_LIBRARIES (tsolid01   ${LIBS})
	SET_TARGET_PROPERTIES (tsolid01   PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tpoisson01 tpoisson01.cpp)
	ADD_TEST              (tpoisson01 tpoisson01)
	TARGET_LINK_LIBRARIES (tpoisson01 ${LIBS})
	SET_TARGET_PROPERTIES (tpoisson01 PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tlaplace01 tlaplace01.cpp)
	ADD_TEST              (tlaplace01 tlaplace01)
	TARGET_LINK_LIBRARIES (tlaplace01 ${LIBS})
	SET_TARGET_PROPERTIES (tlaplace01 PROPERTIES COMPILE_FLAGS "${FLAGS}")

ENDIF(BLITZ_FOUND AND LINALG_FOUND)

IF(IGRAPH_FOUND)
	ADD_EXECUTABLE        (ttree01 ttree01.cpp)
	ADD_TEST              (ttree01 ttree01)
	TARGET_LINK_LIBRARIES (ttree01 ${LIBS})
ENDIF(IGRAPH_FOUND)

IF(CGAL_FOUND)
	ADD_EXECUTABLE        (talphashape01 talphashape01.cpp)
	ADD_TEST              (talphashape01 talphashape01)
	TARGET_LINK_LIBRARIES (talphashape01 ${LIBS})
	SET_TARGET_PROPERTIES (talphashape01 PROPERTIES COMPILE_FLAGS "${FLAGS}")
ENDIF(CGAL_FOUND)

SUBDIRS(midasgts reddy)
