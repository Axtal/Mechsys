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

	ADD_EXECUTABLE        (tstatic121 tstatic121.cpp)
	ADD_TEST              (tstatic121 tstatic121)
	TARGET_LINK_LIBRARIES (tstatic121 ${LIBS})
	SET_TARGET_PROPERTIES (tstatic121 PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (tstatic122 tstatic122.cpp)
	ADD_TEST              (tstatic122 tstatic122)
	TARGET_LINK_LIBRARIES (tstatic122 ${LIBS})
	SET_TARGET_PROPERTIES (tstatic122 PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (texam1 texam1.cpp)
	ADD_TEST              (texam1 texam1)
	TARGET_LINK_LIBRARIES (texam1 ${LIBS})
	SET_TARGET_PROPERTIES (texam1 PROPERTIES COMPILE_FLAGS "${FLAGS}")

	ADD_EXECUTABLE        (texam5 texam5.cpp)
	ADD_TEST              (texam5 texam5)
	TARGET_LINK_LIBRARIES (texam5 ${LIBS})
	SET_TARGET_PROPERTIES (texam5 PROPERTIES COMPILE_FLAGS "${FLAGS}")

ENDIF(BLITZ_FOUND AND LINALG_FOUND)
