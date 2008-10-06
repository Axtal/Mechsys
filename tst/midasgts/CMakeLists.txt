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

ENDIF(BLITZ_FOUND AND LINALG_FOUND)
