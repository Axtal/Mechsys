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

SET(EXES
  test_distances
  test_domain
  test_27cubes
  test_tets
  test_dynamics
  test_voro01
  test_voro02
  test_voro03
  triaxial_test
  cylinder_test
  figure)

SET(TESTS
  test_distances)
  #test_domain
  #test_dynamics)

IF(BLITZ_FOUND AND VORO_FOUND)
	FOREACH(var ${EXES})
		ADD_EXECUTABLE        (${var} "${var}.cpp")
		TARGET_LINK_LIBRARIES (${var} ${LIBS})
	ENDFOREACH(var)
ENDIF(BLITZ_FOUND AND VORO_FOUND)

IF(BLITZ_FOUND AND VORO_FOUND)
	FOREACH(var ${TESTS})
        ADD_TEST (${var} ${var})
	ENDFOREACH(var)
ENDIF(BLITZ_FOUND AND VORO_FOUND)
