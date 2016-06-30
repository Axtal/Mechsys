########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2009 Sergio Galindo                                    #
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

if (NOT A_USE_HDF5)
    message ("tdem tests need HDF5. Please set A_USE_HDF5 on")
endif (NOT A_USE_HDF5)

SET(EXES
  nonconvex
  ttt
  test_periodic
  test_distances
  test_tetrahedra
  test_domain
  test_dynamics
  test_mesh
  test_read
  test_write
  test_beam
  cylinder_test
  brazil_test
  figure
  test_01
  test_02
  GSD
)

SET(TESTS
  test_distances)
  #test_domain
  #test_dynamics)

FOREACH(var ${EXES})
    ADD_EXECUTABLE        (${var} "${var}.cpp")
    TARGET_LINK_LIBRARIES (${var} ${LIBS})
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
ENDFOREACH(var)

FOREACH(var ${TESTS})
    ADD_TEST (${var} ${var})
ENDFOREACH(var)

#SUBDIRS(newtests)
