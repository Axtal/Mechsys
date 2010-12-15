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

SET(TESTS
    fig_05_02
    fig_05_11
    fig_05_15
    fig_05_17
    fig_05_24)

SET(PROGS
    fig_06_09
    fig_06_35
    fig_06_38
    fig_09_08
    fig_11_01
    fig_11_01_new
    fig_11_04
    fig_11_04_new
    fig_11_19)

FOREACH(var ${TESTS})
    ADD_EXECUTABLE        (${var} "${var}.cpp")
    ADD_TEST              (${var} ${var})
    TARGET_LINK_LIBRARIES (${var} ${LIBS})
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
ENDFOREACH(var)

FOREACH(var ${PROGS})
    ADD_EXECUTABLE        (${var} "${var}.cpp")
    TARGET_LINK_LIBRARIES (${var} ${LIBS})
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
ENDFOREACH(var)
