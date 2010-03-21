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
  adv_ex_03_01
  ex_01_04
  ex_01_05
  ex_01_06
  ex_04_02
  ex_04_06
  ex_04_08
  ex_04_11
  ex_05_03
  ex_06_22
  ex_07_07
  ex_07_09
  ex_08_01)

FOREACH(var ${TESTS})
    ADD_EXECUTABLE        (${var} "${var}.cpp")
    ADD_TEST              (${var} ${var})
    TARGET_LINK_LIBRARIES (${var} ${LIBS})
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}")
ENDFOREACH(var)
