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

SET(MECHSYS_INCLUDE_SEARCH_PATH
  $ENV{HOME}/mechsys)

SET(MECHSYS_MODULES_SEARCH_PATH
  $ENV{HOME}/mechsys)

FIND_PATH(MECHSYS_MECHSYS_H mechsys/gui/wxmyapp.h   ${MECHSYS_INCLUDE_SEARCH_PATH})
FIND_PATH(MECHSYS_MODULES   Modules/FindBLITZ.cmake ${MECHSYS_MODULES_SEARCH_PATH})

SET(MECHSYS_FOUND 1)
FOREACH(var MECHSYS_MECHSYS_H MECHSYS_MODULES)
  IF(NOT ${var})
	SET(MECHSYS_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(MECHSYS_FOUND)
  SET(MECHSYS_INCLUDE_DIRS ${MECHSYS_MECHSYS_H})
  SET(MECHSYS_SOURCE_DIR   ${MECHSYS_MODULES_SEARCH_PATH})
ENDIF(MECHSYS_FOUND)
