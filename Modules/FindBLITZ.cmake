#####################################################################################
# MechSys - A C++ library to simulate (Continuum) Mechanical Systems                #
# Copyright (C) 2005 Dorival de Moraes Pedroso <dorival.pedroso at gmail.com>       #
# Copyright (C) 2005 Raul Dario Durand Farfan  <raul.durand at gmail.com>           #
#                                                                                   #
# This file is part of MechSys.                                                     #
#                                                                                   #
# MechSys is free software; you can redistribute it and/or modify it under the      #
# terms of the GNU General Public License as published by the Free Software         #
# Foundation; either version 2 of the License, or (at your option) any later        #
# version.                                                                          #
#                                                                                   #
# MechSys is distributed in the hope that it will be useful, but WITHOUT ANY        #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A   #
# PARTICULAR PURPOSE. See the GNU General Public License for more details.          #
#                                                                                   #
# You should have received a copy of the GNU General Public License along with      #
# MechSys; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, #
# Fifth Floor, Boston, MA 02110-1301, USA                                           #
#####################################################################################

SET(BLITZ_INCLUDE_SEARCH_PATH
  $ENV{HOME}/pkg/blitz-0.9
  /usr/include
  /usr/local/include)

IF(WIN32)
  SET(BLITZ_INCLUDE_SEARCH_PATH
    ${BLITZ_INCLUDE_SEARCH_PATH}
    C:/msys/1.0/local/include
    D:/msys/1.0/local/include)
ENDIF(WIN32)

FIND_PATH(BLITZ_TINYVEC_ET_H blitz/tinyvec-et.h ${BLITZ_INCLUDE_SEARCH_PATH})
FIND_PATH(BLITZ_TINYMAT_H    blitz/tinymat.h    ${BLITZ_INCLUDE_SEARCH_PATH})

SET(BLITZ_FOUND 1)
FOREACH(var BLITZ_TINYVEC_ET_H BLITZ_TINYMAT_H)
	IF(NOT ${var})
		SET(BLITZ_FOUND 0)
	ENDIF(NOT ${var})
ENDFOREACH(var)

IF(BLITZ_FOUND)
	SET(BLITZ_INCLUDE_DIRS ${BLITZ_TINYVEC_ET_H})
ENDIF(BLITZ_FOUND)
