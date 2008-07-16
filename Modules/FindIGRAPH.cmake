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

SET(IGRAPH_INCLUDE_SEARCH_PATH
  /usr/include
  /usr/local/include)

SET(IGRAPH_LIBRARY_SEARCH_PATH
  /usr/lib
  /usr/local/lib)

FIND_PATH(IGRAPH_IGRAPH_H igraph/igraph.h ${IGRAPH_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(IGRAPH_IGRAPH NAMES igraph PATHS ${IGRAPH_LIBRARY_SEARCH_PATH})

SET(IGRAPH_FOUND 1)
FOREACH(var IGRAPH_IGRAPH_H IGRAPH_IGRAPH)
  IF(NOT ${var})
	SET(IGRAPH_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(IGRAPH_FOUND)
  SET(IGRAPH_INCLUDE_DIRS ${IGRAPH_IGRAPH_H})
  SET(IGRAPH_LIBRARIES    ${IGRAPH_IGRAPH})
ENDIF(IGRAPH_FOUND)
