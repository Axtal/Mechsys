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

SET(LAPACK_INCLUDE_SEARCH_PATH
  $ENV{HOME}/pkg/lapack-3.2.1
  /usr/include
  /usr/local/include)

SET(LAPACK_LIBRARY_SEARCH_PATH
  $ENV{HOME}/pkg/lapack-3.2.1
  /usr/lib
  /usr/local/lib)

FIND_LIBRARY(LAPACK_LAPACK NAMES lapack PATHS ${LAPACK_LIBRARY_SEARCH_PATH})

SET(LAPACK_FOUND 1)
FOREACH(var LAPACK_LAPACK)
  IF(NOT ${var})
	SET(LAPACK_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(LAPACK_FOUND)
  SET(LAPACK_LIBRARIES    ${LAPACK_LAPACK})
ENDIF(LAPACK_FOUND)
