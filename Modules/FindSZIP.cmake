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

SET(SZIP_LIBRARY_SEARCH_PATH
  $ENV{HOME}/pkg/szip-2.1/szip/lib/)

FIND_LIBRARY(SZIP        NAMES sz        PATHS ${SZIP_LIBRARY_SEARCH_PATH})

SET(SZIP_FOUND 1)
FOREACH(var SZIP)
  IF(NOT ${var})
	SET(SZIP_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(SZIP_FOUND)
  SET(SZIP_LIBRARIES    ${SZIP} ${SZIP})
ENDIF(SZIP_FOUND)
