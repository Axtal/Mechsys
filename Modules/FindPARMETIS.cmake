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

SET(PARMETIS_INCLUDE_SEARCH_PATH
    #$ENV{MECHSYS_ROOT}/pkg/parmetis-4.0.2/include
    #$ENV{HOME}/pkg/parmetis-4.0.2/include
    #$ENV{MECHSYS_ROOT}/pkg/parmetis-4.0.2/metis/include
    #$ENV{HOME}/pkg/parmetis-4.0.2/metis/include
  $ENV{MECHSYS_ROOT}/pkg/ParMetis-3.2.0/
  $ENV{HOME}/pkg/ParMetis-3.2.0/
  $ENV{MECHSYS_ROOT}/pkg/ParMetis-3.2.0/METISLib
  $ENV{HOME}/pkg/ParMetis-3.2.0/METISLib)

SET(PARMETIS_LIBRARY_SEARCH_PATH
    #$ENV{MECHSYS_ROOT}/pkg/parmetis-4.0.2
    #$ENV{HOME}/pkg/parmetis-4.0.2)
  $ENV{MECHSYS_ROOT}/pkg/ParMetis-3.2.0
  $ENV{HOME}/pkg/ParMetis-3.2.0)

FIND_PATH(PARMETIS_METIS_H metis.h ${PARMETIS_INCLUDE_SEARCH_PATH})
FIND_LIBRARY(PARMETIS_METIS NAMES metis PATHS ${PARMETIS_LIBRARY_SEARCH_PATH})

FIND_PATH(PARMETIS_PARMETIS_H parmetis.h ${PARMETIS_INCLUDE_SEARCH_PATH})
FIND_LIBRARY(PARMETIS_PARMETIS NAMES parmetis PATHS ${PARMETIS_LIBRARY_SEARCH_PATH})

SET(PARMETIS_FOUND 1)
FOREACH(var PARMETIS_METIS_H PARMETIS_PARMETIS_H PARMETIS_METIS PARMETIS_PARMETIS)
  IF(NOT ${var})
	SET(PARMETIS_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(PARMETIS_FOUND)
  SET(PARMETIS_INCLUDE_DIRS ${PARMETIS_METIS_H} ${PARMETIS_PARMETIS_H})
  SET(PARMETIS_LIBRARIES    ${PARMETIS_METIS} ${PARMETIS_PARMETIS})
ENDIF(PARMETIS_FOUND)
