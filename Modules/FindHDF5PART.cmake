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

SET(HDF5PART_INCLUDE_SEARCH_PATH
  $ENV{HOME}/pkg/H5Part-1.6.0)

SET(HDF5PART_LIBRARY_SEARCH_PATH
  $ENV{HOME}/pkg/H5Part-1.6.0/src/)

FIND_PATH(HDF5PART_H    src/H5Part.h     ${HDF5PART_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(HDF5PART    NAMES H5Part    PATHS ${HDF5PART_LIBRARY_SEARCH_PATH})

SET(HDF5PART_FOUND 1)
FOREACH(var HDF5PART_H HDF5PART )
  IF(NOT ${var})
	SET(HDF5PART_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(HDF5PART_FOUND)
  SET(HDF5PART_INCLUDE_DIRS ${HDF5PART_H})
  SET(HDF5PART_LIBRARIES    ${HDF5PART})
ENDIF(HDF5PART_FOUND)
