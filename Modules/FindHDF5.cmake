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

SET(HDF5_INCLUDE_SEARCH_PATH
  /usr/include
  /usr/local/include)

SET(HDF5_LIBRARY_SEARCH_PATH
  /usr/lib
  /usr/local/lib)

FIND_PATH(HDF5_HDF5_H    hdf5.h    ${HDF5_INCLUDE_SEARCH_PATH})
FIND_PATH(HDF5_HDF5_HL_H hdf5_hl.h ${HDF5_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(HDF5_HDF5    NAMES hdf5    PATHS ${HDF5_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(HDF5_HDF5_HL NAMES hdf5_hl PATHS ${HDF5_LIBRARY_SEARCH_PATH})

SET(HDF5_FOUND 1)
FOREACH(var HDF5_HDF5_H HDF5_HDF5_HL_H)
    IF(NOT ${var})
        SET(HDF5_FOUND 0)
    ENDIF(NOT ${var})
ENDFOREACH(var)

IF(HDF5_FOUND)
    SET(HDF5_INCLUDE_DIRS ${HDF5_HDF5_H} ${HDF5_HDF5_HL_H})
    SET(HDF5_LIBRARIES    ${HDF5_HDF5} ${HDF5_HDF5_HL})
ENDIF(HDF5_FOUND)
