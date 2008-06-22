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

PROJECT (MECHSYS)

CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

OPTION(MAKE_VERBOSE       "Show additional messages during compilation/link?" OFF)
OPTION(MAKE_ALL_WARNINGS  "Make with all warnings (-Wall)"                    ON )
OPTION(MAKE_DEBUG_SYMBOLS "Make with debug symbols (-g3)"                     ON )
OPTION(MAKE_OPTIMIZED     "Make optimized (-O3)"                              OFF)

IF(MAKE_VERBOSE)
	SET(CMAKE_VERBOSE_MAKEFILE TRUE)
ENDIF(MAKE_VERBOSE)

IF(MAKE_ALL_WARNINGS)
	ADD_DEFINITIONS(-Wall)
ENDIF(MAKE_ALL_WARNINGS)

IF(MAKE_DEBUG_SYMBOLS)
	ADD_DEFINITIONS(-g3)
ENDIF(MAKE_DEBUG_SYMBOLS)

IF(MAKE_OPTIMIZED)
	ADD_DEFINITIONS(-O3)
ENDIF(MAKE_OPTIMIZED)

ADD_DEFINITIONS(-fmessage-length=0) # Each error message will appear on a single line; no line-wrapping will be done.

INCLUDE_DIRECTORIES(${MECHSYS_SOURCE_DIR}/lib/)

SUBDIRS (tst)
