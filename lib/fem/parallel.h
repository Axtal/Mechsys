/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef MECHSYS_FEM_PARALLEL_H
#define MECHSYS_FEM_PARALLEL_H

namespace FEM
{

/* Structure with data for parallel computation */
struct PData
{
	int         MyID;       ///< ID (rank) of this processor
	int         nProcs;     ///< The number of processors
	int         nDOFs;      ///< Current total number of DOFs
	int         nDOFsMin;   ///< Mininum total number of DOFs
	int         nDOFsMax;   ///< Maximum total number of DOFs
	int         MyNumEqs;   ///< The number of equations of this processor
	int         MyMinEq;    ///< The mininum equation ID of this processor
	int         MyMaxEq;    ///< The maximum equation ID of this processor
	Array<int>  MyElements; ///< The ID of the elements of this processor
	Array<int>  AllNumEqs;  ///< The number of equations of all processors
	Array<int>  AllMinEq;   ///< All min equation ID of all processors
	Array<int>  OutMyElems; ///< Indexes inside MyElements of the elemens to output
}; // struct PData

}; // namespace FEM

#endif // MECHSYS_FEM_PARALLEL_H
