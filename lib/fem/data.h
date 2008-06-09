/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
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

#ifndef MECHSYS_FEM_DATA_H
#define MECHSYS_FEM_DATA_H

// MechSys
#include "util/array.h"
#include "fem/node.h"
#include "fem/element.h"

namespace FEM
{

// Geometry type
int GeometryType = 3; ///< Geometry type:  1:1D, 2:2D(plane-strain), 3:3D, 4:Axis-symmetric, 5:2D(plane-stress)

// Solver constants
int    GFE_nSI   = 1;
double GFE_Resid = 0.0;
int    GME_maxSI = 200;
double GME_DTOL  = 1.0e-2;
double GME_dTini = 0.001;
double GME_mMin  = 0.01;
double GME_mMax  = 10;
double GME_mCoef = 0.7;
double GME_ZTOL  = 1.0e-5;
bool   GME_Cconv = true;   ///< Check convergence ?
double GME_Rerr  = 0.0;    ///< Relative error

// Required for parallel processing
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

// Global methods

inline Node * AddNode (double X, double Y, double Z=0.0)
{
	Node * tmp = new Node;
	tmp->Initialize (Nodes.Size(), X, Y, Z);
	Nodes.Push(tmp);
	return tmp;
}

inline Element * AddElem (String const & Type, bool IsActive=true)
{
	Element * tmp = AllocElement(Type);
	tmp->SetID(Elems.Size());
	tmp->SetGeometryType (GeometryType);
	if (IsActive) tmp->Activate  ();
	else          tmp->Deactivate();
	Elems.Push(tmp);
	return tmp;
}

}; // namespace FEM

#endif // MECHSYS_FEM_DATA_H
