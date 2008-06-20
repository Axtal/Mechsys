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

// STL
#include <iostream>
#include <fstream>

// MechSys
#include "fem/node.h"
#include "fem/element.h"
#include "util/array.h"
#include "util/numstreams.h"
#include "util/exception.h"

#define VTU_NEWLINE(I,K,N,KMAX,OF) if (K>KMAX) { OF<<(I<N-1?"\n        ":"\n"); K=0; } else if (I==N-1) { OF<<"\n"; }

namespace FEM
{

// Problem
int Dim = 3; ///< Space dimension

// Solver constants
int    GFE_nSI   = 1;      ///< Global FE solver: number of sub-increments
double GFE_Resid = 0.0;    ///< Global FE solver:  
int    GME_maxSI = 200;    ///< Global ME solver:
double GME_DTOL  = 1.0e-2; ///< Global ME solver:
double GME_dTini = 0.001;  ///< Global ME solver:  
double GME_mMin  = 0.01;   ///< Global ME solver:  
double GME_mMax  = 10;     ///< Global ME solver:  
double GME_mCoef = 0.7;    ///< Global ME solver:  
double GME_ZTOL  = 1.0e-5; ///< Global FE solver:
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

inline Element * AddElem (char const * Type, bool IsActive=true)
{
	Element * tmp = AllocElement(Type);
	tmp->SetID  (Elems.Size());
	tmp->SetDim (Dim);
	if (IsActive) tmp->Activate  ();
	else          tmp->Deactivate();
	Elems.Push(tmp);
	return tmp;
}

inline void _write_elem_val (size_t ne, size_t nfmax, char const * Key, std::ostringstream & oss)
{
	oss << "        <DataArray type=\"Float32\" Name=\""<< Key <<"\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	size_t k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		double val = 0.0;
		try { val = Elems[i]->Val(Key); } catch (Exception * e) { delete e; }
		oss << (k==0?"  ":" ") << val;
		k++;
		VTU_NEWLINE (i,k,ne,nfmax,oss);
	}
	oss << "        </DataArray>\n";
}

inline void WriteVTUEquilib (char const * FileName)
{
	// Open File
	std::ofstream      of(FileName, std::ios::out);
	std::ostringstream oss;

	// Data
	size_t nn = Nodes.Size(); // Number of Nodes
	size_t ne = Elems.Size(); // Number of Elements

	// Constants
	size_t          nimax = 40;        // number of integers in a line
	size_t          nfmax = 12;        // number of floats in a line
	Util::NumStream nsflo = Util::_8s; // number format for floats

	// Header
	oss << "<?xml version=\"1.0\"?>\n";
	oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	oss << "  <UnstructuredGrid>\n";
	oss << "    <Piece NumberOfPoints=\"" << nn << "\" NumberOfCells=\"" << ne << "\">\n";

	// Nodes: coordinates
	oss << "      <Points>\n";
	oss << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	size_t k = 0; oss << "        ";
	for (size_t i=0; i<nn; ++i)
	{
		oss << "  " << nsflo << Nodes[i]->X() << " ";
		oss <<         nsflo << Nodes[i]->Y() << " ";
		oss <<         nsflo << Nodes[i]->Z();
		k++;
		VTU_NEWLINE (i,k,nn,nfmax/3,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </Points>\n";

	// Elements: connectivity, offsets, types
	oss << "      <Cells>\n";
	oss << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		String con;  Elems[i]->VTKConnect(con);
		oss << "  " << con;
		k++;
		VTU_NEWLINE (i,k,ne,nimax/Elems[i]->nNodes(),oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	size_t ossfset = 0;
	for (size_t i=0; i<ne; ++i)
	{
		ossfset += Elems[i]->nNodes();
		oss << (k==0?"  ":" ") << ossfset;
		k++;
		VTU_NEWLINE (i,k,ne,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		oss << (k==0?"  ":" ") << Elems[i]->VTKCellType();
		k++;
		VTU_NEWLINE (i,k,ne,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </Cells>\n";

	// Data -- nodes
	oss << "      <PointData Vectors=\"TheVectors\">\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "displ" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<nn; ++i)
	{
		oss << " " << nsflo << (Nodes[i]->HasVar("ux") ? Nodes[i]->Val("ux"): 0.0) << " ";
		oss <<        nsflo << (Nodes[i]->HasVar("uy") ? Nodes[i]->Val("uy"): 0.0) << " ";
		oss <<        nsflo << (Nodes[i]->HasVar("uz") ? Nodes[i]->Val("uz"): 0.0) << " ";
		k++;
		VTU_NEWLINE (i,k,nn,nfmax/3,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "force" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<nn; ++i)
	{
		oss << " " << nsflo << (Nodes[i]->HasVar("fx") ? Nodes[i]->Val("fx"): 0.0) << " ";
		oss <<        nsflo << (Nodes[i]->HasVar("fy") ? Nodes[i]->Val("fy"): 0.0) << " ";
		oss <<        nsflo << (Nodes[i]->HasVar("fz") ? Nodes[i]->Val("fz"): 0.0) << " ";
		k++;
		VTU_NEWLINE (i,k,nn,nfmax/3,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </PointData>\n";

	// Data -- elements
	oss << "      <CellData Scalars=\"TheScalars\">\n";
	oss << "        <DataArray type=\"Float32\" Name=\"active\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		oss << (k==0?"  ":" ") << Elems[i]->IsActive();
		k++;
		VTU_NEWLINE (i,k,ne,nfmax,oss);
	}
	oss << "        </DataArray>\n";
	_write_elem_val(ne, nfmax, "Sx", oss);
	_write_elem_val(ne, nfmax, "Sy", oss);
	_write_elem_val(ne, nfmax, "Sz", oss);
	_write_elem_val(ne, nfmax, "Sxy",oss);
	_write_elem_val(ne, nfmax, "Syz",oss);
	_write_elem_val(ne, nfmax, "Szx",oss);
	_write_elem_val(ne, nfmax, "p",  oss);
	_write_elem_val(ne, nfmax, "q",  oss);
	_write_elem_val(ne, nfmax, "Ex", oss);
	_write_elem_val(ne, nfmax, "Ey", oss);
	_write_elem_val(ne, nfmax, "Ez", oss);
	_write_elem_val(ne, nfmax, "Exy",oss);
	_write_elem_val(ne, nfmax, "Eyz",oss);
	_write_elem_val(ne, nfmax, "Ezx",oss);
	_write_elem_val(ne, nfmax, "Ev", oss);
	_write_elem_val(ne, nfmax, "Ed", oss);
	oss << "      </CellData>\n";

	// Bottom
	oss << "    </Piece>\n";
	oss << "  </UnstructuredGrid>\n";
	oss << "</VTKFile>" << std::endl;

	// Write to file
	of << oss.str();
	of.close();
}

}; // namespace FEM

#endif // MECHSYS_FEM_DATA_H
