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
#include <cfloat>  // for DBL_EPSILON

// MechSys
#include "fem/node.h"
#include "fem/element.h"
#include "util/array.h"
#include "util/numstreams.h"
#include "util/exception.h"
#include "mesh/structured.h"

#ifndef VTU_NEWLINE_DEFINED
  #define VTU_NEWLINE_DEFINED
  #define VTU_NEWLINE(I,K,N,KMAX,OF) if (K>KMAX) { OF<<(I<N-1?"\n        ":"\n"); K=0; } else if (I==N-1) { OF<<"\n"; }
#endif

namespace FEM
{

/* Geometry */
class Geom
{
public:
	Geom (int nDim) : _dim(nDim) {}
	void              SetNNodes (size_t nNodes)                                   { for (size_t i=0; i<_nodes.Size(); ++i) { if (_nodes[i]!=NULL) delete _nodes[i]; } _nodes.Resize(nNodes); _nodes = NULL; }
	void              SetNElems (size_t nElems)                                   { for (size_t i=0; i<_elems.Size(); ++i) { if (_elems[i]!=NULL) delete _elems[i]; } _elems.Resize(nElems); _elems = NULL; }
	size_t            NNodes    () const                                          { return _nodes.Size(); }
	size_t            NElems    () const                                          { return _elems.Size(); }
	Node            * SetNode   (size_t i, double X, double Y, double Z=0.0)      { if (_nodes[i]==NULL) { _nodes[i] = new Node;           } _nodes[i]->Initialize (i,X,Y,Z);  return _nodes[i]; }
	Element         * SetElem   (size_t i, char const * Type, bool IsActive=true) { if (_elems[i]==NULL) { _elems[i] = AllocElement(Type); } _elems[i]->SetID      (i); _elems[i]->SetDim (_dim); if (IsActive) _elems[i]->Activate(); else _elems[i]->Deactivate(); return _elems[i]; }
	Node            * Nod       (size_t i)                                        { return _nodes[i]; }
	Element         * Ele       (size_t i)                                        { return _elems[i]; }
	Node      const * Nod       (size_t i) const                                  { return _nodes[i]; }
	Element   const * Ele       (size_t i) const                                  { return _elems[i]; }
	Array<Node*>    & Nodes     ()                                                { return _nodes;    }
	Array<Element*> & Elems     ()                                                { return _elems;    }
private:
	int             _dim;
	Array<Node*>    _nodes;
	Array<Element*> _elems;
}; // class Geom

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

inline void AddNodesElems (Mesh::Structured const * MStruct, char const * ElementType, FEM::Geom * G)
{
	// Nodes
	size_t nn = MStruct->Verts().Size();
	G->SetNNodes (nn);
	for (size_t i=0; i<nn; ++i)
		G->SetNode (i, MStruct->Verts()[i]->C(0), MStruct->Verts()[i]->C(1), (MStruct->Is3D() ? MStruct->Verts()[i]->C(2) : 0.0));

	// Elements
	size_t ne = MStruct->Elems().Size();
	G->SetNElems (ne);
	for (size_t i=0; i<ne; ++i)
	{
		// New elements
		FEM::Element * e = G->SetElem (i, ElementType);

		// Connectivity
		Mesh::Elem * me = MStruct->Elems()[i];
		for (size_t j=0; j<me->V.Size(); ++j)
			e->SetNode (j, G->Nod(me->V[j]->MyID));
	}
}

inline void SetNodeBrys (Mesh::Structured const * MStruct, Array<double> const * X, Array<double> const * Y, Array<double> const * Z, Array<char const *> const * Vars, Array<double> const * Values, FEM::Geom * G, double DistTol=sqrt(DBL_EPSILON))
{
	/* Ex.:
	 *                Coords   Tags    Vars   Values
	 *       [0]       x y z    -10    "ux"      0.0
	 *       [1]       x y z    -20    "fy"    222.0
	 *                 ...                    
	 *    [nNodes-1]   x y z    -30    "uy"      0.0
	 */
	for (size_t i=0; i<MStruct->VertsBry().Size(); ++i) // loop over elements on boundary
	{
		Mesh::Vertex * mv = MStruct->VertsBry()[i];
		for (size_t j=0; j<X->Size(); ++j)
		{
			double d = sqrt(pow((*X)[j]-mv->C(0),2.0) + pow((*Y)[j]-mv->C(1),2.0) + (MStruct->Is3D() ? pow((*Z)[j]-mv->C(2),2.0) : 0.0));
			if (d<DistTol)
			{
				FEM::Node * n = G->Nod(mv->MyID);
				n->Bry ((*Vars)[j], (*Values)[j]);
			}
		}
	}
}

inline void SetFaceBrys (Mesh::Structured const * MStruct, Array<int> const * Tags, Array<char const *> const * Vars, Array<double> const * Values, FEM::Geom * G)
{
	/* Ex.:
	 *                 Tags    Vars   Values
	 *         [0]      -10    "ux"      0.0
	 *         [1]      -20    "fy"    222.0
	 *         ...                    
	 *      [nTags-1]   -30    "uy"      0.0
	 */
	for (size_t i=0; i<MStruct->ElemsBry().Size(); ++i) // loop over elements on boundary
	{
		Mesh::Elem * me = MStruct->ElemsBry()[i];
		for (int j=0; j<me->ETags.Size(); ++j) // j is the local_edge_id
		{
			int tag = me->ETags(j);
			if (tag<0)
			{
				int idx = Tags->Find(tag);
				if (idx>=0)
				{
					FEM::Element * e = G->Ele(me->MyID);
					e->Bry ((*Vars)[idx], (*Values)[idx], j);
				}
				else throw new Fatal("FEM::ApplyEdgeBry: Could not find tag==%d inside Tags array. This tag is set in Elem.ID==%d",tag,me->MyID);
			}
		}
	}
}

inline void _write_elem_val (FEM::Geom const & G, size_t ne, size_t nfmax, char const * Key, std::ostringstream & oss)
{
	oss << "        <DataArray type=\"Float32\" Name=\""<< Key <<"\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	size_t k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		double val = 0.0;
		try { val = G.Ele(i)->Val(Key); } catch (Exception * e) { delete e; }
		oss << (k==0?"  ":" ") << val;
		k++;
		VTU_NEWLINE (i,k,ne,nfmax,oss);
	}
	oss << "        </DataArray>\n";
}

inline void WriteVTUEquilib (FEM::Geom const & G, char const * FileName)
{
	// Open File
	std::ofstream      of(FileName, std::ios::out);
	std::ostringstream oss;

	// Data
	size_t nn = G.NNodes(); // Number of Nodes
	size_t ne = G.NElems(); // Number of Elements

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
		oss << "  " << nsflo << G.Nod(i)->X() << " ";
		oss <<         nsflo << G.Nod(i)->Y() << " ";
		oss <<         nsflo << G.Nod(i)->Z();
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
		String con;  G.Ele(i)->VTKConnect(con);
		oss << "  " << con;
		k++;
		VTU_NEWLINE (i,k,ne,nimax/G.Ele(i)->nNodes(),oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	size_t ossfset = 0;
	for (size_t i=0; i<ne; ++i)
	{
		ossfset += G.Ele(i)->nNodes();
		oss << (k==0?"  ":" ") << ossfset;
		k++;
		VTU_NEWLINE (i,k,ne,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		oss << (k==0?"  ":" ") << G.Ele(i)->VTKCellType();
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
		oss << " " << nsflo << (G.Nod(i)->HasVar("ux") ? G.Nod(i)->Val("ux"): 0.0) << " ";
		oss <<        nsflo << (G.Nod(i)->HasVar("uy") ? G.Nod(i)->Val("uy"): 0.0) << " ";
		oss <<        nsflo << (G.Nod(i)->HasVar("uz") ? G.Nod(i)->Val("uz"): 0.0) << " ";
		k++;
		VTU_NEWLINE (i,k,nn,nfmax/3,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "force" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<nn; ++i)
	{
		oss << " " << nsflo << (G.Nod(i)->HasVar("fx") ? G.Nod(i)->Val("fx"): 0.0) << " ";
		oss <<        nsflo << (G.Nod(i)->HasVar("fy") ? G.Nod(i)->Val("fy"): 0.0) << " ";
		oss <<        nsflo << (G.Nod(i)->HasVar("fz") ? G.Nod(i)->Val("fz"): 0.0) << " ";
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
		oss << (k==0?"  ":" ") << G.Ele(i)->IsActive();
		k++;
		VTU_NEWLINE (i,k,ne,nfmax,oss);
	}
	oss << "        </DataArray>\n";
	_write_elem_val(G, ne, nfmax, "Sx", oss);
	_write_elem_val(G, ne, nfmax, "Sy", oss);
	_write_elem_val(G, ne, nfmax, "Sz", oss);
	_write_elem_val(G, ne, nfmax, "Sxy",oss);
	_write_elem_val(G, ne, nfmax, "Syz",oss);
	_write_elem_val(G, ne, nfmax, "Szx",oss);
	_write_elem_val(G, ne, nfmax, "p",  oss);
	_write_elem_val(G, ne, nfmax, "q",  oss);
	_write_elem_val(G, ne, nfmax, "Ex", oss);
	_write_elem_val(G, ne, nfmax, "Ey", oss);
	_write_elem_val(G, ne, nfmax, "Ez", oss);
	_write_elem_val(G, ne, nfmax, "Exy",oss);
	_write_elem_val(G, ne, nfmax, "Eyz",oss);
	_write_elem_val(G, ne, nfmax, "Ezx",oss);
	_write_elem_val(G, ne, nfmax, "Ev", oss);
	_write_elem_val(G, ne, nfmax, "Ed", oss);
	oss << "      </CellData>\n";

	// Bottom
	oss << "    </Piece>\n";
	oss << "  </UnstructuredGrid>\n";
	oss << "</VTKFile>" << std::endl;

	// Write to file
	of << oss.str();
	of.close();
}

inline void WriteVTK (FEM::Geom const & G, char const * FileName)
{
	// Filter elements
	Array<Element const*> act_elems; // Array for active elements
	for (size_t i=0; i<G.NElems(); ++i)
	{
		if (G.Ele(i)->IsActive()) // Only active elements are considered
			act_elems.Push(G.Ele(i));
	}

	// Data
	size_t n_nodes = G.NNodes();       // Number of Nodes
	size_t n_elems = act_elems.Size(); // Number of Elements
	std::map<String, int>  index_map;  // Map to associate labels with indexes

	// Get all possible labels from elements
	for (size_t i_elem=0; i_elem<n_elems; ++i_elem)
	{
		Array<String>  elem_labels;
		act_elems[i_elem]->GetLabels(elem_labels);
		int n_labels = elem_labels.Size();
		for (int j_label=0; j_label<n_labels; ++j_label)
		{
			String & current_label = elem_labels[j_label];
			if (index_map.find(current_label)==index_map.end())
				index_map[current_label] = index_map.size()-1; // add a new entry
		}
	}

	// Collect nodal values
	size_t n_comps = index_map.size();
	LinAlg::Matrix<double> values(n_nodes, n_comps); values.SetValues(0.0);  // Matrix for nodal values collected from elements
	LinAlg::Matrix<size_t> refs  (n_nodes, n_comps); refs  .SetValues(0);    // Matrix for nodal references of variablea
	for (size_t i_elem=0; i_elem<n_elems; ++i_elem)
	{
		LinAlg::Matrix<double> elem_values;
		Array<String>          elem_labels;
		act_elems[i_elem]->OutNodes(elem_values, elem_labels);
		int                    n_labels     = elem_labels.Size();
		int                    n_elem_nodes = act_elems[i_elem]->nNodes();
		for (int j_label=0; j_label<n_labels; ++j_label)
		{
			String & current_label = elem_labels[j_label];
			int index = index_map[current_label];
			for (int j_node=0; j_node<n_elem_nodes; ++j_node)
			{
				int node_number            = act_elems[i_elem]->GetNode(j_node)->GetID();
				values(node_number,index) += elem_values(j_node, j_label);  // accumulate values
				refs  (node_number,index) ++;                               // update references number
			}
		}
	}
	
	// Compute average values
	for (size_t i=0; i<n_nodes; i++)
	for (size_t j=0; j<n_comps; j++)
	{
		if   (refs(i,j)!=0) values(i,j) /= refs(i,j);
		else                values(i,j)  = 0.0;
	}
	
	// Total number of CELLS data = Sum (1 + nNodes); 1 for the numPts label
	int n_data = 0;
	for (size_t i=0; i<n_elems; ++i) n_data += 1 + act_elems[i]->nNodes();

	// Define variables for displacements
	const String UX = "ux";
	const String UY = "uy";
	const String UZ = "uz";
	
	// Define variables for velocity
	const String VX = "vx";
	const String VY = "vy";
	const String VZ = "vz";

	// Check if exists data about displacements
	bool has_disp = false;
	if (index_map.find(UX)!=index_map.end()) has_disp = true;
	if (index_map.find(UY)!=index_map.end()) has_disp = true;
	if (index_map.find(UZ)!=index_map.end()) has_disp = true;

	// Check if exists data about velocities
	bool has_vel = false;
	if (index_map.find(VX)!=index_map.end()) has_vel  = true;
	if (index_map.find(VY)!=index_map.end()) has_vel  = true;
	if (index_map.find(VZ)!=index_map.end()) has_vel  = true;

	// Structure for output and float number output format
	std::ostringstream oss;
	Util::NumStream nsflo = Util::_8s; // number format for floats

	// Write Legacy VTK file header
	oss << "# vtk DataFile Version 3.0"     << std::endl;
	oss << "MechSys/FEM - "                 << std::endl;
	oss << "ASCII"                          << std::endl;
	oss << "DATASET UNSTRUCTURED_GRID"      << std::endl;
	oss <<                                     std::endl;

	// Node coordinates
	oss << "POINTS " << n_nodes << " float" << std::endl;
	for (size_t i=0; i<n_nodes; ++i)
		oss << nsflo << G.Nod(i)->X() << nsflo << G.Nod(i)->Y() << nsflo << G.Nod(i)->Z() << std::endl;
	oss << std::endl;

	// Elements connectivities
	oss << "CELLS "<< n_elems << " " << n_data << std::endl;
	for (size_t i=0; i<n_elems; ++i)
	{
		String connect; act_elems[i]->VTKConnect(connect);
		oss << act_elems[i]->nNodes() << " " << connect << std::endl;
	}
	oss << std::endl;

	// Cell types
	oss << "CELL_TYPES " << n_elems << std::endl;
	for (size_t i=0; i<n_elems; ++i)
		oss << act_elems[i]->VTKCellType() << std::endl;
	oss << std::endl;

	oss << "POINT_DATA " << n_nodes << std::endl;

	// Vectors
	if (has_disp)
	{
		oss << "VECTORS " << "Disp float" << std::endl;
		for (size_t j=0; j<n_nodes; ++j)
		{
			oss << nsflo << values(j, index_map[UX]) << nsflo << values(j, index_map[UY])
			    << nsflo << ((index_map.find(UZ)==index_map.end())?0.0:values(j, index_map[UZ])) << std::endl;
		}
		oss << std::endl;
	}
	if (has_vel)
	{
		oss << "VECTORS " << "Vel float" << std::endl;
		for (size_t j=0; j<n_nodes; ++j)
		{
			oss << nsflo << values(j, index_map[VX]) << nsflo << values(j, index_map[VY])
			    << nsflo << ((index_map.find(VZ)==index_map.end())?0.0:values(j, index_map[VZ])) << std::endl; 
		}
		oss << std::endl;
	}

	// Scalars
	std::map<String,int>::iterator iter;
	for (iter=index_map.begin(); iter!=index_map.end(); iter++)
	{
		oss << "SCALARS " << iter->first << " float 1" << std::endl;
		oss << "LOOKUP_TABLE default" << std::endl;
		for (size_t j=0; j<n_nodes; ++j)
			oss << nsflo << values(j,iter->second) << std::endl;
		oss << std::endl;
	}

	// Create file and copy contents of 'oss' into it
	std::ofstream ofile;
	ofile.open (FileName, std::ios::out);
	ofile << oss.str();
	ofile.close();
}


}; // namespace FEM

#endif // MECHSYS_FEM_DATA_H
