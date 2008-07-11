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

#ifndef MECHSYS_FEM_FUNCTIONS_H
#define MECHSYS_FEM_FUNCTIONS_H

// STL
#include <iostream>
#include <fstream>
#include <cfloat>  // for DBL_EPSILON

// MechSys
#include "fem/node.h"
#include "fem/element.h"
#include "fem/geometry.h"
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

inline void AddNodesElems (Mesh::Structured const * M, char const * ElementType, FEM::Geom * G)
{
	// Nodes
	size_t nn = M->Verts().Size();
	G->SetNNodes (nn);
	for (size_t i=0; i<nn; ++i)
		G->SetNode (i, M->Verts()[i]->C(0), M->Verts()[i]->C(1), (M->Is3D() ? M->Verts()[i]->C(2) : 0.0));

	// Elements
	size_t ne = M->Elems().Size();
	G->SetNElems (ne);
	for (size_t i=0; i<ne; ++i)
	{
		// New elements
		FEM::Element * e = G->SetElem (i, ElementType);

		// Connectivity
		Mesh::Elem * me = M->Elems()[i];
		for (size_t j=0; j<me->V.Size(); ++j)
			e->Connect (j, G->Nod(me->V[j]->MyID));
	}
}

inline void SetNodeBrys (Mesh::Structured const * M, Array<double> const * X, Array<double> const * Y, Array<double> const * Z, Array<char const *> const * Vars, Array<double> const * Values, FEM::Geom * G, double DistTol=sqrt(DBL_EPSILON))
{
	/* Ex.:
	 *                Coords   Tags    Vars   Values
	 *       [0]       x y z    -10    "ux"      0.0
	 *       [1]       x y z    -20    "fy"    222.0
	 *                 ...                    
	 *    [nNodes-1]   x y z    -30    "uy"      0.0
	 */
	for (size_t i=0; i<M->VertsBry().Size(); ++i) // loop over elements on boundary
	{
		Mesh::Vertex * mv = M->VertsBry()[i];
		for (size_t j=0; j<X->Size(); ++j)
		{
			double d = sqrt(pow((*X)[j]-mv->C(0),2.0) + pow((*Y)[j]-mv->C(1),2.0) + (M->Is3D() ? pow((*Z)[j]-mv->C(2),2.0) : 0.0));
			if (d<DistTol)
			{
				FEM::Node * n = G->Nod(mv->MyID);
				n->Bry ((*Vars)[j], (*Values)[j]);
			}
		}
	}
}

inline void SetFaceBrys (Mesh::Structured const * M, Array<int> const * Tags, Array<char const *> const * Vars, Array<double> const * Values, FEM::Geom * G)
{
	/* Ex.:
	 *                 Tags    Vars   Values
	 *         [0]      -10    "ux"      0.0
	 *         [1]      -20    "fy"    222.0
	 *         ...                    
	 *      [nTags-1]   -30    "uy"      0.0
	 */
	for (size_t i=0; i<M->ElemsBry().Size(); ++i) // loop over elements on boundary
	{
		Mesh::Elem * me = M->ElemsBry()[i];
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
	size_t nn = G.nNodes(); // Number of Nodes
	size_t ne = G.nElems(); // Number of Elements

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
	for (size_t i=0; i<G.nElems(); ++i)
	{
		if (G.Ele(i)->IsActive()) // Only active elements are considered
			act_elems.Push(G.Ele(i));
	}

	// Data
	size_t n_nodes = G.nNodes();       // Number of Nodes
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
				int node_number            = act_elems[i_elem]->Nod(j_node)->GetID();
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


#ifdef USE_BOOST_PYTHON
// {

namespace boopy = boost::python;

void PyAddNodesElems   (PyMeshStruct const & MS, boopy::str const & ElemType, PyGeom & G) { FEM::AddNodesElems (MS.GetMesh(), boopy::extract<char const *>(ElemType)(), G.GetGeom()); }
void PySetNodeBrys     (PyMeshStruct const & MS, boopy::list & X, boopy::list & Y, boopy::list & Z, boopy::list & Vars, boopy::list & Values, PyGeom & G)
{
	for (size_t i=0; i<MS.GetMesh()->VertsBry().Size(); ++i) // loop over elements on boundary
	{
		Mesh::Vertex * mv = MS.GetMesh()->VertsBry()[i];
		for (int j=0; j<len(X); ++j)
		{
			double x = boopy::extract<double>(X[j])();
			double y = boopy::extract<double>(Y[j])();
			double z = boopy::extract<double>(Z[j])();
			double d = sqrt(pow(x-mv->C(0),2.0) + pow(y-mv->C(1),2.0) + (MS.GetMesh()->Is3D() ? pow(z-mv->C(2),2.0) : 0.0));
			if (d<sqrt(DBL_EPSILON))
			{
				FEM::Node * n = G.GetGeom()->Nod(mv->MyID);
				n->Bry (boopy::extract<char const *>(Vars[j])(), boopy::extract<double>(Values[j])());
			}
		}
	}
}
void PySetFaceBrys     (PyMeshStruct const & MS, boopy::list & Tags, boopy::list & Vars, boopy::list & Values, PyGeom & G)
{
	for (size_t i=0; i<MS.GetMesh()->ElemsBry().Size(); ++i) // loop over elements on boundary
	{
		Mesh::Elem * me = MS.GetMesh()->ElemsBry()[i];
		for (int j=0; j<me->ETags.Size(); ++j) // j is the local_edge_id
		{
			int tag = me->ETags(j);
			if (tag<0)
			{
				long idx = Tags.index(tag);
				if (idx>=0)
				{
					FEM::Element * e = G.GetGeom()->Ele(me->MyID);
					e->Bry (boopy::extract<char const*>(Vars[idx])(), boopy::extract<double>(Values[idx])(), j);
				}
				else throw new Fatal("PySetFaceBrys: Could not find tag==%d inside Tags array. This tag is set in Elem.ID==%d",tag,me->MyID);
			}
		}
	}
}
void PyWriteVTUEquilib (PyGeom const & G, boopy::str const & FileName) { FEM::WriteVTUEquilib((*G.GetGeom()), boopy::extract<char const *>(FileName)()); }
void PyWriteVTK        (PyGeom const & G, boopy::str const & FileName) { FEM::WriteVTK       ((*G.GetGeom()), boopy::extract<char const *>(FileName)()); }

void PySetGeom (PyMeshStruct const & M, boopy::list const & NodesBrys, boopy::list const & FacesBrys, boopy::list const & ElemsAtts, PyGeom & G)
{
	/* Example:
	 *          NodesBrys = [ [(0.,1.,2.), (0.,0.,2.)], ['ux','fy'], [0., -1.] ]
	 *          FacesBrys = [ [-10, -20], ['uy', 'fy'], [0.0, -1] ]
	 *          ElemAtts  = [ [-1, -2], ['Quad4PStrain', 'Tri6PStrain'], ['LinElastic', 'LinElastic'], ['E=10 nu=0.2', 'E=20 nu=0.3'], ['Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0', 'Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0'] ]
	 *   where:
	 *          NodesBrys = [ [(x1,y1,z1), (x2,y2,z2)], [key1, key2], [val1, val2] ]   // (z is optional)
	 *          FacesBrys = [ [tag1, tag2], [key1, key2], [val1, val2] ]
	 *          ElemAtts  = [ [tag1, tag2], [type1, type2], [mdl1, mdl2], [prms1, prms2], [inis1, inis2] ]
	 */

	// 3D mesh?
	bool is3d = M.GetMesh()->Is3D();

	// Extract list with nodes boundaries
	int         nb_nxyz = 0;
	int         nb_nkey = 0;
	int         nb_nval = 0;
	boopy::list nb_xyzs;
	boopy::list nb_keys;
	boopy::list nb_vals;
	if (len(NodesBrys)==3) // 3 sublists: [[x,y,z], [key], [val]]
	{
		nb_nxyz = len(NodesBrys[0]);
		nb_nkey = len(NodesBrys[1]);
		nb_nval = len(NodesBrys[2]);
		if (nb_nxyz==nb_nkey && nb_nkey==nb_nval)
		{
			nb_xyzs = boopy::extract<boopy::list>(NodesBrys[0])();
			nb_keys = boopy::extract<boopy::list>(NodesBrys[1])();
			nb_vals = boopy::extract<boopy::list>(NodesBrys[2])();
		}
		else throw new Fatal("PySetGeom: Each sublist of NodesBrys list must have the same number of items.\n\tExample: [ [(x1,y1,z1), (x2,y2,z2)] , [key1, key2] , [val1, val2] ]  where (z is optional).");
	}
	else if (len(NodesBrys)!=0) throw new Fatal("PySetGeom: The list with node boundaries must have 3 sublists.\n\tExample: [ [(x1,y1,z1), (x2,y2,z2)], [key1, key2], [val1, val2] ]  where (z is optional).");

	// Extract list with faces boundaries
	int         fb_ntag = 0;
	int         fb_nkey = 0;
	int         fb_nval = 0;
	boopy::list fb_tags;
	boopy::list fb_keys;
	boopy::list fb_vals;
	if (len(FacesBrys)==3) // 3 sublists: [[tag], [key], [val]]
	{
		fb_ntag = len(FacesBrys[0]);
		fb_nkey = len(FacesBrys[1]);
		fb_nval = len(FacesBrys[2]);
		if (fb_ntag==fb_nkey && fb_nkey==fb_nval)
		{
			fb_tags = boopy::extract<boopy::list>(FacesBrys[0])();
			fb_keys = boopy::extract<boopy::list>(FacesBrys[1])();
			fb_vals = boopy::extract<boopy::list>(FacesBrys[2])();
		}
		else throw new Fatal("PySetGeom: Each sublist of FacesBrys list must have the same number of items.\n\tExample: [ [tag1, tag2] , [key1, key2] , [val1, val2] ].");
	}
	else if (len(FacesBrys)!=0) throw new Fatal("PySetGeom: The list with face boundaries must have 3 sublists.\n\tExample: [ [tag1, tag2], [key1, key2], [val1, val2] ].");

	// Extract list with elements attributes
	int         lb_ntag = 0; // (l) stands for eLement
	int         lb_ntyp = 0;
	int         lb_nmdl = 0;
	int         lb_nprm = 0;
	int         lb_nini = 0;
	boopy::list lb_tags;
	boopy::list lb_typs;
	boopy::list lb_mdls;
	boopy::list lb_prms;
	boopy::list lb_inis;
	if (len(ElemsAtts)==5) // 5 sublists: [[tag], [type], [model], [prms], [inis]]
	{
		lb_ntag = len(ElemsAtts[0]);
		lb_ntyp = len(ElemsAtts[1]);
		lb_nmdl = len(ElemsAtts[2]);
		lb_nprm = len(ElemsAtts[3]);
		lb_nini = len(ElemsAtts[4]);
		if (lb_ntag==lb_ntyp && lb_ntyp==lb_nmdl && lb_nmdl==lb_nprm && lb_nprm==lb_nini)
		{
			lb_tags = boopy::extract<boopy::list>(ElemsAtts[0])();
			lb_typs = boopy::extract<boopy::list>(ElemsAtts[1])();
			lb_mdls = boopy::extract<boopy::list>(ElemsAtts[2])();
			lb_prms = boopy::extract<boopy::list>(ElemsAtts[3])();
			lb_inis = boopy::extract<boopy::list>(ElemsAtts[4])();
		}
		else throw new Fatal("PySetGeom: Each sublist of ElemsAtts list must have the same number of items.\n\tExample: [ [tag1, tag2], [typ1, typ2], [mdl1, mdl2], [prms1, prms2], [inis1, ini2s] ].");
	}
	else throw new Fatal("PySetGeom: The list with elements attributes must have 5 sublists.\n\tExample: [ [tag1, tag2], [typ1, typ2], [mdl1, mdl2], [prms1, prms2], [inis1, inis2] ].");

	// Set nodes
	size_t nn = M.GetMesh()->Verts().Size();
	G.GetGeom()->SetNNodes (nn);
	for (size_t i=0; i<nn; ++i) // loop over all vertices
	{
		// Current vertex
		Mesh::Vertex const * v = M.GetMesh()->Verts()[i];

		// New node
		G.GetGeom()->SetNode (i, v->C(0), v->C(1), (is3d ? v->C(2) : 0.0));
	}

	// Set elements
	size_t ne = M.GetMesh()->Elems().Size();
	G.GetGeom()->SetNElems (ne);
	for (size_t i=0; i<ne; ++i)
	{
		// Current mesh element
		Mesh::Elem const * me = M.GetMesh()->Elems()[i];

		// Set element
		bool found = false;
		for (int j=0; j<lb_ntag; ++j)
		{
			if (me->Tag==boopy::extract<int>(lb_tags[j])())
			{
				// New finite element
				found = true;
				FEM::Element * fe = G.GetGeom()->SetElem (i, boopy::extract<char const *>(lb_typs[j])());

				// Set connectivity
				for (size_t k=0; k<me->V.Size(); ++k)
					fe->Connect (k, G.GetGeom()->Nod(me->V[k]->MyID));

				// Set parameters and initial values
				fe->SetModel (boopy::extract<char const *>(lb_mdls[j])(), boopy::extract<char const *>(lb_prms[j])(), boopy::extract<char const*>(lb_inis[j])());
			}
		}
		if (found==false) throw new Fatal("PySetGeom: Could not find Tag==%d for Element %d in the ElemsAtts list",me->Tag,i);
	}

	// Set faces boundaries
	if (fb_ntag>0)
	{
		for (size_t i=0; i<M.GetMesh()->ElemsBry().Size(); ++i) // loop over all elements on boundary
		{
			// Current mesh element
			Mesh::Elem const * me = M.GetMesh()->ElemsBry()[i];

			for (int j=0; j<me->ETags.Size(); ++j) // j is the local face id
			{
				int tag = me->ETags(j);
				if (tag<0) // this element has a face tag
				{
					bool found = false;
					for (int k=0; k<fb_ntag; ++k)
					{
						if (tag==boopy::extract<int>(fb_tags[k])())
						{
							found = true;
							G.GetGeom()->Ele(me->MyID)->Bry (boopy::extract<char const *>(fb_keys[k])(), boopy::extract<double>(fb_vals[k])(), j);
							break;
						}
					}
					if (found==false) throw new Fatal("PySetGeom: Could not find Tag==%d for Face %d of Element %d in the FacesBrys list",tag,j,me->MyID);
				}
			}
		}
	}

	// Set nodes boundaries
	if (nb_nxyz>0)
	{
		for (size_t i=0; i<M.GetMesh()->VertsBry().Size(); ++i) // loop over all vertices on boundary
		{
			// Current vertex
			Mesh::Vertex const * v = M.GetMesh()->VertsBry()[i];

			for (int j=0; j<nb_nxyz; ++j)
			{
				double x =         boopy::extract<double>(boopy::extract<boopy::tuple>(nb_xyzs[j])()[0])();
				double y =         boopy::extract<double>(boopy::extract<boopy::tuple>(nb_xyzs[j])()[1])();
				double z = (is3d ? boopy::extract<double>(boopy::extract<boopy::tuple>(nb_xyzs[j])()[2])() : 0.0);
				double d = sqrt(pow(x - v->C(0),2.0) + pow(y - v->C(1),2.0) + (is3d ? pow(z - v->C(2),2.0) : 0.0));
				if (d<sqrt(DBL_EPSILON))
					G.GetGeom()->Nod(v->MyID)->Bry (boopy::extract<char const *>(nb_keys[j])(), boopy::extract<double>(nb_vals[j])());
			}
		}
	}

	std::cout << "ok" << std::endl;
}

// }
#endif // USE_BOOST_PYTHON


#endif // MECHSYS_FEM_FUNCTIONS_H
