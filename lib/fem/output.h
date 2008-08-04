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

#ifndef MECHSYS_OUTPUT_H
#define MECHSYS_OUTPUT_H

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
#include "mesh/mesh.h"
#include "mesh/structured.h"

#ifndef VTU_NEWLINE_DEFINED
  #define VTU_NEWLINE_DEFINED
  #define VTU_NEWLINE(I,K,N,KMAX,OF) if (K>KMAX) { OF<<(I<N-1?"\n        ":"\n"); K=0; } else if (I==N-1) { OF<<"\n"; }
#endif

namespace Output
{

class VTU
{
public:
	// Constructor
	VTU () : _nimax(40), _nfmax(12), _nn(0), _ne(0), _nsflo(Util::_8s) {}

	// Methods
	void Equilib (FEM::Geom const * G, char const * FileName); ///< Write a ParaView-VTU file for Equilibrium problems
	void Heat    (FEM::Geom const * G, char const * FileName); ///< Write a ParaView-VTU file for Heat problems

private:
	// Data
	size_t                  _nimax; ///< number of integers in a line
	size_t                  _nfmax; ///< number of floats in a line
	size_t                  _nn;    ///< Number of Nodes
	size_t                  _ne;    ///< Number of Elements
	Util::NumStream         _nsflo; ///< number format for floats
	FEM::Geom       const * _g;     ///< Geometry structure

	// Methods
	void _write_header   (                  std::ostringstream & oss) const;
	void _write_geometry (                  std::ostringstream & oss) const;
	void _write_elem_val (char const * Key, std::ostringstream & oss) const;
	void _write_bottom   (char const * FN,  std::ostringstream & oss) const;

}; // class VTU


inline void VTK (FEM::Geom const & G, char const * FileName)
{
	// Filter elements
	Array<FEM::Element const*> act_elems; // Array for active elements
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
		int                    n_elem_nodes = act_elems[i_elem]->NNodes();
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
	
	// Total number of CELLS data = Sum (1 + NNodes); 1 for the numPts label
	int n_data = 0;
	for (size_t i=0; i<n_elems; ++i) n_data += 1 + act_elems[i]->NNodes();

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
		oss << act_elems[i]->NNodes() << " " << connect << std::endl;
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


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void VTU::Equilib(FEM::Geom const * G, char const * FileName)
{
	// Set data
	_nn = G->NNodes(); // Number of Nodes
	_ne = G->NElems(); // Number of Elements
	_g  = G;           // geometry

	// Output structure
	std::ostringstream oss;

	// Header
	_write_header (oss);

	// Geometry
	_write_geometry (oss);

	// Data -- nodes
	oss << "      <PointData Vectors=\"TheVectors\">\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "displ" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	size_t k = 0; oss << "        ";
	for (size_t i=0; i<_nn; ++i)
	{
		oss << " " << _nsflo << (_g->Nod(i)->HasVar("ux") ? _g->Nod(i)->Val("ux"): 0.0) << " ";
		oss <<        _nsflo << (_g->Nod(i)->HasVar("uy") ? _g->Nod(i)->Val("uy"): 0.0) << " ";
		oss <<        _nsflo << (_g->Nod(i)->HasVar("uz") ? _g->Nod(i)->Val("uz"): 0.0) << " ";
		k++;
		VTU_NEWLINE (i,k,_nn,_nfmax/3,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "force" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<_nn; ++i)
	{
		oss << " " << _nsflo << (_g->Nod(i)->HasVar("fx") ? _g->Nod(i)->Val("fx"): 0.0) << " ";
		oss <<        _nsflo << (_g->Nod(i)->HasVar("fy") ? _g->Nod(i)->Val("fy"): 0.0) << " ";
		oss <<        _nsflo << (_g->Nod(i)->HasVar("fz") ? _g->Nod(i)->Val("fz"): 0.0) << " ";
		k++;
		VTU_NEWLINE (i,k,_nn,_nfmax/3,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </PointData>\n";

	// Data -- elements
	oss << "      <CellData Scalars=\"TheScalars\">\n";
	oss << "        <DataArray type=\"Float32\" Name=\"active\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<_ne; ++i)
	{
		oss << (k==0?"  ":" ") << _g->Ele(i)->IsActive();
		k++;
		VTU_NEWLINE (i,k,_ne,_nfmax,oss);
	}
	oss << "        </DataArray>\n";
	_write_elem_val("Sx", oss);
	_write_elem_val("Sy", oss);
	_write_elem_val("Sz", oss);
	_write_elem_val("Sxy",oss);
	_write_elem_val("Syz",oss);
	_write_elem_val("Szx",oss);
	_write_elem_val("p",  oss);
	_write_elem_val("q",  oss);
	_write_elem_val("Ex", oss);
	_write_elem_val("Ey", oss);
	_write_elem_val("Ez", oss);
	_write_elem_val("Exy",oss);
	_write_elem_val("Eyz",oss);
	_write_elem_val("Ezx",oss);
	_write_elem_val("Ev", oss);
	_write_elem_val("Ed", oss);
	oss << "      </CellData>\n";

	// Bottom and write to file
	_write_bottom (FileName, oss);
}

inline void VTU::Heat(FEM::Geom const * G, char const * FileName)
{
	// Set data
	_nn = G->NNodes(); // Number of Nodes
	_ne = G->NElems(); // Number of Elements
	_g  = G;           // geometry

	// Output structure
	std::ostringstream oss;

	// Header
	_write_header (oss);

	// Geometry
	_write_geometry (oss);

	// Data -- nodes
	oss << "      <PointData Scalars=\"TheScalars\">\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "temp" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	size_t k = 0; oss << "        ";
	for (size_t i=0; i<_nn; ++i)
	{
		oss << (k==0?"  ":" ") << _g->Nod(i)->Val("T");
		k++;
		VTU_NEWLINE (i,k,_nn,_nfmax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "heat" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<_nn; ++i)
	{
		oss << (k==0?"  ":" ") << _g->Nod(i)->Val("F");
		k++;
		VTU_NEWLINE (i,k,_nn,_nfmax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </PointData>\n";

	// Bottom and write to file
	_write_bottom (FileName, oss);
}


/* privete */

inline void VTU::_write_header(std::ostringstream & oss) const
{
	// Header
	oss << "<?xml version=\"1.0\"?>\n";
	oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	oss << "  <UnstructuredGrid>\n";
	oss << "    <Piece NumberOfPoints=\"" << _nn << "\" NumberOfCells=\"" << _ne << "\">\n";
}

inline void VTU::_write_geometry(std::ostringstream & oss) const
{
	// Nodes: coordinates
	oss << "      <Points>\n";
	oss << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	size_t k = 0; oss << "        ";
	for (size_t i=0; i<_nn; ++i)
	{
		oss << "  " << _nsflo << _g->Nod(i)->X() << " ";
		oss <<         _nsflo << _g->Nod(i)->Y() << " ";
		oss <<         _nsflo << _g->Nod(i)->Z();
		k++;
		VTU_NEWLINE (i,k,_nn,_nfmax/3,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </Points>\n";

	// Elements: connectivity, offsets, types
	oss << "      <Cells>\n";
	oss << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<_ne; ++i)
	{
		String con;  _g->Ele(i)->VTKConnect(con);
		oss << "  " << con;
		k++;
		VTU_NEWLINE (i,k,_ne,_nimax/_g->Ele(i)->NNodes(),oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	size_t offset = 0;
	for (size_t i=0; i<_ne; ++i)
	{
		offset += _g->Ele(i)->NNodes();
		oss << (k==0?"  ":" ") << offset;
		k++;
		VTU_NEWLINE (i,k,_ne,_nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<_ne; ++i)
	{
		oss << (k==0?"  ":" ") << _g->Ele(i)->VTKCellType();
		k++;
		VTU_NEWLINE (i,k,_ne,_nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </Cells>\n";
}

inline void VTU::_write_elem_val(char const * Key, std::ostringstream & oss) const
{
	oss << "        <DataArray type=\"Float32\" Name=\""<< Key <<"\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	size_t k = 0; oss << "        ";
	for (size_t i=0; i<_ne; ++i)
	{
		double val = 0.0;
		try { val = _g->Ele(i)->Val(Key); } catch (Exception * e) { delete e; }
		oss << (k==0?"  ":" ") << val;
		k++;
		VTU_NEWLINE (i,k,_ne,_nfmax,oss);
	}
	oss << "        </DataArray>\n";
}

inline void VTU::_write_bottom(char const * FN, std::ostringstream & oss) const
{
	// Bottom
	oss << "    </Piece>\n";
	oss << "  </UnstructuredGrid>\n";
	oss << "</VTKFile>" << std::endl;

	// Write to file
	std::ofstream of(FN, std::ios::out);
	of << oss.str();
	of.close();
}


}; // namespace Output

#ifdef USE_BOOST_PYTHON
// {

namespace BPy = boost::python;

void PyWriteVTUEquilib (FEM::Geom const & G, BPy::str const & FileName) { Output::VTU o; o.Equilib(G, BPy::extract<char const *>(FileName)()); }
void PyWriteVTK        (FEM::Geom const & G, BPy::str const & FileName) { Output::VTK             (G, BPy::extract<char const *>(FileName)()); }

// }
#endif // USE_BOOST_PYTHON

#endif // MECHSYS_OUTPUT_H
