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

#ifndef OUT_NEWLINE_DEFINED
  #define OUT_NEWLINE_DEFINED
  #define OUT_NEWLINE(I,K,N,KMAX,OF) if (K>KMAX) { OF<<(I<N-1?"\n        ":"\n"); K=0; } else if (I==N-1) { OF<<"\n"; }
#endif

class Output
{

public:
	// Constructor
	Output () : _nimax(40), _nfmax(12), _nn(0), _ne(0), _nsflo(Util::_8s) {}

	// Methods
	void VTK   (FEM::Geom const * G, char const * FileName); ///< Write a ParaView-VTK file
	void VTU   (FEM::Geom const * G, char const * FileName); ///< Write a ParaView-VTU file
	void VTUcg (FEM::Geom const * G, char const * FileName); ///< Write a ParaView-VTU file with element values at CG

private:
	// Data
	size_t                      _nimax; ///< number of integers in a line
	size_t                      _nfmax; ///< number of floats in a line
	size_t                      _nn;    ///< Number of Nodes
	size_t                      _ne;    ///< Number of Elements
	Util::NumStream             _nsflo; ///< number format for floats
	FEM::Geom           const * _g;     ///< Geometry structure
	Array<FEM::Element const *> _aes;   ///< Active elements
	std::map<String, int>       _map;   ///< Map to associate labels with indexes
	LinAlg::Matrix<double>      _vals;  ///< Matrix for nodal values collected from elements

	// Methods
	void _fill_map                ();
	void _calc_nodal_vals         ();
	void _vtk_write_header        (                 std::ostringstream & oss) const;
	void _vtk_write_geometry      (                 std::ostringstream & oss) const;
	void _vtk_write_vals_at_nodes (                 std::ostringstream & oss) const;
	void _vtk_write_bottom        (char const * FN, std::ostringstream & oss) const;
	void _vtu_write_header        (                 std::ostringstream & oss) const;
	void _vtu_write_geometry      (                 std::ostringstream & oss) const;
	void _vtu_write_vals_at_nodes (                 std::ostringstream & oss) const;
	void _vtu_write_vals_at_elems (                 std::ostringstream & oss) const;
	void _vtu_write_bottom        (char const * FN, std::ostringstream & oss) const;

}; // class Output


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void Output::VTK(FEM::Geom const * G, char const * FileName)
{
	// Set array of active elements
	_aes.Resize(0);
	for (size_t i=0; i<G->NElems(); ++i)
		if (G->Ele(i)->IsActive()) _aes.Push (G->Ele(i)); // Only active elements are considered

	// Set data
	_nn = G->NNodes(); // Number of Nodes
	_ne = _aes.Size(); // Number of Elements
	_g  = G;           // geometry

	std::ostringstream oss;    // Output structure
	_fill_map           ();    // fill map
	_calc_nodal_vals    ();    // Calculate averaged nodal values
	_vtk_write_header   (oss); // Header
	_vtk_write_geometry (oss); // Geometry

	// Data -- nodes
	oss << "POINT_DATA " << _nn << std::endl;

	// Data -- nodes -- displacements and forces
	if (_map.find("ux")!=_map.end())
	{
		oss << "VECTORS " << "displ float" << std::endl;
		for (size_t j=0; j<_nn; ++j)
			oss << _nsflo<< _vals(j, _map["ux"]) << _nsflo<< _vals(j, _map["uy"]) << _nsflo<< ((_map.find("uz")==_map.end()) ? 0.0 : _vals(j, _map["uz"])) << "\n";
		oss << "\n";
		oss << "VECTORS " << "force float" << std::endl;
		for (size_t j=0; j<_nn; ++j)
			oss << _nsflo<< _vals(j, _map["fx"]) << _nsflo<< _vals(j, _map["fy"]) << _nsflo<< ((_map.find("fz")==_map.end()) ? 0.0 : _vals(j, _map["uz"])) << "\n";
		oss << "\n";
	}

	// Data -- nodes -- scalars
	_vtk_write_vals_at_nodes (oss);

	// Bottom and write to file
	_vtk_write_bottom (FileName, oss);
}

inline void Output::VTU(FEM::Geom const * G, char const * FileName)
{
	// Set array of active elements
	_aes.Resize(0);
	for (size_t i=0; i<G->NElems(); ++i)
		if (G->Ele(i)->IsActive()) _aes.Push (G->Ele(i)); // Only active elements are considered

	// Set data
	_nn = G->NNodes(); // Number of Nodes
	_ne = _aes.Size(); // Number of Elements
	_g  = G;           // geometry

	std::ostringstream oss;    // Output structure
	_fill_map           ();    // fill map
	_calc_nodal_vals    ();    // Calculate averaged nodal values
	_vtu_write_header   (oss); // Header
	_vtu_write_geometry (oss); // Geometry

	// Data -- nodes -- displacements and forces
	oss << "      <PointData Scalars=\"TheScalars\" Vectors=\"TheVectors\">\n";
	if (_map.find("ux")!=_map.end())
	{
		oss << "        <DataArray type=\"Float32\" Name=\"" << "displ" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
		size_t k = 0; oss << "        ";
		for (size_t i=0; i<_nn; ++i)
		{
			oss << " " << _nsflo << _vals(i, _map["ux"]) << " " << _nsflo<< _vals(i, _map["uy"]) << " " << _nsflo<< ((_map.find("uz")==_map.end()) ? 0.0 : _vals(i, _map["uz"])) << " ";
			k++;
			OUT_NEWLINE (i,k,_nn,_nfmax/3,oss);
		}
		oss << "        </DataArray>\n";
		oss << "        <DataArray type=\"Float32\" Name=\"" << "force" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
		k = 0; oss << "        ";
		for (size_t i=0; i<_nn; ++i)
		{
			oss << " " << _nsflo << _vals(i, _map["fx"]) << " " << _nsflo<< _vals(i, _map["fy"]) << " " << _nsflo<< ((_map.find("fz")==_map.end()) ? 0.0 : _vals(i, _map["fz"])) << " ";
			k++;
			OUT_NEWLINE (i,k,_nn,_nfmax/3,oss);
		}
		oss << "        </DataArray>\n";
	}

	// Data -- nodes -- scalars
	_vtu_write_vals_at_nodes (oss);
	oss << "      </PointData>\n";

	// Bottom and write to file
	_vtu_write_bottom (FileName, oss);
}

inline void Output::VTUcg(FEM::Geom const * G, char const * FileName)
{
	// Set array of active elements
	_aes.Resize(0);
	for (size_t i=0; i<G->NElems(); ++i)
		if (G->Ele(i)->IsActive()) _aes.Push (G->Ele(i)); // Only active elements are considered

	// Set data
	_nn = G->NNodes(); // Number of Nodes
	_ne = _aes.Size(); // Number of Elements
	_g  = G;           // geometry

	std::ostringstream oss;    // Output structure
	_fill_map           ();    // fill map
	_vtu_write_header   (oss); // Header
	_vtu_write_geometry (oss); // Geometry

	// Data -- nodes -- temperature and heat
	oss << "      <PointData Scalars=\"TheScalars\" Vectors=\"TheVectors\">\n";
	if (_map.find("T")!=_map.end())
	{
		oss << "        <DataArray type=\"Float32\" Name=\"" << "T" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
		size_t k = 0; oss << "        ";
		for (size_t i=0; i<_nn; ++i)
		{
			oss << (k==0?"  ":" ") << _g->Nod(i)->Val("T");
			k++;
			OUT_NEWLINE (i,k,_nn,_nfmax,oss);
		}
		oss << "        </DataArray>\n";
		oss << "        <DataArray type=\"Float32\" Name=\"" << "F" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
		k = 0; oss << "        ";
		for (size_t i=0; i<_nn; ++i)
		{
			oss << (k==0?"  ":" ") << _g->Nod(i)->Val("F");
			k++;
			OUT_NEWLINE (i,k,_nn,_nfmax,oss);
		}
		oss << "        </DataArray>\n";
	}

	// Data -- nodes -- displacements and forces
	if (_map.find("ux")!=_map.end())
	{
		oss << "        <DataArray type=\"Float32\" Name=\"" << "displ" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
		size_t k = 0; oss << "        ";
		for (size_t i=0; i<_nn; ++i)
		{
			oss << " " << _nsflo << (_g->Nod(i)->HasVar("ux") ? _g->Nod(i)->Val("ux"): 0.0) << " ";
			oss <<        _nsflo << (_g->Nod(i)->HasVar("uy") ? _g->Nod(i)->Val("uy"): 0.0) << " ";
			oss <<        _nsflo << (_g->Nod(i)->HasVar("uz") ? _g->Nod(i)->Val("uz"): 0.0) << " ";
			k++;
			OUT_NEWLINE (i,k,_nn,_nfmax/3,oss);
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
			OUT_NEWLINE (i,k,_nn,_nfmax/3,oss);
		}
		oss << "        </DataArray>\n";
	}
	oss << "      </PointData>\n";

	// Data -- elements
	_vtu_write_vals_at_elems (oss);

	// Bottom and write to file
	_vtu_write_bottom (FileName, oss);
}


/* private */

inline void Output::_fill_map()
{
	// Get all possible labels from elements
	_map.clear();
	for (size_t i=0; i<_ne; ++i)
	{
		Array<String> labels;
		_aes[i]->GetLabels (labels);
		for (size_t j=0; j<labels.Size(); ++j)
		{
			if (_map.find(labels[j])==_map.end())
			{
				size_t icomp = _map.size();
				_map[labels[j]] = icomp; // _map.size()-1 trick was not working on Mac OS X. Now it's ok.
			}
		}
	}
}

inline void Output::_calc_nodal_vals()
{
	// Collect nodal values
	size_t ncomps = _map.size();
	LinAlg::Matrix<size_t> refs;
	_vals.Resize (_nn, ncomps);  _vals.SetValues (0.0);
	 refs.Resize (_nn, ncomps);   refs.SetValues (  0);
	for (size_t i=0; i<_ne; ++i)
	{
		LinAlg::Matrix<double> values;
		Array<String>          labels;
		_aes[i]->OutNodes (values, labels);
		size_t elem_nnodes = _aes[i]->NNodes();
		for (size_t j=0; j<labels.Size(); ++j)
		{
			int index = _map[labels[j]];
			for (size_t k=0; k<elem_nnodes; ++k)
			{
				int inode           = _aes[i]->Nod(k)->GetID();
				_vals(inode,index) += values(k, j); // accumulate values
				 refs(inode,index)++;               // update references number
			}
		}
	}
	
	// Compute averaged values
	for (size_t i=0; i<_nn;    i++)
	for (size_t j=0; j<ncomps; j++)
	{
		if (refs(i,j)!=0) _vals(i,j) /= refs(i,j);
		else              _vals(i,j)  = 0.0;
	}
}

inline void Output::_vtk_write_header(std::ostringstream & oss) const
{
	// Write Legacy VTK file header
	oss << "# vtk DataFile Version 3.0"     << std::endl;
	oss << "MechSys/FEM"                    << std::endl;
	oss << "ASCII"                          << std::endl;
	oss << "DATASET UNSTRUCTURED_GRID"      << std::endl;
	oss <<                                     std::endl;
}

inline void Output::_vtk_write_geometry(std::ostringstream & oss) const
{
	// Node coordinates
	oss << "POINTS " << _nn << " float" << std::endl;
	for (size_t i=0; i<_nn; ++i)
		oss << _nsflo << _g->Nod(i)->X() << _nsflo << _g->Nod(i)->Y() << _nsflo << _g->Nod(i)->Z() << std::endl;
	oss << std::endl;

	// Total number of CELLS data = Sum (1 + NNodes); 1 for the numPts label
	int n_data = 0;
	for (size_t i=0; i<_ne; ++i) n_data += 1 + _aes[i]->NNodes();

	// Elements connectivities
	oss << "CELLS "<< _ne << " " << n_data << std::endl;
	for (size_t i=0; i<_ne; ++i)
	{
		String connect; _aes[i]->VTKConnect(connect);
		oss << _aes[i]->NNodes() << " " << connect << std::endl;
	}
	oss << std::endl;

	// Cell types
	oss << "CELL_TYPES " << _ne << std::endl;
	for (size_t i=0; i<_ne; ++i)
		oss << _aes[i]->VTKCellType() << std::endl;
	oss << std::endl;
}

inline void Output::_vtk_write_vals_at_nodes(std::ostringstream & oss) const
{
	std::map<String,int>::const_iterator it;
	for (it=_map.begin(); it!=_map.end(); it++)
	{
		oss << "SCALARS " << it->first << " float 1" << std::endl;
		oss << "LOOKUP_TABLE default" << std::endl;
		for (size_t i=0; i<_nn; ++i)
			oss << _nsflo << _vals(i, it->second) << std::endl;
		oss << std::endl;
	}
}

inline void Output::_vtk_write_bottom(char const * FN, std::ostringstream & oss) const
{
	std::ofstream of(FN, std::ios::out);
	of << oss.str();
	of.close();
}

inline void Output::_vtu_write_header(std::ostringstream & oss) const
{
	// Header
	oss << "<?xml version=\"1.0\"?>\n";
	oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	oss << "  <UnstructuredGrid>\n";
	oss << "    <Piece NumberOfPoints=\"" << _nn << "\" NumberOfCells=\"" << _ne << "\">\n";
}

inline void Output::_vtu_write_geometry(std::ostringstream & oss) const
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
		OUT_NEWLINE (i,k,_nn,_nfmax/3,oss);
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
		OUT_NEWLINE (i,k,_ne,_nimax/_g->Ele(i)->NNodes(),oss);
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
		OUT_NEWLINE (i,k,_ne,_nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<_ne; ++i)
	{
		oss << (k==0?"  ":" ") << _g->Ele(i)->VTKCellType();
		k++;
		OUT_NEWLINE (i,k,_ne,_nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </Cells>\n";
}

inline void Output::_vtu_write_vals_at_nodes(std::ostringstream & oss) const
{
	std::map<String,int>::const_iterator it;
	for (it=_map.begin(); it!=_map.end(); it++)
	{
		oss << "        <DataArray type=\"Float32\" Name=\"" << it->first << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
		size_t k = 0; oss << "        ";
		for (size_t i=0; i<_nn; ++i)
		{
			oss << (k==0?"  ":" ") << _vals(i, it->second);
			k++;
			OUT_NEWLINE (i,k,_nn,_nfmax,oss);
		}
		oss << "        </DataArray>\n";
	}
}

inline void Output::_vtu_write_vals_at_elems(std::ostringstream & oss) const
{
	oss << "      <CellData Scalars=\"TheScalars\">\n";
	std::map<String,int>::const_iterator it;
	for (it=_map.begin(); it!=_map.end(); it++)
	{
		oss << "        <DataArray type=\"Float32\" Name=\""<< it->first <<"\" NumberOfComponents=\"1\" format=\"ascii\">\n";
		size_t k = 0; oss << "        ";
		for (size_t i=0; i<_ne; ++i)
		{
			double val = 0.0;
			try { val = _g->Ele(i)->Val(it->first.CStr()); } catch (Exception * e) { delete e; }
			oss << (k==0?"  ":" ") << val;
			k++;
			OUT_NEWLINE (i,k,_ne,_nfmax,oss);
		}
		oss << "        </DataArray>\n";
	}
	oss << "      </CellData>\n";
}

inline void Output::_vtu_write_bottom(char const * FN, std::ostringstream & oss) const
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


#ifdef USE_BOOST_PYTHON
// {

namespace BPy = boost::python;

void PyOutputVTK   (FEM::Geom const & G, BPy::str const & FileName) { Output o; o.VTK  (&G, BPy::extract<char const *>(FileName)()); }
void PyOutputVTU   (FEM::Geom const & G, BPy::str const & FileName) { Output o; o.VTU  (&G, BPy::extract<char const *>(FileName)()); }
void PyOutputVTUcg (FEM::Geom const & G, BPy::str const & FileName) { Output o; o.VTUcg(&G, BPy::extract<char const *>(FileName)()); }

// }
#endif // USE_BOOST_PYTHON

#endif // MECHSYS_OUTPUT_H
