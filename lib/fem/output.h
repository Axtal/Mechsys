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

// Boost::Python
#ifdef USE_BOOST_PYTHON
  #include <boost/python.hpp> // this includes everything
  namespace BPy = boost::python;
#endif

// MechSys
#include "fem/node.h"
#include "fem/element.h"
#include "fem/data.h"
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
	Output (FEM::Data const * D, char const * FileKey) { _initialize (D, FileKey); }

	// Destructor
	~Output ();

	// Methods
	void   Write ();                                  ///< Write output file for the current timestep
	double Val   (int iNode, char const * Key) const; ///< Output value at a specific node. This method must be called only AFTER Write, since the values at nodes must be extrapolated from the values at integration points

#ifdef USE_BOOST_PYTHON
// {
           Output      (FEM::Data const & D, BPy::str const & FileKey) { _initialize (&D, BPy::extract<char const *>(FileKey)()); }
	double PyVal       (int iNode, BPy::str const & Key) const         { return Val (iNode, BPy::extract<char const *>(Key)()); }
	void   PyGetLabels (BPy::list & Labels) const                      { for (std::map<String,int>::const_iterator iter=_map.begin(); iter!=_map.end(); iter++) Labels.append (iter->first.CStr()); }
// }
#endif // USE_BOOST_PYTHON

private:
	// Data
	FEM::Data           const * _data;     ///< The FEM data structure
	size_t                      _nimax;    ///< number of integers in a line
	size_t                      _nfmax;    ///< number of floats in a line
	size_t                      _nn;       ///< Number of Nodes
	size_t                      _ne;       ///< Number of Elements
	Util::NumStream             _nsflo;    ///< number format for floats
	Array<FEM::Element const *> _aes;      ///< Active elements
	std::map<String, int>       _map;      ///< Map to associate labels with indexes==column
	LinAlg::Matrix<double>      _vals;     ///< Matrix for nodal values collected from elements. Rows==iNode, Cols==iLabel
	String                      _file_key; ///< File key for collection
	std::ofstream             * _pvd_file; ///< File for PVD (PavaView) collection
	int                         _idx_file; ///< Increment to add to a file when working with collections

	// Private methods
	void _initialize              (FEM::Data const * D, char const * FileKey);
	void _vtu                     (char const * FileName);
	void _fill_map                ();
	void _calc_nodal_vals         ();
	void _vtu_write_header        (                 std::ostringstream & oss) const;
	void _vtu_write_geometry      (                 std::ostringstream & oss) const;
	void _vtu_write_vals_at_nodes (                 std::ostringstream & oss) const;
	void _vtu_write_vals_at_elems (                 std::ostringstream & oss) const;
	void _vtu_write_bottom        (char const * FN, std::ostringstream & oss) const;

}; // class Output


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline Output::~Output()
{
	if (_pvd_file==NULL) throw new Fatal("Output::~Output: Collection file was not open in constructor");
    (*_pvd_file) << "  </Collection>\n";
    (*_pvd_file) << "</VTKFile>\n";
	_pvd_file->close();
	delete _pvd_file;
}

inline double Output::Val(int iNode, char const * Key) const
{ 
	std::map<String,int>::const_iterator iter = _map.find(Key);
	if (iter==_map.end()) throw new Fatal(_("Output::Val: Could not find key < %s > for output"), Key);
	return _vals(iNode, iter->second);
}

inline void Output::Write()
{
	if (_pvd_file==NULL) throw new Fatal("Output::VTU: Collection file could not be open");
	String buffer;
	buffer.Printf ("    <DataSet timestep=\"%d\" file=\"%s_%d.vtu\" />\n", _idx_file, _file_key.CStr(), _idx_file);
	(*_pvd_file) << buffer;
	buffer.Printf ("%s_%d.vtu", _file_key.CStr(), _idx_file);
	_vtu (buffer.CStr());
	_idx_file++;
}


/* private */

inline void Output::_initialize(FEM::Data const * D, char const * FileKey)
{
	// Variables
	_data  = D;
	_nimax = 40;
	_nfmax = 12;
	_nn    = 0;
	_ne    = 0;
	_nsflo = Util::_8s;

	// Open collection
	String filename;  filename.Printf("%s.pvd", FileKey);
	_pvd_file = new std::ofstream();
	_pvd_file->open (filename.CStr(), std::ios::out);
	(*_pvd_file) << "<?xml version=\"1.0\" ?>\n";
	(*_pvd_file) << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	(*_pvd_file) << "  <Collection>\n";
	_idx_file = 0;
	_file_key = FileKey;
}

inline void Output::_vtu(char const * FileName)
{
	// Set array of active elements
	_aes.Resize(0);
	for (size_t i=0; i<_data->NElems(); ++i)
		if (_data->Ele(i)->IsActive()) _aes.Push (_data->Ele(i)); // Only active elements are considered

	// Set data
	_nn = _data->NNodes(); // Number of Nodes
	_ne = _aes.Size();     // Number of Elements

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
		_aes[i]->CalcDepVars ();
		_aes[i]->OutNodes    (values, labels);
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
		oss << "  " << _nsflo << _data->Nod(i)->Coord(0) << " ";
		oss <<         _nsflo << _data->Nod(i)->Coord(1) << " ";
		oss <<         _nsflo << _data->Nod(i)->Coord(2);
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
		String con;  _aes[i]->VTKConnect(con);
		oss << "  " << con;
		k++;
		OUT_NEWLINE (i,k,_ne,_nimax/_aes[i]->NNodes(),oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	size_t offset = 0;
	for (size_t i=0; i<_ne; ++i)
	{
		offset += _aes[i]->NNodes();
		oss << (k==0?"  ":" ") << offset;
		k++;
		OUT_NEWLINE (i,k,_ne,_nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<_ne; ++i)
	{
		oss << (k==0?"  ":" ") << _aes[i]->VTKCellType();
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
			try { val = _aes[i]->Val(it->first.CStr()); } catch (Exception * e) { delete e; }
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

#endif // MECHSYS_OUTPUT_H
