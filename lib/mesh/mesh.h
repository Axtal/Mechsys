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

#ifndef MECHSYS_MESH_H
#define MECHSYS_MESH_H

// STL
#include <iostream>
#include <fstream>
#include <cfloat>   // for DBL_EPSILON

// Blitz++
#include <blitz/tinyvec-et.h>

// Boost::Python
#ifdef USE_BOOST_PYTHON
  #include <boost/python.hpp> // this includes everything
  namespace BPy = boost::python;
#endif

// MechSys
#include "util/array.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "fem/elems/vtkCellType.h"

#ifndef VTU_NEWLINE_DEFINED
  #define VTU_NEWLINE_DEFINED
  #define VTU_NEWLINE(I,K,N,KMAX,OF) if (K>KMAX) { OF<<(I<N-1?"\n        ":"\n"); K=0; } else if (I==N-1) { OF<<"\n"; }
#endif

using LinAlg::Vector;
using LinAlg::Matrix;
using blitz::TinyVector;

namespace Mesh
{

struct Edge
{
	int L; // Left vertex id
	int R; // Right vertex id
};

struct Elem;

struct Share
{
	Elem* E; ///< The element
	int   N; ///< Local node index. Example: 2D=>0,1,2,3, 3D=>0,1,2,3,4,5,6,7
};

struct Vertex
{
	long              MyID;    ///< ID
	bool              OnBry;   ///< Is on boundary?
	TinyVector<int,3> EdgesID; ///< Local indexes (3) of what edges this vertex is located on (from 0 to 12). -1 => Not on boundary
	TinyVector<int,3> FacesID; ///< Local indexes (3) of what faces this vertex is located on (from 0 to 6) . -1 => Not on boundary
	bool              Dupl;    ///< Is this a duplicated node?
	Vector<double>    C;       ///< X, Y, and Z coordinates
	Array<Share>      Shares;  ///< Shared elements
};

struct Elem
{
	long           MyID;        ///< ID
	int            Tag;         ///< Element tag. Required for setting up of attributes, for example.
	bool           OnBry;       ///< On boundary?
	int            VTKCellType; ///< VTK cell type such as VTK_LINE, VTK_TRIANGLE, VTK_QUAD, VTK_HEXAHEDRON, etc.
	Array<Vertex*> V;           ///< Connectivity
	Vector<int>    ETags;       ///< Edge tags (size==nLocalEdges)
	Vector<int>    FTags;       ///< Face tags (size==nLocalFaces)
};

class Generic
{
public:
	// Constructor
	Generic (double Tol=sqrt(DBL_EPSILON)) : _tol(Tol), _is_3d(false) {} ///< Tol is the tolerance to regard two vertices as coincident

	// Destructor
	virtual ~Generic () { _erase(); }

	// Methods
	void WriteVTU (char const * FileName) const;

	// Access methods
	bool                   Is3D     () const { return _is_3d;     } ///< Is 3D mesh ?
	Array<Vertex*> const & Verts    () const { return _verts;     } ///< Access all vertices
	Array<Elem*>   const & Elems    () const { return _elems;     } ///< Access all elements
	Array<Elem*>   const & ElemsBry () const { return _elems_bry; } ///< Access all elements on boundary
	Array<Vertex*> const & VertsBry () const { return _verts_bry; } ///< Access all vertices on boundary

#ifdef USE_BOOST_PYTHON
// {
	void PyWriteVTU   (BPy::str const & FileName) { WriteVTU (BPy::extract<char const *>(FileName)()); }
	void PyGetVerts   (BPy::list & V) const;
	void PyGetElems   (BPy::list & E) const;
	void PyGetEleTags (BPy::list & T) const;
	void PyGetETags   (BPy::list & T) const;
	void PySetNVerts  (size_t NumVerts);
	void PySetNElems  (size_t NumElems);
	void PySetVert2D  (int i, bool IsOnBry, double X, double Y);
	void PySetVert3D  (int i, bool IsOnBry, double X, double Y, double Z);
	void PySetElem    (int i, int Tag, bool IsOnBry, int VTKCellType, BPy::list const & Conn, BPy::list const & EdgeTags);
// }
#endif

protected:
	// Data
	double         _tol;       ///< Tolerance to remove duplicate nodes
	bool           _is_3d;     ///< Is 3D mesh?
	Array<Vertex*> _verts;     ///< Vertices
	Array<Elem*>   _elems;     ///< Elements
	Array<Elem*>   _elems_bry; ///< Elements on boundary
	Array<Vertex*> _verts_bry; ///< Vertices on boundary

	// Methods that MAY be overloaded (otherwise, these work for linear elements such as Rods VTK_LINE)
	virtual void _vtk_con          (Elem const * E, String & Connect) const; ///< Returns a string with the connectivites (global vertices IDs) of an element
	virtual void _erase            ();                                       ///< Erase current mesh (deallocate memory)
	virtual int  _edge_to_lef_vert (int EdgeLocalID) const { return 0; }     ///< Returns the local left vertex ID for a given Local Edge ID
	virtual int  _edge_to_rig_vert (int EdgeLocalID) const { return 1; }     ///< Returns the local right vertex ID for a given Local Edge ID

}; // class Generic


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void Generic::WriteVTU(char const * FileName) const
{
	// Results
	std::ostringstream oss;

	// Data
	size_t nn = _verts.Size(); // Number of Nodes
	size_t ne = _elems.Size(); // Number of Elements

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
		oss << "  " << nsflo <<         _verts[i]->C(0) << " ";
		oss <<         nsflo <<         _verts[i]->C(1) << " ";
		oss <<         nsflo << (_is_3d?_verts[i]->C(2):0.0);
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
		String con; _vtk_con (_elems[i], con);
		oss << "  " << con;
		k++;
		VTU_NEWLINE (i,k,ne,nimax/_elems[i]->V.Size(),oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	size_t offset = 0;
	for (size_t i=0; i<ne; ++i)
	{
		offset += _elems[i]->V.Size();
		oss << (k==0?"  ":" ") << offset;
		k++;
		VTU_NEWLINE (i,k,ne,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		oss << (k==0?"  ":" ") << _elems[i]->VTKCellType;
		k++;
		VTU_NEWLINE (i,k,ne,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </Cells>\n";

	// Data -- nodes
	oss << "      <PointData Scalars=\"TheScalars\">\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "onbry" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<nn; ++i)
	{
		oss << (k==0?"  ":" ") << _verts[i]->OnBry;
		k++;
		VTU_NEWLINE (i,k,nn,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "local_edge_id" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<nn; ++i)
	{
		oss << (k==0?"  ":" ") << _verts[i]->EdgesID(0) << " " << _verts[i]->EdgesID(1) << " " << _verts[i]->EdgesID(2) << " ";
		k++;
		VTU_NEWLINE (i,k,nn,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "local_face_id" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<nn; ++i)
	{
		oss << (k==0?"  ":" ") << _verts[i]->FacesID(0) << " " << _verts[i]->FacesID(1) << " " << _verts[i]->FacesID(2) << " ";
		k++;
		VTU_NEWLINE (i,k,nn,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "shares" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<nn; ++i)
	{
		oss << (k==0?"  ":" ") << _verts[i]->Shares.Size();
		k++;
		VTU_NEWLINE (i,k,nn,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </PointData>\n";

	// Data -- elements
	oss << "      <CellData Scalars=\"TheScalars\">\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "onbry" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		oss << (k==0?"  ":" ") << _elems[i]->OnBry;
		k++;
		VTU_NEWLINE (i,k,ne,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "tag" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		oss << (k==0?"  ":" ") << _elems[i]->Tag;
		k++;
		VTU_NEWLINE (i,k,ne,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </CellData>\n";

	// Bottom
	oss << "    </Piece>\n";
	oss << "  </UnstructuredGrid>\n";
	oss << "</VTKFile>" << std::endl;

	// Write to file
	std::ofstream of(FileName, std::ios::out);
	of << oss.str();
	of.close();
}

#ifdef USE_BOOST_PYTHON
// {

inline void Generic::PyGetVerts(BPy::list & V) const
{
	if (Is3D())
	{
		for (size_t i=0; i<_verts.Size(); ++i)
			V.append (BPy::make_tuple(_verts[i]->C(0), _verts[i]->C(1), _verts[i]->C(2)));
	}
	else
	{
		for (size_t i=0; i<_verts.Size(); ++i)
			V.append (BPy::make_tuple(_verts[i]->C(0), _verts[i]->C(1), 0.0));
	}
}

inline void Generic::PyGetElems(BPy::list & E) const
{
	for (size_t i=0; i<_elems.Size(); ++i)
	{
		BPy::list conn;
		for (size_t j=0; j<_elems[i]->V.Size(); ++j)
			conn.append (_elems[i]->V[j]->MyID);
		E.append (conn);
	}
}

inline void Generic::PyGetEleTags(BPy::list & Tags) const
{
	/* Returns a list of tuples: [(int,int, (int,int), ..., num of elems with tags]
	 *
	 *   Each tuple has two values:  element ID, element tag
	 */
	for (size_t i=0; i<_elems.Size(); ++i)
		Tags.append (BPy::make_tuple(_elems[i]->MyID, _elems[i]->Tag));
}

inline void Generic::PyGetETags(BPy::list & Tags) const
{
	/* Returns a list of tuples: [(int,int,int,int), (int,int,int,int), ..., num of elems with tags]
	 *
	 *   Each tuple has three values: (eid, L, R, tag)
	 *
	 *   where:  eid: element ID
	 *           L:   global ID of the left vertex on edge
	 *           R:   global ID of the right vertex on edge
	 *           tag: edge tag
	 */
	for (size_t i=0; i<_elems_bry.Size(); ++i) // elements on boundary
	{
		for (int j=0; j<_elems_bry[i]->ETags.Size(); ++j) // j is the local_edge_id
		{
			int tag = _elems_bry[i]->ETags(j);
			if (tag<0)
			{
				int eid = _elems_bry[i]->MyID;
				int L   = _elems_bry[i]->V[_edge_to_lef_vert(j)]->MyID;
				int R   = _elems_bry[i]->V[_edge_to_rig_vert(j)]->MyID;
				Tags.append (BPy::make_tuple(eid, L, R, tag));
			}
		}
	}
}

inline void Generic::PySetNVerts(size_t NumVerts)
{
	_erase();
	_verts.Resize(NumVerts);
	_verts.SetValues(NULL);
	_verts_bry.Resize(0);
}

inline void Generic::PySetNElems(size_t NumElems)
{
	for (size_t i=0; i<_elems.Size(); ++i) if (_elems[i]!=NULL) delete _elems[i];
	_elems.Resize(NumElems);
	_elems.SetValues(NULL);
	_elems_bry.Resize(0);
}

inline void Generic::PySetVert2D(int i, bool IsOnBry, double X, double Y)
{
	// Set _verts
	if (_verts[i]==NULL) _verts[i] = new Vertex;
	_verts[i]->MyID    = i;
	_verts[i]->OnBry   = IsOnBry;
	_verts[i]->EdgesID = -1;
	_verts[i]->FacesID = -1;
	_verts[i]->Dupl    = false;
	_verts[i]->C.Resize(2);
	_verts[i]->C = X, Y;

	// Set _verts_bry
	if (IsOnBry) _verts_bry.Push (_verts[i]);
}

inline void Generic::PySetVert3D(int i, bool IsOnBry, double X, double Y, double Z)
{
	// Set _verts
	if (_verts[i]==NULL) _verts[i] = new Vertex;
	_verts[i]->MyID    = i;
	_verts[i]->OnBry   = IsOnBry;
	_verts[i]->EdgesID = -1;
	_verts[i]->FacesID = -1;
	_verts[i]->Dupl    = false;
	_verts[i]->C.Resize(3);
	_verts[i]->C = X, Y, Z;

	// Set _verts_bry
	if (IsOnBry) _verts_bry.Push (_verts[i]);
}

inline void Generic::PySetElem(int i, int Tag, bool IsOnBry, int VTKCellType, BPy::list const & Conn, BPy::list const & EdgeTags)
{
	// Set _elems
	if (_elems[i]==NULL) _elems[i] = new Elem;
	_elems[i]->MyID        = i;
	_elems[i]->Tag         = Tag;
	_elems[i]->OnBry       = IsOnBry;
	_elems[i]->VTKCellType = VTKCellType;

	// Set connectivity
	int nverts = len(Conn);
	_elems[i]->V.Resize(nverts);
	for (int j=0; j<nverts; ++j)
		_elems[i]->V[j] = _verts[BPy::extract<int>(Conn[j])()];

	// Set _elems_bry
	if (IsOnBry) _elems_bry.Push(_elems[i]);

	// Set edge tags
	int netags = len(EdgeTags);
	if (netags>0)
	{
		_elems[i]->ETags.Resize(netags);
		for (int j=0; j<netags; ++j)
			_elems[i]->ETags(j) = BPy::extract<int>(EdgeTags[j])();
	}
}

// }
#endif // USE_BOOST_PYTHON


/* private */

inline void Generic::_vtk_con(Elem const * E, String & Connect) const
{
	Connect.Printf("%d %d",E->V[0]->MyID,
	                       E->V[1]->MyID);
}

inline void Generic::_erase()
{
	for (size_t i=0; i<_verts.Size(); ++i) if (_verts[i]!=NULL) delete _verts[i]; // it is only necessary to delete nodes in _verts array
	for (size_t i=0; i<_elems.Size(); ++i) if (_elems[i]!=NULL) delete _elems[i]; // it is only necessary to delete elems in _elems array
	_is_3d = false;
	_verts      .Resize(0);
	_elems      .Resize(0);
	_elems_bry  .Resize(0);
	_verts_bry  .Resize(0);
}

/* output */
std::ostream & operator<< (std::ostream & os, Mesh::Vertex const & V)
{
	os << "[" << V.MyID << "] " << V.OnBry << " " << V.Dupl;
	if (V.C.Size()==2) os << " (" << Util::_8_4<<V.C(0) << Util::_8_4<<V.C(1) << ")";
	else               os << " (" << Util::_8_4<<V.C(0) << Util::_8_4<<V.C(1) << Util::_8_4<<V.C(2) << ")";
	os << " {" << V.EdgesID(0) << "," << V.EdgesID(1) << "," << V.EdgesID(2) << "}";
	os << " {" << V.FacesID(0) << "," << V.FacesID(1) << "," << V.FacesID(2) << "}";
	return os;
}
std::ostream & operator<< (std::ostream & os, Mesh::Elem const & E)
{
	os << "[" << E.MyID << "] " << E.Tag << " " << E.OnBry;
	os << " {"; for (int i=0; i<E.ETags.Size(); ++i) os << E.ETags(i) << (i==E.ETags.Size()-1?"":","); os << "}";
	os << " {"; for (int i=0; i<E.FTags.Size(); ++i) os << E.FTags(i) << (i==E.FTags.Size()-1?"":","); os << "}";
	os << "\n";
	for (size_t i=0; i<E.V.Size(); ++i)
		if (E.V[i]!=NULL) os << "   " << (*E.V[i]) << "\n";
	return os;
}
std::ostream & operator<< (std::ostream & os, Mesh::Generic const & G)
{
	for (size_t i=0; i<G.Elems().Size(); ++i)
		if (G.Elems()[i]!=NULL) os << (*G.Elems()[i]);
	return os;
}

}; // namespace Mesh

#endif // MECHSYS_MESH_H
