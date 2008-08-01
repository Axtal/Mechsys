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
	Generic (double Tol=sqrt(DBL_EPSILON)) : _tol(Tol), _is_3d(false), _is_o2(false) {} ///< Tol is the tolerance to regard two vertices as coincident

	// Destructor
	virtual ~Generic () { _erase(); }

	// Methods
	void WriteVTU (char const * FileName) const;

	// Set methods
	void SetO2      (bool IsO2=true) { _is_o2=IsO2; }                               ///< (Un)set quadratic elements
	void SetNVerts  (size_t NumVerts);                                              ///< Erase old mesh and set number of vertices
	void SetNElems  (size_t NumElems);                                              ///< Set number of elements
	void SetVert2D  (int i, bool IsOnBry, double X, double Y);                      ///< Set 2D vertex
	void SetVert3D  (int i, bool IsOnBry, double X, double Y, double Z);            ///< Set 3D vertex
	void SetElem    (int i, int Tag, bool IsOnBry, int VTKCellType, size_t NVerts); ///< Set element
	void SetElemCon (int i, int j, size_t iVert);                                   ///< Set element connectivity

	// Access methods
	bool                   Is3D     () const { return _is_3d;     } ///< Is 3D mesh ?
	Array<Vertex*> const & Verts    () const { return _verts;     } ///< Access all vertices
	Array<Elem*>   const & Elems    () const { return _elems;     } ///< Access all elements
	Array<Elem*>   const & ElemsBry () const { return _elems_bry; } ///< Access all elements on boundary
	Array<Vertex*> const & VertsBry () const { return _verts_bry; } ///< Access all vertices on boundary

#ifdef USE_BOOST_PYTHON
// {
	void   PyWriteVTU  (BPy::str const & FileName) { WriteVTU (BPy::extract<char const *>(FileName)()); }
	size_t PyGetVerts  (BPy::list & Verts) const; ///< return the number of vertices
	size_t PyGetEdges  (BPy::list & Edges) const; ///< return the number of edges
	size_t PyGetElems  (BPy::dict & Elems) const; ///< return the number of elements
	void   PySetElem   (int i, int Tag, bool IsOnBry, int VTKCellType, BPy::list const & Conn, BPy::list const & EdgeTags);
// }
#endif

protected:
	// Data
	double         _tol;       ///< Tolerance to remove duplicate nodes
	bool           _is_3d;     ///< Is 3D mesh?
	bool           _is_o2;     ///< Is quadratic element?
	Array<Vertex*> _verts;     ///< Vertices
	Array<Elem*>   _elems;     ///< Elements
	Array<Elem*>   _elems_bry; ///< Elements on boundary
	Array<Vertex*> _verts_bry; ///< Vertices on boundary

	// Methods that MAY be overloaded (otherwise, these work for linear elements such as Rods VTK_LINE)
	virtual void _vtk_con          (Elem const * E, String & Connect) const; ///< Returns a string with the connectivites (global vertices IDs) of an element
	virtual void _erase            ();                                       ///< Erase current mesh (deallocate memory)
	virtual int  _edge_to_lef_vert (int EdgeLocalID) const { return 0; }     ///< Returns the local left vertex ID for a given Local Edge ID
	virtual int  _edge_to_rig_vert (int EdgeLocalID) const { return 1; }     ///< Returns the local right vertex ID for a given Local Edge ID

private:
	// Private methods
	int _nedges (int VTKCellType) const;

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

inline void Generic::SetNVerts(size_t NumVerts)
{
	_erase            ();
	_verts.Resize     (NumVerts);
	_verts.SetValues  (NULL);
	_verts_bry.Resize (0);
}

inline void Generic::SetNElems(size_t NumElems)
{
	for (size_t i=0; i<_elems.Size(); ++i) if (_elems[i]!=NULL) delete _elems[i];
	_elems.Resize     (NumElems);
	_elems.SetValues  (NULL);
	_elems_bry.Resize (0);
}

inline void Generic::SetVert2D(int i, bool IsOnBry, double X, double Y)
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

inline void Generic::SetVert3D(int i, bool IsOnBry, double X, double Y, double Z)
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

inline void Generic::SetElem(int i, int Tag, bool IsOnBry, int VTKCellType, size_t NVerts)
{
	// Set _elems
	if (_elems[i]==NULL) _elems[i] = new Elem;
	_elems[i]->MyID        = i;
	_elems[i]->Tag         = Tag;
	_elems[i]->OnBry       = IsOnBry;
	_elems[i]->VTKCellType = VTKCellType;
	_elems[i]->V.Resize (NVerts);

	// Set _elems_bry
	if (IsOnBry) _elems_bry.Push (_elems[i]);
}

inline void Generic::SetElemCon(int i, int j, size_t iVert)
{
	_elems[i]->V[j] = _verts[iVert];
}


#ifdef USE_BOOST_PYTHON
// {

inline size_t Generic::PyGetVerts(BPy::list & Verts) const
{
	/* Out:
	 *      Verts = [(x,y,z), (x,y,z), ... num verts]
	 */
	if (Is3D())
	{
		for (size_t i=0; i<_verts.Size(); ++i)
			Verts.append (BPy::make_tuple(_verts[i]->C(0), _verts[i]->C(1), _verts[i]->C(2)));
	}
	else
	{
		for (size_t i=0; i<_verts.Size(); ++i)
			Verts.append (BPy::make_tuple(_verts[i]->C(0), _verts[i]->C(1), 0.0));
	}
	return len(Verts);
}

inline size_t Generic::PyGetEdges(BPy::list & Edges) const
{
	/* Out:
	 *      E = [[v1,v2], [v1,v2], ... num edges]
	 */
	for (size_t i=0; i<_elems.Size(); ++i)
	{
		for (int j=0; j<_nedges(_elems[i]->VTKCellType); ++j)
		{
			BPy::list pair;
			pair.append  (_elems[i]->V[_edge_to_lef_vert(j)]->MyID);
			pair.append  (_elems[i]->V[_edge_to_rig_vert(j)]->MyID);
			Edges.append (pair);
		}
	}
	return len(Edges);
}

inline size_t Generic::PyGetElems(BPy::dict & Elems) const
{
	/* Out:
	 *    
	 *    Elems = {
	 *      'tags'    : [t1,t2,t3, ... num elements]
	 *      'onbs'    : [ 1, 0, 1, ... num elements]
	 *      'vtks'    : [ 9, 9,12, ... num elements]
	 *      'cons'    : {id0:[node0, node1, node3, num nodes in element 0],
	 *                   id1:[n0,n1,n2,n3,n4,n5,n6,n7], ... num elems}
	 *      'etags'   : {id0:[tag_edge_0, tag_edge_1, tag_edge_2, tag_edge_3],
	 *                   id1:[t0,t1,t2,t3], ... num elems]}
	 *      'ftags'   : {id0:[tag_face_0, tag_face_1, tag_face_2, tag_face_3, tag_face_4, tag_face_5],
	 *                   id1:[t0,t1,t2,t3,t4,t5], ... num elems]}
	 *      'etags_g' : {(L,R):tag, (L,R):tag, ... num edge tags]}
	 *    }
	 */
	BPy::list tags,  onbs,  vtks;
	BPy::dict cons, etags, ftags, etags_g;
	for (size_t i=0; i<_elems.Size(); ++i)
	{
		// Data
		tags.append (_elems[i]->Tag);
		onbs.append (_elems[i]->OnBry ? 1 : 0);
		vtks.append (_elems[i]->VTKCellType);

		// Connectivities
		BPy::list co;
		for (size_t j=0; j<_elems[i]->V.Size(); ++j)
			co.append (_elems[i]->V[j]->MyID);
		cons[i] = co;

		// ETags and Global ETags
		BPy::list et;
		for (int j=0; j<_elems[i]->ETags.Size(); ++j) // j is the local edge id
		{
			int tag = _elems[i]->ETags(j);
			et.append (tag);
			if (tag<0)
			{
				int L = _elems[i]->V[_edge_to_lef_vert(j)]->MyID;
				int R = _elems[i]->V[_edge_to_rig_vert(j)]->MyID;
				etags_g[BPy::make_tuple(L, R)] = tag;
			}
		}
		etags[i] = et;

		// FTags
		BPy::list ft;
		for (int j=0; j<_elems[i]->FTags.Size(); ++j) // j is the local face id
		{
			int tag = _elems[i]->FTags(j);
			ft.append (tag);
		}
		ftags[i] = ft;
	}
	Elems["tags" ]   = tags;
	Elems["onbs" ]   = onbs;
	Elems["vtks" ]   = vtks;
	Elems["cons" ]   = cons;
	Elems["etags"]   = etags;
	Elems["etags_g"] = etags_g;
	Elems["ftags"]   = ftags;
	return _elems.Size();
}

inline void Generic::PySetElem(int i, int Tag, bool IsOnBry, int VTKCellType, BPy::list const & Conn, BPy::list const & EdgeTags)
{
	// Set Elements & Connectivity
	int nverts = len(Conn);
	SetElem (i, Tag, IsOnBry, VTKCellType, nverts);
	for (int j=0; j<nverts; ++j)
		SetElemCon (i, j, BPy::extract<int>(Conn[j])());

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


/* protected */

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


/* private */

inline int Generic::_nedges(int VTKCellType) const
{
	switch (VTKCellType)
	{
		case VTK_LINE:                 { return  1; }
		case VTK_TRIANGLE:             { return  3; }
		case VTK_QUAD:                 { return  4; }
		case VTK_TETRA:                { return  6; }
		case VTK_HEXAHEDRON:           { return 12; }
		case VTK_QUADRATIC_EDGE:       { return  2; }
		case VTK_QUADRATIC_TRIANGLE:   { return  6; }
		case VTK_QUADRATIC_QUAD:       { return  8; }
		case VTK_QUADRATIC_TETRA:      { return 12; }
		case VTK_QUADRATIC_HEXAHEDRON: { return 24; }
		default: throw new Fatal("Generic::_nedges: VTKCellType==%d is invalid", VTKCellType);
	}
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
