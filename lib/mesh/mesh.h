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
	size_t L; // Left vertex id
	size_t R; // Right vertex id
};

struct Face
{
	size_t v0;
	size_t v1;
	size_t v2;
	size_t v3;
	size_t v4;
	size_t v5;
	size_t v6;
	size_t v7;
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
	Generic () : _is_3d(false), _is_o2(false) {}

	// Destructor
	virtual ~Generic () { _erase(); }

	// Methods
	void WriteVTU (char const * FileName) const; ///< Write output file for ParaView

	// Set methods
	virtual void Set3D       (bool Is3D=true) { _is_3d=Is3D; }                      ///< (Un)set 3D mesh
	virtual void SetO2       (bool IsO2=true) { _is_o2=IsO2; }                      ///< (Un)set quadratic elements
	virtual void SetNVerts   (size_t NumVerts);                                     ///< Erase old mesh and set number of vertices
	virtual void SetNElems   (size_t NumElems);                                     ///< Set number of elements
	virtual void SetVert     (int i, bool IsOnBry, double X, double Y, double Z=0); ///< Set vertex
	virtual void SetElem     (int i, int Tag, bool IsOnBry, int VTKCellType);       ///< Set element
	virtual void SetElemCon  (int i, int j, size_t iVert);                          ///< Set element connectivity
	virtual void SetElemETag (int i, int j, int Tag);                               ///< Set element's edge tag
	virtual void SetElemFTag (int i, int j, int Tag);                               ///< Set element's face tag

	// Get methods
	        bool   Is3D            ()                   const { return _is_3d;                  } ///< Is 3D mesh ?
	virtual size_t NVerts          ()                   const { return _verts.Size();           } ///< Return the number of vertices
	virtual size_t NVertsBry       ()                   const { return _verts_bry.Size();       } ///< Return the number of elements on boundary
	virtual size_t NElems          ()                   const { return _elems.Size();           } ///< Return the number of elements
	virtual size_t NElemsBry       ()                   const { return _elems_bry.Size();       } ///< Return the number of elements on boundary
	virtual long   VertBry         (size_t i)           const { return _verts_bry[i]->MyID;     } ///< Return the ID of a vertex on boundary
	virtual long   ElemBry         (size_t i)           const { return _elems_bry[i]->MyID;     } ///< Return the ID of an element on boundary
	virtual bool   IsVertOnBry     (size_t i)           const { return _verts[i]->OnBry;        } ///< Return whether a vertex in on boundary or not
	virtual double VertX           (size_t i)           const { return _verts[i]->C(0);         } ///< Return the X coordinate of a vertex i
	virtual double VertY           (size_t i)           const { return _verts[i]->C(1);         } ///< Return the Y coordinate of a vertex i
	virtual double VertZ           (size_t i)           const { return _verts[i]->C(2);         } ///< Return the Z coordinate of a vertex i
	virtual int    ElemTag         (size_t i)           const { return _elems[i]->Tag;          } ///< Return the Tag of a element i
	virtual bool   IsElemOnBry     (size_t i)           const { return _elems[i]->OnBry;        } ///< Return whether an element is on boundary or not
	virtual int    ElemVTKCellType (size_t i)           const { return _elems[i]->VTKCellType;  } ///< Return the VTKCellType of an element i
	virtual size_t ElemNVerts      (size_t i)           const { return _elems[i]->V.Size();     } ///< Return the number of vertices of an element
	virtual size_t ElemCon         (size_t i, size_t j) const { return _elems[i]->V[j]->MyID;   } ///< Return the vertex j of an element i
	virtual size_t ElemNETags      (size_t i)           const { return _elems[i]->ETags.Size(); } ///< Return the number of edge tags of an element i
	virtual size_t ElemNFTags      (size_t i)           const { return _elems[i]->FTags.Size(); } ///< Return the number of face tags of an element i
	virtual int    ElemETag        (size_t i, size_t j) const { return _elems[i]->ETags(j);     } ///< Return the edge tag j of an element i
	virtual int    ElemFTag        (size_t i, size_t j) const { return _elems[i]->FTags(j);     } ///< Return the face tag j of an element i

#ifdef USE_BOOST_PYTHON
// {
	void   PyWriteVTU (BPy::str const & FileName) { WriteVTU (BPy::extract<char const *>(FileName)()); }
	size_t PyGetVerts (BPy::list & Verts) const;
	size_t PyGetEdges (BPy::list & Edges) const;
	void   PyGetETags (BPy::dict & ETags) const;
	void   PyGetFTags (BPy::dict & FTags) const;
	size_t PyGetElems (BPy::dict & Elems) const;
// }
#endif // USE_BOOST_PYTHON

protected:
	// Data
	bool           _is_3d;     ///< Is 3D mesh?
	bool           _is_o2;     ///< Is quadratic element?
	Array<Vertex*> _verts;     ///< Vertices
	Array<Elem*>   _elems;     ///< Elements
	Array<Elem*>   _elems_bry; ///< Elements on boundary
	Array<Vertex*> _verts_bry; ///< Vertices on boundary

	// Private methods to be overloaded
	virtual void   _erase            ();                                                  ///< Erase current mesh (deallocate memory)
	virtual size_t _edge_to_lef_vert (size_t EdgeLocalID)             const { return 0; } ///< Returns the local left vertex ID for a given Local Edge ID
	virtual size_t _edge_to_rig_vert (size_t EdgeLocalID)             const { return 0; } ///< Returns the local right vertex ID for a given Local Edge ID
	virtual void   _face_to_verts    (size_t FaceLocalID, Array<size_t> & Verts) const {} ///< Returns the local face vertex IDs for a given Local Edge ID

private:
	// Private methods
	void   _vtk_con (size_t i, String & Connect) const; ///< Returns a string with the connectivites (global vertices IDs) of an element
	size_t _nverts  (int VTKCellType)            const; ///< Returns the number of vertices of a VTKCell
	size_t _nedges  (int VTKCellType)            const; ///< Returns the number of edges of a VTKCell
	size_t _nfaces  (int VTKCellType)            const; ///< Returns the number of faces of a VTKCell

}; // class Generic


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void Generic::WriteVTU(char const * FileName) const
{
	// Results
	std::ostringstream oss;

	// Data
	size_t nn = NVerts(); // Number of Nodes
	size_t ne = NElems(); // Number of Elements

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
		oss << "  " << nsflo <<         VertX(i) << " ";
		oss <<         nsflo <<         VertY(i) << " ";
		oss <<         nsflo << (_is_3d?VertZ(i):0.0);
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
		String con; _vtk_con (i, con);
		oss << "  " << con;
		k++;
		VTU_NEWLINE (i,k,ne,nimax/ElemNVerts(i),oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	size_t offset = 0;
	for (size_t i=0; i<ne; ++i)
	{
		offset += ElemNVerts(i);
		oss << (k==0?"  ":" ") << offset;
		k++;
		VTU_NEWLINE (i,k,ne,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		oss << (k==0?"  ":" ") << ElemVTKCellType(i);
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
		oss << (k==0?"  ":" ") << IsVertOnBry(i);
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
		oss << (k==0?"  ":" ") << IsElemOnBry(i);
		k++;
		VTU_NEWLINE (i,k,ne,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "tag" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		oss << (k==0?"  ":" ") << ElemTag(i);
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

inline void Generic::SetVert(int i, bool IsOnBry, double X, double Y, double Z)
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

inline void Generic::SetElem(int i, int Tag, bool IsOnBry, int VTKCellType)
{
	// Set _elems
	if (_elems[i]==NULL) _elems[i] = new Elem;
	_elems[i]->MyID        = i;
	_elems[i]->Tag         = Tag;
	_elems[i]->OnBry       = IsOnBry;
	_elems[i]->VTKCellType = VTKCellType;
	_elems[i]->V.Resize (_nverts(VTKCellType));

	// Resize ETags and FTags
	int ned = _nedges (VTKCellType);
	int nfa = _nfaces (VTKCellType);
	if (ned>0) _elems[i]->ETags.Resize (ned);
	if (nfa>0) _elems[i]->FTags.Resize (nfa);

	// Set _elems_bry
	if (IsOnBry) _elems_bry.Push (_elems[i]);
}

inline void Generic::SetElemCon(int i, int j, size_t iVert)
{
	_elems[i]->V[j] = _verts[iVert];
}

inline void Generic::SetElemETag(int i, int j, int Tag)
{
	_elems[i]->ETags(j) = Tag;
}

inline void Generic::SetElemFTag(int i, int j, int Tag)
{
	_elems[i]->FTags(j) = Tag;
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
		for (size_t i=0; i<NVerts(); ++i)
			Verts.append (BPy::make_tuple(VertX(i), VertY(i), VertZ(i)));
	}
	else
	{
		for (size_t i=0; i<NVerts(); ++i)
			Verts.append (BPy::make_tuple(VertX(i), VertY(i), 0.0));
	}
	return len(Verts);
}

inline size_t Generic::PyGetEdges(BPy::list & Edges) const
{
	/* Out:
	 *      E = [[v1,v2], [v1,v2], ... num edges]
	 */
	for (size_t i=0; i<NElems(); ++i)
	{
		for (size_t j=0; j<_nedges(ElemVTKCellType(i)); ++j)
		{
			BPy::list pair;
			pair.append  (ElemCon(i, _edge_to_lef_vert(j)));
			pair.append  (ElemCon(i, _edge_to_rig_vert(j)));
			Edges.append (pair);
		}
	}
	return len(Edges);
}

inline void Generic::PyGetETags(BPy::dict & ETags) const
{
	/* Out:
	 *       ETags = {(L,R):tag, (L,R):tag, ... num edge tags]}
	 */
	for (size_t i=0; i<NElems(); ++i)
	{
		for (size_t j=0; j<ElemNETags(i); ++j) // j is the local edge id
		{
			if (ElemETag(i,j)<0)
			{
				int L = ElemCon(i, _edge_to_lef_vert(j));
				int R = ElemCon(i, _edge_to_rig_vert(j));
				ETags[BPy::make_tuple(L, R)] = ElemETag(i,j);
			}
		}
	}
}

inline void Generic::PyGetFTags(BPy::dict & FTags) const
{
	/* Out:
	 *       FTags = {'v1_v2_v3_v4_...':tag, 'v1_v2_v3_v4_...':tag, ... num face tags]}
	 */
	for (size_t i=0; i<NElems(); ++i)
	{
		for (size_t j=0; j<ElemNFTags(i); ++j) // j is the local face id
		{
			if (ElemFTag(i,j)<0)
			{
				Array<size_t> fv; // face verts
				_face_to_verts (j, fv);
				String key; for (size_t k=0; k<fv.Size(); ++k) key.Printf("%s_%d", fv[k]);
				FTags[key] = ElemFTag(i,j);
			}
		}
	}
}

inline size_t Generic::PyGetElems(BPy::dict & Elems) const
{
	/* Out:
	 *    
	 *    Elems = {
	 *      'tags'    : [t1,t2,t3, ... num elements]
	 *      'onbs'    : [ 1, 0, 1, ... num elements]
	 *      'vtks'    : [ 9, 9,12, ... num elements]
	 *      'cons'    : {'id0':[node0, node1, node3, num nodes in element 0],
	 *                   'id1':[n0,n1,n2,n3,n4,n5,n6,n7], ... num elems}
	 *      'etags'   : {'id0':[tag_edge_0, tag_edge_1, tag_edge_2, tag_edge_3],
	 *                   'id1':[t0,t1,t2,t3], ... num elems]}
	 *      'ftags'   : {'id0':[tag_face_0, tag_face_1, tag_face_2, tag_face_3, tag_face_4, tag_face_5],
	 *                   'id1':[t0,t1,t2,t3,t4,t5], ... num elems]}
	 *    }
	 */
	BPy::list tags,  onbs,  vtks;
	BPy::dict cons, etags, ftags;
	for (size_t i=0; i<NElems(); ++i)
	{
		// Data
		tags.append (ElemTag(i));
		onbs.append (IsElemOnBry(i) ? 1 : 0);
		vtks.append (ElemVTKCellType(i));

		// Connectivities
		BPy::list co;
		for (size_t j=0; j<ElemNVerts(i); ++j)
			co.append (ElemCon(i,j));
		cons[BPy::str(i)] = co;

		// ETags
		BPy::list et;
		for (size_t j=0; j<ElemNETags(i); ++j) // j is the local edge id
		{
			int tag = ElemETag(i,j);
			et.append (tag);
		}
		etags[BPy::str(i)] = et;

		// FTags
		if (Is3D())
		{
			BPy::list ft;
			for (size_t j=0; j<ElemNFTags(i); ++j) // j is the local face id
			{
				int tag = ElemFTag(i,j);
				ft.append (tag);
			}
			ftags[BPy::str(i)] = ft;
		}
	}
	Elems["tags" ] = tags;
	Elems["onbs" ] = onbs;
	Elems["vtks" ] = vtks;
	Elems["cons" ] = cons;
	Elems["etags"] = etags;
	Elems["ftags"] = ftags;
	return NElems();
}

// }
#endif // USE_BOOST_PYTHON


/* private */

inline void Generic::_erase()
{
	for (size_t i=0; i<_verts.Size(); ++i) if (_verts[i]!=NULL) delete _verts[i]; // it is only necessary to delete nodes in _verts array
	for (size_t i=0; i<_elems.Size(); ++i) if (_elems[i]!=NULL) delete _elems[i]; // it is only necessary to delete elems in _elems array
	_is_3d = false;
	if (_verts    .Size()>0) _verts      .Resize(0);
	if (_elems    .Size()>0) _elems      .Resize(0);
	if (_elems_bry.Size()>0) _elems_bry  .Resize(0);
	if (_verts_bry.Size()>0) _verts_bry  .Resize(0);
}

inline void Generic::_vtk_con(size_t i, String & Connect) const
{
	Connect = "";
	for (size_t j=0; j<ElemNVerts(i); j++)
		Connect.Printf ("%s %d ", Connect.CStr(), ElemCon(i,j));
}

inline size_t Generic::_nverts(int VTKCellType) const
{
	switch (VTKCellType)
	{
		case VTK_LINE:                 { return  2; }
		case VTK_TRIANGLE:             { return  3; }
		case VTK_QUAD:                 { return  4; }
		case VTK_TETRA:                { return  4; }
		case VTK_HEXAHEDRON:           { return  8; }
		case VTK_QUADRATIC_EDGE:       { return  3; }
		case VTK_QUADRATIC_TRIANGLE:   { return  6; }
		case VTK_QUADRATIC_QUAD:       { return  8; }
		case VTK_QUADRATIC_TETRA:      { return 10; }
		case VTK_QUADRATIC_HEXAHEDRON: { return 20; }
		default: throw new Fatal("Generic::_nverts: VTKCellType==%d is invalid", VTKCellType);
	}
}

inline size_t Generic::_nedges(int VTKCellType) const
{
	switch (VTKCellType)
	{
		case VTK_LINE:                 { return  1; }
		case VTK_TRIANGLE:             { return  3; }
		case VTK_QUAD:                 { return  4; }
		case VTK_TETRA:                { return  6; }
		case VTK_HEXAHEDRON:           { return 12; }
		case VTK_QUADRATIC_EDGE:       { return  1; }
		case VTK_QUADRATIC_TRIANGLE:   { return  3; }
		case VTK_QUADRATIC_QUAD:       { return  4; }
		case VTK_QUADRATIC_TETRA:      { return  6; }
		case VTK_QUADRATIC_HEXAHEDRON: { return 12; }
		default: throw new Fatal("Generic::_nedges: VTKCellType==%d is invalid", VTKCellType);
	}
}

inline size_t Generic::_nfaces(int VTKCellType) const
{
	switch (VTKCellType)
	{
		case VTK_LINE:                 { return 0; }
		case VTK_TRIANGLE:             { return 0; }
		case VTK_QUAD:                 { return 0; }
		case VTK_TETRA:                { return 4; }
		case VTK_HEXAHEDRON:           { return 6; }
		case VTK_QUADRATIC_EDGE:       { return 0; }
		case VTK_QUADRATIC_TRIANGLE:   { return 0; }
		case VTK_QUADRATIC_QUAD:       { return 0; }
		case VTK_QUADRATIC_TETRA:      { return 4; }
		case VTK_QUADRATIC_HEXAHEDRON: { return 6; }
		default: throw new Fatal("Generic::_nfaces: VTKCellType==%d is invalid", VTKCellType);
	}
}


/* output */

std::ostream & operator<< (std::ostream & os, Mesh::Generic const & G)
{
	for (size_t i=0; i<G.NElems(); ++i)
	{
		os << "[" << i << "] " << G.ElemTag(i) << " " << G.IsElemOnBry(i);
		os << " {"; for (size_t j=0; j<G.ElemNETags(i); ++j) os << G.ElemETag(i,j) << (j==G.ElemNETags(i)-1?"":","); os << "}";
		os << " {"; for (size_t j=0; j<G.ElemNFTags(i); ++j) os << G.ElemFTag(i,j) << (j==G.ElemNFTags(i)-1?"":","); os << "}";
		os << "\n";
		for (size_t j=0; j<G.ElemNVerts(i); ++j)
		{
			size_t k = G.ElemCon (i, j);
			os << "   [" << k << "] " << G.IsVertOnBry(k);
			if (G.Is3D()) os << " (" << Util::_8_4<< G.VertX(k) << Util::_8_4<< G.VertY(k) << Util::_8_4<< G.VertZ(k) << ")";
			else          os << " (" << Util::_8_4<< G.VertX(k) << Util::_8_4<< G.VertY(k) << ")";
			os << "\n";
		}
	}
	return os;
}

}; // namespace Mesh

#endif // MECHSYS_MESH_H
