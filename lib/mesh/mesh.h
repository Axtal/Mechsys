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
  //#include <boost/python.hpp> // this includes everything
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

using std::cout;
using std::endl;

namespace Mesh
{

struct Edge
{
	size_t L; // Left vertex id
	size_t R; // Right vertex id
};

struct Face
{
	size_t I0; // Edge or Vertex local id
	size_t I1;
	size_t I2;
	size_t I3;
	size_t I4;
	size_t I5;
	size_t I6;
	size_t I7;
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
	// Constants
	static Edge TRI_EDGE2VERT[]; ///< Triangle: Map from local edge ID to local vertex ID
	static Face TET_FACE2VERT[]; ///< Tetrahedron: Map from local face ID to local vertex ID
	static Face TET_FACE2EDGE[]; ///< Tetrahedron: Map from local face ID to local edge ID
	static Edge QUA_EDGE2VERT[]; ///< Quad: Map from local edge ID to local vertex ID
	static Face HEX_FACE2VERT[]; ///< Hex: Map from local face ID to local vertex IDs
	static Face HEX_FACE2EDGE[]; ///< Hex: Map from local face ID to local edge IDs

	// Constructor
	Generic (bool Is3D) : _is_3d(Is3D), _is_o2(false) {}

	// Destructor
	virtual ~Generic () { _erase(); }

	// Methods
	void WriteVTU (char const * FileName) const; ///< Write output file for ParaView

	// Methods
	size_t EdgeToLef (size_t iElem, size_t EdgeLocalID) const; ///< Returns the GLOBAL left vertex ID for a given Local Edge ID and Element # i
	size_t EdgeToRig (size_t iElem, size_t EdgeLocalID) const; ///< Returns the GLOBAL right vertex ID for a given Local Edge ID and Element # i

	// Set methods
	virtual void SetO2       (bool IsO2=true) { _is_o2=IsO2; }                      ///< (Un)set quadratic elements
	virtual void SetNVerts   (size_t NumVerts);                                     ///< Erase old mesh and set number of vertices
	virtual void SetNElems   (size_t NumElems);                                     ///< Set number of elements
	virtual void SetVert     (int i, bool IsOnBry, double X, double Y, double Z=0); ///< Set vertex
	virtual void SetElem     (int i, int Tag, bool IsOnBry, int VTKCellType);       ///< Set element
	virtual void SetElemCon  (int i, int j, size_t iVert);                          ///< Set element connectivity
	virtual void SetElemETag (int i, int j, int Tag);                               ///< Set element's edge tag
	virtual void SetElemFTag (int i, int j, int Tag);                               ///< Set element's face tag

	// Beams and Joints
	void   SetETagsBeams (size_t NBeamETags, ...);                ///< Set what edge tags represent Beams
	size_t nETagsBeams   () const { return _etags_beams.Size(); } ///< Return the number of edge tags that represent Beams
	bool   IsETagBeam    (int ETag) const;                        ///< Find whether an EdgeTag represents a Beam or not

	// Get methods
	        void   ElemCentre      (size_t i, double & X, double & Y, double & Z) const;          ///< Calculate centre of an element i
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
	void   PyWriteVTU    (BPy::str const & FileName) { WriteVTU (BPy::extract<char const *>(FileName)()); }
	size_t PyGetVerts    (BPy::list & Verts)      const;
	size_t PyGetVertsBry (BPy::list & VertsOnBry) const;
	size_t PyGetEdges    (BPy::list & Edges)      const;
	void   PyGetETags    (BPy::dict & ETags)      const;
	void   PyGetFTags    (BPy::dict & FTags)      const;
	size_t PyGetElems    (BPy::dict & Elems)      const; ///< Return: Elems = { "ID1":[Tag,VTK,Xcentre,Ycentre,Zcentre], "ID2":[Tag,VTK,X,Y,Z], num elements }
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
	Array<int>     _etags_beams; ///< Edge tags that represent Beams

	// Private methods that MAY be overloaded
	virtual void _erase (); ///< Erase current mesh (deallocate memory)

	// Private methods
	size_t _edge_to_lef_vert (size_t iElem, size_t EdgeLocalID) const;                        ///< Returns the local left vertex ID for a given Local Edge ID of element # iElem
	size_t _edge_to_rig_vert (size_t iElem, size_t EdgeLocalID) const;                        ///< Returns the local right vertex ID for a given Local Edge ID of element # iElem
	void   _face_to_verts    (size_t iElem, size_t FaceLocalID, Array<size_t> & Verts) const; ///< Returns the local vertex IDs for a given Local Face ID
	void   _face_to_edges    (size_t iElem, size_t FaceLocalID, Array<size_t> & Edges) const; ///< Returns the local edge IDs for a given Local Face ID
	void   _vtk_con          (size_t iElem, String & Connect) const;                          ///< Returns a string with the connectivites (global vertices IDs) of an element

	// Private methods
	size_t _nverts  (int VTKCellType) const; ///< Returns the number of vertices of a VTKCell
	size_t _nedges  (int VTKCellType) const; ///< Returns the number of edges of a VTKCell
	size_t _nfaces  (int VTKCellType) const; ///< Returns the number of faces of a VTKCell

	// Beam
	void _add_beams ();

}; // class Generic

Edge Generic::TRI_EDGE2VERT[]= {{ 0, 1 },
                                { 1, 2 },
                                { 0, 2 }};

Face Generic::TET_FACE2VERT[]= {{  0,  2,  3,  6,  9,  7 },  // Face # 0 => Vertices 0,2,3...
                                {  0,  1,  3,  4,  8,  7 },  // Face # 1
                                {  0,  1,  2,  4,  5,  6 },  // Face # 2
                                {  1,  2,  3,  5,  9,  8 }}; // Face # 3

Face Generic::TET_FACE2EDGE[]= {{  0,  3,  5,  6,  9, 11 },  // Face # 0 => Edges 0,4,8
                                {  1,  3,  4,  7,  9, 10 },  // Face # 1
                                {  0,  1,  2,  6,  7,  8 },  // Face # 2
                                {  2,  4,  5,  8, 10, 11 }}; // Face # 3

Edge Generic::QUA_EDGE2VERT[]= {{ 0, 3 },  // Edge #  0
                                { 1, 2 },  // Edge #  1
                                { 0, 1 },  // Edge #  2
                                { 2, 3 },  // Edge #  3
                                { 4, 7 },  // Edge #  4
                                { 5, 6 },  // Edge #  5
                                { 4, 5 },  // Edge #  6
                                { 6, 7 },  // Edge #  7
                                { 0, 4 },  // Edge #  8
                                { 1, 5 },  // Edge #  9
                                { 2, 6 },  // Edge # 10
                                { 3, 7 }}; // Edge # 11

Face Generic::HEX_FACE2VERT[]= {{  0,  3,  7,  4, 11, 19, 15, 16 },  // Face # 0 => Vertices 0,3,7...
                                {  1,  2,  6,  5,  9, 18, 13, 17 },  // Face # 1
                                {  0,  1,  5,  4,  8, 17, 12, 16 },  // Face # 2
                                {  2,  3,  7,  6, 10, 19, 14, 18 },  // Face # 3
                                {  0,  1,  2,  3,  8,  9, 10, 11 },  // Face # 4
                                {  4,  5,  6,  7, 12, 13, 14, 15 }}; // Face # 5

Face Generic::HEX_FACE2EDGE[]= {{  0,  4,  8, 11, 12, 16, 20, 23 },  // Face # 0 => Edges 0,4,8
                                {  1,  5,  9, 10, 13, 17, 21, 22 },  // Face # 1
                                {  2,  6,  8,  9, 14, 18, 20, 21 },  // Face # 2
                                {  3,  7, 10, 11, 15, 19, 22, 23 },  // Face # 3
                                {  0,  1,  2,  3, 12, 13, 14, 15 },  // Face # 4
                                {  4,  5,  6,  7, 16, 17, 18, 19 }}; // Face # 5

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

inline size_t Generic::EdgeToLef (size_t iElem, size_t EdgeLocalID) const
{
	return ElemCon(iElem, _edge_to_lef_vert(iElem, EdgeLocalID));
}

inline size_t Generic::EdgeToRig (size_t iElem, size_t EdgeLocalID) const
{
	return ElemCon(iElem, _edge_to_rig_vert(iElem, EdgeLocalID));
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

inline void Generic::SetETagsBeams(size_t NBeamETags, ...)
{
	_etags_beams.Resize(NBeamETags);
	va_list arg_list;
	va_start (arg_list, NBeamETags); // initialize arg_list with parameters AETER NBeamETags
	for (size_t i=0; i<NBeamETags; ++i)
		_etags_beams[i] = va_arg (arg_list, int);
	va_end (arg_list);
}

inline bool Generic::IsETagBeam(int ETag) const
{
	return (_etags_beams.Find(ETag)>=0);
}

inline void Generic::ElemCentre (size_t i, double & X, double & Y, double & Z) const
{
	if (ElemVTKCellType(i)==VTK_LINE)
	{
		X =           (VertX(ElemCon(i,0))+VertX(ElemCon(i,1)))/2.0;
		Y =           (VertY(ElemCon(i,0))+VertY(ElemCon(i,1)))/2.0;
		Z = (_is_3d ? (VertZ(ElemCon(i,0))+VertZ(ElemCon(i,1)))/2.0 : 0.0);
	}
	else if (ElemVTKCellType(i)==VTK_TRIANGLE || ElemVTKCellType(i)==VTK_QUADRATIC_TRIANGLE)
	{
		X =           (VertX(ElemCon(i,0))+VertX(ElemCon(i,1))+VertX(ElemCon(i,2)))/3.0;
		Y =           (VertY(ElemCon(i,0))+VertY(ElemCon(i,1))+VertY(ElemCon(i,2)))/3.0;
		Z = (_is_3d ? (VertZ(ElemCon(i,0))+VertZ(ElemCon(i,1))+VertZ(ElemCon(i,2)))/3.0 : 0.0);
	}
	else if (ElemVTKCellType(i)==VTK_QUAD || ElemVTKCellType(i)==VTK_QUADRATIC_QUAD)
	{
		X =           (VertX(ElemCon(i,0))+VertX(ElemCon(i,2)))/2.0;
		Y =           (VertY(ElemCon(i,0))+VertY(ElemCon(i,2)))/2.0;
		Z = (_is_3d ? (VertZ(ElemCon(i,0))+VertZ(ElemCon(i,2)))/2.0 : 0.0);
	}
	//else if (ElemVTKCellType(i)==VTK_TETRA)
	else if (ElemVTKCellType(i)==VTK_HEXAHEDRON || ElemVTKCellType(i)==VTK_QUADRATIC_HEXAHEDRON)
	{
		X =           (VertX(ElemCon(i,0))+VertX(ElemCon(i,6)))/2.0;
		Y =           (VertY(ElemCon(i,0))+VertY(ElemCon(i,6)))/2.0;
		Z = (_is_3d ? (VertZ(ElemCon(i,0))+VertZ(ElemCon(i,6)))/2.0 : 0.0);
	}
	//else if (ElemVTKCellType(i)==VTK_QUADRATIC_EDGE)
	//else if (ElemVTKCellType(i)==VTK_QUADRATIC_TETRA)
	else throw new Fatal("Generic::ElemCentre: VTKCellType==%d is invalid (not implemented yet)", ElemVTKCellType(i));
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

inline size_t Generic::PyGetVertsBry(BPy::list & VertsOnBry) const
{
	/* Out:
	 *      Verts = [id0, id1, ... num verts on bry]
	 */
	for (size_t i=0; i<NVertsBry(); ++i)
		VertsOnBry.append (VertBry(i));
	return NVertsBry();
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
			pair.append  (ElemCon(i, _edge_to_lef_vert(i, j)));
			pair.append  (ElemCon(i, _edge_to_rig_vert(i, j)));
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
				int L = ElemCon(i, _edge_to_lef_vert(i, j));
				int R = ElemCon(i, _edge_to_rig_vert(i, j));
				ETags[BPy::make_tuple(L, R)] = ElemETag(i,j);
			}
		}
	}
}

inline void Generic::PyGetFTags(BPy::dict & FTags) const
{
	/* Out:
	 *                 edg0   edg1            edg1  edg2
	 *       FTags = {'v1,v2_v3,v4'...:tag, 'v1,v2_v3,v4'...:tag, ... num face tags]}
	 */
	for (size_t i=0; i<NElems(); ++i)
	{
		for (size_t j=0; j<ElemNFTags(i); ++j) // j is the local face id
		{
			if (ElemFTag(i,j)<0)
			{
				Array<size_t> fe; // face-edges
				_face_to_edges (i, j, fe);
				String key;
				for (size_t k=0; k<fe.Size(); ++k)
				{
					if (k==0) key.Printf(   "%d,%d",             ElemCon(i, _edge_to_lef_vert(i, fe[k])), ElemCon(i, _edge_to_rig_vert(i, fe[k])));
					else      key.Printf("%s_%d,%d", key.CStr(), ElemCon(i, _edge_to_lef_vert(i, fe[k])), ElemCon(i, _edge_to_rig_vert(i, fe[k])));
				}
				FTags[key.CStr()] = ElemFTag(i,j);
			}
		}
	}
}

inline size_t Generic::PyGetElems(BPy::dict & Elems) const
{
	for (size_t i=0; i<NElems(); ++i)
	{
		// Calcuate Centre
		double X=0.0;
		double Y=0.0;
		double Z=0.0;
		ElemCentre (i, X, Y, Z);

		// Add to dict
		BPy::list info;
		info.append (ElemTag(i));
		info.append (ElemVTKCellType(i));
		info.append (X);
		info.append (Y);
		info.append (Z);
		String id;  id.Printf("%d",i);
		Elems[id.CStr()] = info;
	}
	return NElems();
}

// }
#endif // USE_BOOST_PYTHON


/* private */

inline void Generic::_erase()
{
	for (size_t i=0; i<_verts.Size(); ++i) if (_verts[i]!=NULL) delete _verts[i]; // it is only necessary to delete nodes in _verts array
	for (size_t i=0; i<_elems.Size(); ++i) if (_elems[i]!=NULL) delete _elems[i]; // it is only necessary to delete elems in _elems array
	if (_verts    .Size()>0) _verts      .Resize(0);
	if (_elems    .Size()>0) _elems      .Resize(0);
	if (_elems_bry.Size()>0) _elems_bry  .Resize(0);
	if (_verts_bry.Size()>0) _verts_bry  .Resize(0);
}

inline size_t Generic::_edge_to_lef_vert (size_t iElem, size_t EdgeLocalID) const
{
	switch (ElemVTKCellType(iElem))
	{
		case VTK_LINE:                 { return  0; }
		case VTK_TRIANGLE:             { return  TRI_EDGE2VERT[EdgeLocalID].L; }
		case VTK_QUAD:                 { return  QUA_EDGE2VERT[EdgeLocalID].L; }
		case VTK_TETRA:                { throw new Fatal("Generic::_edge_to_lef_vert: Method not available for Tetrahedrons"); }
		case VTK_HEXAHEDRON:           { return  QUA_EDGE2VERT[EdgeLocalID].L; }
		case VTK_QUADRATIC_EDGE:       { return  0; }
		case VTK_QUADRATIC_TRIANGLE:   { return  TRI_EDGE2VERT[EdgeLocalID].L; }
		case VTK_QUADRATIC_QUAD:       { return  QUA_EDGE2VERT[EdgeLocalID].L; }
		case VTK_QUADRATIC_TETRA:      { throw new Fatal("Generic::_edge_to_lef_vert: Method not available for Quadratic Tetrahedrons"); }
		case VTK_QUADRATIC_HEXAHEDRON: { throw new Fatal("Generic::_edge_to_lef_vert: Method not available for Quadratic Hexahedrons"); }
		default: throw new Fatal("Generic::_edge_to_lef_vert: VTKCellType==%d is invalid (not implemented yet)", ElemVTKCellType(iElem));
	}
}

inline size_t Generic::_edge_to_rig_vert (size_t iElem, size_t EdgeLocalID) const
{
	switch (ElemVTKCellType(iElem))
	{
		case VTK_LINE:                 { return  1; }
		case VTK_TRIANGLE:             { return  TRI_EDGE2VERT[EdgeLocalID].R; }
		case VTK_QUAD:                 { return  QUA_EDGE2VERT[EdgeLocalID].R; }
		case VTK_TETRA:                { throw new Fatal("Generic::_edge_to_rig_vert: Method not available for Tetrahedrons"); }
		case VTK_HEXAHEDRON:           { throw new Fatal("Generic::_edge_to_rig_vert: Method not available for Hexahedrons"); }
		case VTK_QUADRATIC_EDGE:       { return  1; }
		case VTK_QUADRATIC_TRIANGLE:   { return  TRI_EDGE2VERT[EdgeLocalID].R; }
		case VTK_QUADRATIC_QUAD:       { return  QUA_EDGE2VERT[EdgeLocalID].R; }
		case VTK_QUADRATIC_TETRA:      { throw new Fatal("Generic::_edge_to_rig_vert: Method not available for Quadratic Tetrahedrons"); }
		case VTK_QUADRATIC_HEXAHEDRON: { throw new Fatal("Generic::_edge_to_rig_vert: Method not available for Quadratic Hexahedrons"); }
		default: throw new Fatal("Generic::_edge_to_rig_vert: VTKCellType==%d is invalid (not implemented yet)", ElemVTKCellType(iElem));
	}
}

inline void Generic::_face_to_verts(size_t iElem, size_t FaceLocalID, Array<size_t> & Verts) const
{
	if (ElemVTKCellType(iElem)==VTK_TETRA)
	{
		Verts.Resize(3);
		Verts[0] = TET_FACE2VERT[FaceLocalID].I0;
		Verts[1] = TET_FACE2VERT[FaceLocalID].I1;
		Verts[2] = TET_FACE2VERT[FaceLocalID].I2;
	}
	else if (ElemVTKCellType(iElem)==VTK_QUADRATIC_TETRA)
	{
		Verts.Resize(6);
		Verts[0] = TET_FACE2VERT[FaceLocalID].I0;
		Verts[1] = TET_FACE2VERT[FaceLocalID].I1;
		Verts[2] = TET_FACE2VERT[FaceLocalID].I2;
		Verts[3] = TET_FACE2VERT[FaceLocalID].I3;
		Verts[4] = TET_FACE2VERT[FaceLocalID].I4;
		Verts[5] = TET_FACE2VERT[FaceLocalID].I5;
	}
	else if (ElemVTKCellType(iElem)==VTK_HEXAHEDRON)
	{
		Verts.Resize(4);
		Verts[0] = HEX_FACE2VERT[FaceLocalID].I0;
		Verts[1] = HEX_FACE2VERT[FaceLocalID].I1;
		Verts[2] = HEX_FACE2VERT[FaceLocalID].I2;
		Verts[3] = HEX_FACE2VERT[FaceLocalID].I3;
	}
	else if (ElemVTKCellType(iElem)==VTK_QUADRATIC_HEXAHEDRON)
	{
		Verts.Resize(8);
		Verts[0] = HEX_FACE2VERT[FaceLocalID].I0;
		Verts[1] = HEX_FACE2VERT[FaceLocalID].I1;
		Verts[2] = HEX_FACE2VERT[FaceLocalID].I2;
		Verts[3] = HEX_FACE2VERT[FaceLocalID].I3;
		Verts[4] = HEX_FACE2VERT[FaceLocalID].I4;
		Verts[5] = HEX_FACE2VERT[FaceLocalID].I5;
		Verts[6] = HEX_FACE2VERT[FaceLocalID].I6;
		Verts[7] = HEX_FACE2VERT[FaceLocalID].I7;
	}
	else throw new Fatal("Generic::_face_to_verts: Method not available for VTKCellType==%d",ElemVTKCellType(iElem));
}

inline void Generic::_face_to_edges(size_t iElem, size_t FaceLocalID, Array<size_t> & Edges) const
{
	if (ElemVTKCellType(iElem)==VTK_TETRA)
	{
		Edges.Resize(3);
		Edges[0] = TET_FACE2EDGE[FaceLocalID].I0;
		Edges[1] = TET_FACE2EDGE[FaceLocalID].I1;
		Edges[2] = TET_FACE2EDGE[FaceLocalID].I2;
	}
	else if (ElemVTKCellType(iElem)==VTK_QUADRATIC_TETRA)
	{
		Edges.Resize(6);
		Edges[0] = TET_FACE2EDGE[FaceLocalID].I0;
		Edges[1] = TET_FACE2EDGE[FaceLocalID].I1;
		Edges[2] = TET_FACE2EDGE[FaceLocalID].I2;
		Edges[3] = TET_FACE2EDGE[FaceLocalID].I3;
		Edges[4] = TET_FACE2EDGE[FaceLocalID].I4;
		Edges[5] = TET_FACE2EDGE[FaceLocalID].I5;
	}
	else if (ElemVTKCellType(iElem)==VTK_HEXAHEDRON)
	{
		Edges.Resize(4);
		Edges[0] = HEX_FACE2EDGE[FaceLocalID].I0;
		Edges[1] = HEX_FACE2EDGE[FaceLocalID].I1;
		Edges[2] = HEX_FACE2EDGE[FaceLocalID].I2;
		Edges[3] = HEX_FACE2EDGE[FaceLocalID].I3;
	}
	else if (ElemVTKCellType(iElem)==VTK_QUADRATIC_HEXAHEDRON)
	{
		Edges.Resize(8);
		Edges[0] = HEX_FACE2EDGE[FaceLocalID].I0;
		Edges[1] = HEX_FACE2EDGE[FaceLocalID].I1;
		Edges[2] = HEX_FACE2EDGE[FaceLocalID].I2;
		Edges[3] = HEX_FACE2EDGE[FaceLocalID].I3;
		Edges[4] = HEX_FACE2EDGE[FaceLocalID].I4;
		Edges[5] = HEX_FACE2EDGE[FaceLocalID].I5;
		Edges[6] = HEX_FACE2EDGE[FaceLocalID].I6;
		Edges[7] = HEX_FACE2EDGE[FaceLocalID].I7;
	}
	else throw new Fatal("Generic::_face_to_edges: Method not available for VTKCellType==%d",ElemVTKCellType(iElem));
}

inline void Generic::_vtk_con(size_t i, String & Connect) const
{
	if (ElemVTKCellType(i)==VTK_LINE)
	{
		Connect.Printf("%d %d",ElemCon(i,0), ElemCon(i,1));
	}
	else if (ElemVTKCellType(i)==VTK_QUADRATIC_EDGE)
	{
		throw new Fatal("Generic::_vtk_con: Quadratic Linear elements are not available");
	}
	else if (ElemVTKCellType(i)==VTK_TRIANGLE)
	{
		Connect.Printf("%d %d %d",ElemCon(i,0), ElemCon(i,1), ElemCon(i,2));
	}
	else if (ElemVTKCellType(i)==VTK_QUADRATIC_TRIANGLE)
	{
		Connect.Printf("%d %d %d %d %d %d",ElemCon(i,0), ElemCon(i,1), ElemCon(i,2),
		                                   ElemCon(i,3), ElemCon(i,4), ElemCon(i,5));
	}
	else if (ElemVTKCellType(i)==VTK_QUAD)
	{
		Connect.Printf("%d %d %d %d",ElemCon(i,0), ElemCon(i,1), ElemCon(i,2), ElemCon(i,3));
	}
	else if (ElemVTKCellType(i)==VTK_QUADRATIC_QUAD)
	{
		Connect.Printf("%d %d %d %d %d %d %d %d",
		               ElemCon(i,0), ElemCon(i,1), ElemCon(i,2), ElemCon(i,3),
		               ElemCon(i,4), ElemCon(i,5), ElemCon(i,6), ElemCon(i,7));
	}
	else if (ElemVTKCellType(i)==VTK_TETRA)
	{
		throw new Fatal("Generic::_vtk_con: Tetrahedrons elements are not available");
	}
	else if (ElemVTKCellType(i)==VTK_QUADRATIC_TETRA)
	{
		throw new Fatal("Generic::_vtk_con: Quadratic Tetrahedrons elements are not available");
	}
	else if (ElemVTKCellType(i)==VTK_HEXAHEDRON)
	{
		Connect.Printf("%d %d %d %d %d %d %d %d",ElemCon(i,0), ElemCon(i,1), ElemCon(i,2), ElemCon(i,3),
		                                         ElemCon(i,4), ElemCon(i,5), ElemCon(i,6), ElemCon(i,7));
	}
	else if (ElemVTKCellType(i)==VTK_QUADRATIC_HEXAHEDRON)
	{
		Connect.Printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
		               ElemCon(i, 0), ElemCon(i, 1), ElemCon(i, 2), ElemCon(i, 3),
		               ElemCon(i, 4), ElemCon(i, 5), ElemCon(i, 6), ElemCon(i, 7),
		               ElemCon(i, 8), ElemCon(i, 9), ElemCon(i,10), ElemCon(i,11),
		               ElemCon(i,12), ElemCon(i,13), ElemCon(i,14), ElemCon(i,15),
		               ElemCon(i,16), ElemCon(i,17), ElemCon(i,18), ElemCon(i,19));
	}
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

inline void Generic::_add_beams()
{
	// Add Beam elements
	if (nETagsBeams()>0)
	{
		for (size_t i=0; i<NElems(); ++i)
		{
			for (size_t j=0; j<ElemNETags(i); ++j)
			{
				int etag = ElemETag (i,j);
				if (etag<0)
				{
					if (IsETagBeam(etag))
					{
						cout << etag << endl;
					}
				}
			}
		}
		//for (int p=0; p<e->ETags.Size(); ++p) // p is EdgeLocalID
		//{
			//int etag = e->ETags(m);
			//if (IsETagBeam(etag))
			//{
				// Add new Beam element
				//size_t l = EdgeToLef (_elems.Size()-1, p);
				//size_t m = EdgeToMid (_elems.Size()-1, p);
				//size_t r = EdgeToRig (_elems.Size()-1, p);
				//Elem * beam = new Elem;
				//e->MyID = _elems.Size(); // id
				//e->Tag  = etag;          // tag
				// connectivity
				//e->V.Resize((_is_o2?3:8));
			//}
		//}
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
