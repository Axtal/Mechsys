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

#ifndef MECHSYS_MESH_STRUCTURED_H
#define MECHSYS_MESH_STRUCTURED_H

/* LOCAL indexes of Vertices, Edges, and Faces

  2D:
                Vertices                               Edges
   y
   |        3      6      2                             y+
   +--x      @-----@-----@                        +----(3)----+
             |           |                        |           |
             |           |                        |           |
           7 @           @ 5                  x- (0)         (1) x+
             |           |                        |           |
             |           |                        |           |
             @-----@-----@                        +----(2)----+
            0      4      1                             y-

  3D:
                  Vertices                             Edges                              Faces
    z
    |           4        15        7
   ,+--y         @_______@________@                 +_______(4)______+                 +________________+ 
 x'            ,'|              ,'|               ,'|              ,'|               ,'|              ,'| 
          12 @'  |         14 ,'  |             ,'  |            ,'  |             ,'  |  ___       ,'  | 
           ,'    |16        ,@    |19         (6)  (8)         (7)   |           ,'    |,'5,'  [0],'    | 
     5   ,'      @      6 ,'      @         ,'      |        ,'     (11)       ,'      |~~~     ,'      | 
       @'=======@=======@'        |       +'==========(5)==+'        |       +'===============+'  ,'|   | 
       |      13 |      |         |       |         |      |         |       |   ,'|   |      |   |3|   | 
       |         |      |  11     |       |         |      |         |       |   |2|   |      |   |,'   | 
    17 |       0 @______|_@_______@       |         +______|_(0)_____+       |   |,'   +______|_________+ 
       @       ,'       @       ,' 3     (9)      ,'       |       ,'        |       ,'       |       ,'  
       |   8 @'      18 |     ,'          |    (2)        (10)   ,'          |     ,' [1]  ___|     ,'    
       |   ,'           |   ,@ 10         |   ,'           |   (3)           |   ,'      ,'4,'|   ,'      
       | ,'             | ,'              | ,'             | ,'              | ,'        ~~~  | ,'        
       @_______@________@'                +______(1)_______+'                +________________+'          
     1         9         2
*/

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
#include "util/lineparser.h"
#include "util/tree.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/laexpr.h"
#include "mesh/mesh.h"

using LinAlg::Vector;
using LinAlg::Matrix;
using blitz::TinyVector;

namespace Mesh
{

	/* TODO:
	 *        1) Add additional checks, especially for 3D meshes
	 *        2) Remove boundary information for nodes duplicated not on boundary
	 *        3) Add boundary marks evaluation
	 *        4) Extend to second-order (o2) elements
	 *        5) Add quality check and improvement
	 */

struct Pair { int L; int R; };

class Block
{
public:
	/* 2D: 8 nodes => C. Resize(2, 8)
	 *                Wx.Resize(nDivX)
	 *                Wy.Resize(nDivY)
	 *          _                         _
	 *    C  = |  x0 x1 x2 x3 x4 x5 x6 x7  |
	 *         |_ y0 y1 y2 y3 y4 y5 y6 y7 _|
	 *
	 * 3D: 20 nodes => C. Resize(3, 20)
	 *                 Wx.Resize(nDivX)
	 *                 Wy.Resize(nDivY)
	 *                 Wz.Resize(nDivZ)
	 *         _                                         _
	 *        |  x0 x1 x2 x3 x4 x5 x6 x7 ... x17 x18 x19  |
	 *    C = |  y0 y1 y2 y3 y4 y5 y6 y7 ... y17 y18 y19  |
	 *        |_ z0 z1 z2 z3 z4 z5 z6 z7 ... z17 z18 z19 _|
	 */
	// Constants
	static Pair Edge2Face[]; ///< Map from local edge ID to a pair of faces (local ID) that share this edge

	// Constructor
	Block () : _n_div_x(0), _n_div_y(0), _n_div_z(0), _tag(-1), _has_etags(false), _has_ftags(false) { Set2D(); }

	// Set methods
	void             Set2D  ();                                           ///< Set 2D
	void             Set3D  ();                                           ///< Set 2D
	void             SetTag (int Tag) { _tag = Tag; }                     ///< Set the tag to be replicated to all elements generated inside this block
	Matrix<double> & C      ()        { return _c;  }                     ///< Set coordinates
	Vector<int>    & ETags  ()        { _has_etags=true; return _etags; } ///< Set edge tags
	Vector<int>    & FTags  ()        { _has_ftags=true; return _ftags; } ///< Set face tags
	void             SetWx  (char const * Str);                           ///< Set x weights and the number of divisions along x. Ex: Str = "1 1 1 2 2 2 3 3 3"
	void             SetWy  (char const * Str);                           ///< Set y weights and the number of divisions along y. Ex: Str = "1 2 3 4"
	void             SetWz  (char const * Str);                           ///< Set z weights and the number of divisions along z. Ex: Str = "1 1 1 1"

	// Access methods
	int    Tag        ()         const { return _tag;          }
	bool   Is3D       ()         const { return _is_3d;        }
	int    nDivX      ()         const { return _n_div_x;      }
	int    nDivY      ()         const { return _n_div_y;      }
	int    nDivZ      ()         const { return _n_div_z;      }
	double SumWeightX ()         const { return _sum_weight_x; }
	double SumWeightY ()         const { return _sum_weight_y; }
	double SumWeightZ ()         const { return _sum_weight_z; }
	double Wx         (int iDiv) const { return _wx[iDiv];     }
	double Wy         (int iDiv) const { return _wy[iDiv];     }
	double Wz         (int iDiv) const { return _wz[iDiv];     }

	/* Find in which edge a vertex will be located, for given the indexes
	 * of the natural coordinates corresponding to each division.
	 *
	 * Ex.:  i=0 => r=-1    i=nDivX() => r=+1
	 *       j=0 => s=-1    j=nDivY() => s=+1
	 *       k=0 => t=-1    k=nDivZ() => t=+1
	 */
	void FindLocalEdgesFacesID (int i, int j, int k, Vertex * V) const;

	// Apply tags to Edges or Faces on boundary
	void ApplyTags (Elem * E) const;

	// Access the coordinates of all 8 or 20 nodes
	Matrix<double> const & C() const { return _c; }

	// Check if all arrays/vectors/matrices have the right sizes
	void Alright() const;

#ifdef USE_BOOST_PYTHON
	void PySet2D    (int Tag, BPy::list const & C, BPy::list const & Wx, BPy::list const & Wy);
	void PySet3D    (int               Tag,     ///< Tag to be replicated to elements
	                 BPy::list const & Verts,   ///< [(x1,y1,z1), (x2,y2,z2), ... 8 or 20 vertices]
	                 BPy::list const & Wx,      ///< Weights along x
	                 BPy::list const & Wy,      ///< Weights along y
	                 BPy::list const & Wz,      ///< Weights along z
	                 long              OrigID,  ///< ID of the vertex at origin
	                 long              XPlusID, ///< ID of the vertex at the positive corner of the local x-axis
	                 long              YPlusID, ///< ID of the vertex at the positive corner of the local y-axis
	                 long              ZPlusID, ///< ID of the vertex at the positive corner of the local z-axis
	                 BPy::list const & Edges,   ///< [(v1,v2), (v1,v2), ... 12 or 24 edges]
	                 BPy::dict const & ETags,   ///< {(v1,v2):tag1, (v3,v4):tag2, ... num edges with tags}
	                 BPy::dict const & FTags);  ///< {(v1,v2..v8):tag1, (v3,v4..v8):tag2, ... num faces with tags}
	void PySetETags (BPy::list const & Tags);
	void PySetFTags (BPy::list const & Tags);
#endif

private:
	Matrix<double> _c;            ///< X, Y, and Z coordinates of all 8 or 20 nodes of this block
	Array <double> _wx;           ///< X weights
	Array <double> _wy;           ///< Y weights
	Array <double> _wz;           ///< Z weights
	int            _n_div_x;      ///< number of divisions along X
	int            _n_div_y;      ///< number of divisions along Y
	int            _n_div_z;      ///< number of divisions along Z
	double         _sum_weight_x; ///< sum of weights along X
	double         _sum_weight_y; ///< sum of weights along Y
	double         _sum_weight_z; ///< sum of weights along Z
	bool           _is_3d;        ///< Is 3D block?
	int            _tag;          ///< A tag to be inherited by all elements generated inside this block
	bool           _has_etags;    ///< Has edge tags ?
	bool           _has_ftags;    ///< Has face tags ?
	Vector<int>    _etags;        ///< Edges tags: size = 2D:4, 3D:12 or 24
	Vector<int>    _ftags;        ///< Faces tags: size = 2D:0, 3D: 6

	void _set_wx();
	void _set_wy();
	void _set_wz();

}; // class Block

class Structured : public virtual Mesh::Generic
{
public:
	// Constants
	static Edge Edge2Vert[]; ///< Map from local edge ID to local vertex ID
	static Face Face2Vert[]; ///< Map from local face ID to local vertex IDs
	static Face Face2Edge[]; ///< Map from local face ID to local edge IDs

	// Constructor
	Structured (double Tol=sqrt(DBL_EPSILON)) : _tol(Tol) {} ///< Tol is the tolerance to regard two vertices as coincident

	// Destructor
	~Structured () { _erase(); }

	// Methods
	size_t Generate (Array<Block*> const & Blocks); ///< Returns the number of elements. Boundary marks are set first for Faces, then Edges, then Vertices (if any)

#ifdef USE_BOOST_PYTHON
	size_t PyGenerate (BPy::list const & ListOfMeshBlock);
#endif

private:
	// Data
	double         _tol;         ///< Tolerance to remove duplicate nodes
	Array<Vertex*> _verts_d;     ///< Vertices (with duplicates)
	Array<Vertex*> _verts_d_bry; ///< Vertices on boundary (with duplicates)
	Vector<double> _s;           ///< Current shape (interpolation) values, computed just after _shape(r,s,t)

	// Private methods
	void _shape_2d (double r, double s);
	void _shape_3d (double r, double s, double t);

	// Overloaded private methods
	void   _vtk_con          (Elem const * E, String & Connect) const;
	void   _erase            ();
	size_t _edge_to_lef_vert (size_t EdgeLocalID) const { return Edge2Vert[EdgeLocalID].L; }
	size_t _edge_to_rig_vert (size_t EdgeLocalID) const { return Edge2Vert[EdgeLocalID].R; }
	void   _face_to_verts    (size_t FaceLocalID, Array<size_t> & Verts) const;
	void   _face_to_edges    (size_t FaceLocalID, Array<size_t> & Edges) const;

}; // class Structured

Edge Structured::Edge2Vert[]= {{ 0, 3 },  // Edge #  0
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

Face Structured::Face2Vert[]= {{  0,  3,  7,  4, 11, 19, 15, 16 },  // Face # 0 => Vertices 0,3,7...
                               {  1,  2,  6,  5,  9, 18, 13, 17 },  // Face # 1
                               {  0,  1,  5,  4,  8, 17, 12, 16 },  // Face # 2
                               {  2,  3,  7,  6, 10, 19, 14, 18 },  // Face # 3
                               {  0,  1,  2,  3,  8,  9, 10, 11 },  // Face # 4
                               {  4,  5,  6,  7, 12, 13, 14, 15 }}; // Face # 5

Face Structured::Face2Edge[]= {{  0,  4,  8, 11, 12, 16, 20, 23 },  // Face # 0 => Edges 0,4,8
                               {  1,  5,  9, 10, 13, 17, 21, 22 },  // Face # 1
                               {  2,  6,  8,  9, 14, 18, 20, 21 },  // Face # 2
                               {  3,  7, 10, 11, 15, 19, 22, 23 },  // Face # 3
                               {  0,  1,  2,  3, 12, 13, 14, 15 },  // Face # 4
                               {  4,  5,  6,  7, 16, 17, 18, 19 }}; // Face # 5

Pair Block::Edge2Face[] = {{ 0, 4},  // Edge #  0 is shared among Faces # 0 and 4
                           { 1, 4},  // Edge #  1 ...
                           { 2, 4},  // Edge #  2
                           { 3, 4},  // Edge #  3
                           { 0, 5},  // Edge #  4
                           { 1, 5},  // Edge #  5
                           { 2, 5},  // Edge #  6
                           { 3, 5},  // Edge #  7
                           { 0, 2},  // Edge #  8
                           { 2, 1},  // Edge #  9
                           { 1, 3},  // Edge # 10
                           { 3, 0},  // Edge # 11
                           { 0, 4},  // Edge # 12
                           { 1, 4},  // Edge # 13
                           { 2, 4},  // Edge # 14
                           { 3, 4},  // Edge # 15
                           { 0, 5},  // Edge # 16
                           { 1, 5},  // Edge # 17
                           { 2, 5},  // Edge # 18
                           { 3, 5},  // Edge # 19
                           { 0, 2},  // Edge # 20
                           { 2, 1},  // Edge # 21
                           { 1, 3},  // Edge # 22
                           { 3, 0}}; // Edge # 23


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Methods ----------------------------------------------------------------------------------------------- Block

inline void Block::Set2D()
{
	_is_3d = false;
	_c.Resize        (2,8);
	_etags.Resize    (4);
	_etags.SetValues (0);
}

inline void Block::Set3D()
{
	_is_3d = true;
	_c    .Resize    (3,20);
	_etags.Resize    (24);
	_ftags.Resize    (6);
	_etags.SetValues (0);
	_ftags.SetValues (0);
}

inline void Block::SetWx(char const * Str)
{
	LineParser lp(Str);
	lp.ToArray (_wx);
	_set_wx ();
}

inline void Block::SetWy(char const * Str)
{
	LineParser lp(Str);
	lp.ToArray (_wy);
	_set_wy ();
}

inline void Block::SetWz(char const * Str)
{
	LineParser lp(Str);
	lp.ToArray (_wz);
	_set_wz ();
}

inline void Block::_set_wx()
{
	_n_div_x      = _wx.Size();
	_sum_weight_x = 0.0;
	for (int i=0; i<_n_div_x; ++i) _sum_weight_x += _wx[i];
	_wx.Push(0.0); // extra value just to allow looping over weights
}

inline void Block::_set_wy()
{
	_n_div_y      = _wy.Size();
	_sum_weight_y = 0.0;
	for (int i=0; i<_n_div_y; ++i) _sum_weight_y += _wy[i];
	_wy.Push(0.0); // extra value just to allow looping over weights
}

inline void Block::_set_wz()
{
	_n_div_z      = _wz.Size();
	_sum_weight_z = 0.0;
	for (int i=0; i<_n_div_z; ++i) _sum_weight_z += _wz[i];
	_wz.Push(0.0); // extra value just to allow looping over weights
}

inline void Block::FindLocalEdgesFacesID(int i, int j, int k, Vertex * V) const
{
	V->OnBry   = false;      // not on boundary
	V->EdgesID = -1, -1, -1; // not on any edge
	V->FacesID = -1, -1, -1; // not on any face
	if (_is_3d)
	{
		// Edges IDs
		int m = 0; // position in the EdgesID array
		if (k==0) // bottom
		{
			if (i==0)                   { V->EdgesID(m)= 0; V->OnBry=true; m++; } // behind
			if (i==_n_div_x)            { V->EdgesID(m)= 1; V->OnBry=true; m++; } // front
			if (j==0)                   { V->EdgesID(m)= 2; V->OnBry=true; m++; } // left
			if (j==_n_div_y)            { V->EdgesID(m)= 3; V->OnBry=true; m++; } // right
		}
		else if (k==_n_div_z) // top
		{
			if (i==0)                   { V->EdgesID(m)= 4; V->OnBry=true; m++; } // behind
			if (i==_n_div_x)            { V->EdgesID(m)= 5; V->OnBry=true; m++; } // front
			if (j==0)                   { V->EdgesID(m)= 6; V->OnBry=true; m++; } // left
			if (j==_n_div_y)            { V->EdgesID(m)= 7; V->OnBry=true; m++; } // right
		}
		if (i==0        && j==0       ) { V->EdgesID(m)= 8; V->OnBry=true; m++; } // vertical: left-behind
		if (i==_n_div_x && j==0       ) { V->EdgesID(m)= 9; V->OnBry=true; m++; } // vertical: left-front
		if (i==_n_div_x && j==_n_div_y) { V->EdgesID(m)=10; V->OnBry=true; m++; } // vertical: right-front
		if (i==0        && j==_n_div_y) { V->EdgesID(m)=11; V->OnBry=true; m++; } // vertical: right-behind

		// Faces IDs
		int n = 0; // position in the FacesID array
		if (i==0)        { V->FacesID(n)=0; V->OnBry=true; n++; } // behind
		if (i==_n_div_x) { V->FacesID(n)=1; V->OnBry=true; n++; } // front
		if (j==0)        { V->FacesID(n)=2; V->OnBry=true; n++; } // left
		if (j==_n_div_y) { V->FacesID(n)=3; V->OnBry=true; n++; } // right
		if (k==0)        { V->FacesID(n)=4; V->OnBry=true; n++; } // bottom
		if (k==_n_div_z) { V->FacesID(n)=5; V->OnBry=true; n++; } // top
	}
	else
	{
		int m = 0; // position in the EdgesID array
		if (i==0)        { V->EdgesID(m)=0; V->OnBry=true; m++; } // left
		if (i==_n_div_x) { V->EdgesID(m)=1; V->OnBry=true; m++; } // right
		if (j==0)        { V->EdgesID(m)=2; V->OnBry=true; m++; } // bottom
		if (j==_n_div_y) { V->EdgesID(m)=3; V->OnBry=true; m++; } // top
	}
}

inline void Block::ApplyTags(Elem * E) const
{
	if (_is_3d)
	{
		// ETags
		E->ETags.Resize(12);
		if (_has_etags)
		{
			for (size_t i=0; i<E->V.Size(); ++i)
			{
				if (E->V[i]->EdgesID(0)>-1) E->ETags(E->V[i]->EdgesID(0)) = _etags(E->V[i]->EdgesID(0));
				if (E->V[i]->EdgesID(1)>-1) E->ETags(E->V[i]->EdgesID(1)) = _etags(E->V[i]->EdgesID(1));
				if (E->V[i]->EdgesID(2)>-1) E->ETags(E->V[i]->EdgesID(2)) = _etags(E->V[i]->EdgesID(2));
			}
		}
		else E->ETags = 0;

		// FTags
		E->FTags.Resize(6);
		if (_has_ftags==false) { E->FTags = 0; return; }
		for (size_t i=0; i<E->V.Size(); ++i)
		{
			if (E->V[i]->FacesID(0)>-1) E->FTags(E->V[i]->FacesID(0)) = _ftags(E->V[i]->FacesID(0));
			if (E->V[i]->FacesID(1)>-1) E->FTags(E->V[i]->FacesID(1)) = _ftags(E->V[i]->FacesID(1));
			if (E->V[i]->FacesID(2)>-1) E->FTags(E->V[i]->FacesID(2)) = _ftags(E->V[i]->FacesID(2));
		}
	}
	else
	{
		// ETags
		E->ETags.Resize(4);
		if (_has_etags==false) { E->ETags = 0; return; }
		for (size_t i=0; i<E->V.Size(); ++i)
		{
			if (E->V[i]->EdgesID(0)>-1) E->ETags(E->V[i]->EdgesID(0)) = _etags(E->V[i]->EdgesID(0));
			if (E->V[i]->EdgesID(1)>-1) E->ETags(E->V[i]->EdgesID(1)) = _etags(E->V[i]->EdgesID(1));
			if (E->V[i]->EdgesID(2)>-1) E->ETags(E->V[i]->EdgesID(2)) = _etags(E->V[i]->EdgesID(2));
		}
	}
}

inline void Block::Alright() const
{
	if (_is_3d)
	{
		if (_c.Rows()!= 3) throw new Fatal("Block::Alright failed: C matrix of 3D blocks must have 3 rows (C.Rows==%d is invalid)",    _c.Rows());
		if (_c.Cols()!=20) throw new Fatal("Block::Alright failed: C matrix of 3D blocks must have 20 columns (C.Cols==%d is invalid)",_c.Cols());
	}
	else
	{
		if (_c.Rows()!=2) throw new Fatal("Block::Alright failed: C matrix of 2D blocks must have 2 rows (C.Rows==%d is invalid)",   _c.Rows());
		if (_c.Cols()!=8) throw new Fatal("Block::Alright failed: C matrix of 2D blocks must have 8 columns (C.Cols==%d is invalid)",_c.Cols());
	}
}


#ifdef USE_BOOST_PYTHON

inline void Block::PySet2D(int Tag, BPy::list const & C, BPy::list const & Wx, BPy::list const & Wy)
{
	// Basic information
	SetTag (Tag);
	Set2D  ();

	// Read C
	int nrow = BPy::len(C); if (nrow<1) throw new Fatal("Block::PySet2D: Number of rows of C matrix must be greater than 0 (%d is invalid)",nrow);
	int ncol = BPy::len(C[0]);
	for (int i=0; i<nrow; ++i)
	{
		BPy::list row = BPy::extract<BPy::list>(C[i])();
		if (BPy::len(row)!=ncol) throw new Fatal("Block::PySet2D: All rows of C matrix must have the same number of columns (%d is invalid)",BPy::len(row));
		for (int j=0; j<ncol; ++j) _c(i,j) = BPy::extract<double>(row[j])();
	}

	// Read Wx
	int sz_wx = BPy::len(Wx); if (sz_wx<1) throw new Fatal("Block::PySet2D: Number of elements in Wx list must be greater than 0 (%d is invalid)",sz_wx);
	_wx.Resize (sz_wx);
	for (int i=0; i<sz_wx; ++i) _wx[i] = BPy::extract<double>(Wx[i])();
	_set_wx();

	// Read Wy
	int sz_wy = BPy::len(Wy); if (sz_wy<1) throw new Fatal("Block::PySet2D: Number of elements in Wy list must be greater than 0 (%d is invalid)",sz_wy);
	_wy.Resize (sz_wy);
	for (int i=0; i<sz_wy; ++i) _wy[i] = BPy::extract<double>(Wy[i])();
	_set_wy();
}

inline void Block::PySet3D(int               Tag,    
                           BPy::list const & Verts,  
                           BPy::list const & Wx,     
                           BPy::list const & Wy,     
                           BPy::list const & Wz,     
                           long              OrigID, 
                           long              XPlusID,
                           long              YPlusID,
                           long              ZPlusID,
                           BPy::list const & Edges,  
                           BPy::dict const & ETags,  
                           BPy::dict const & FTags)  
{
	// Basic information
	SetTag (Tag);
	Set3D  ();

	// Read Wx
	int sz_wx = BPy::len(Wx); if (sz_wx<1) throw new Fatal("Block::PySet3D: Number of elements in Wx list must be greater than 0 (%d is invalid)",sz_wx);
	_wx.Resize (sz_wx);
	for (int i=0; i<sz_wx; ++i) _wx[i] = BPy::extract<double>(Wx[i])();
	_set_wx();

	// Read Wy
	int sz_wy = BPy::len(Wy); if (sz_wy<1) throw new Fatal("Block::PySet3D: Number of elements in Wy list must be greater than 0 (%d is invalid)",sz_wy);
	_wy.Resize (sz_wy);
	for (int i=0; i<sz_wy; ++i) _wy[i] = BPy::extract<double>(Wy[i])();
	_set_wy();

	// Read Wz
	int sz_wz = BPy::len(Wz); if (sz_wz<1) throw new Fatal("Block::PySet3D: Number of elements in Wz list must be greater than 0 (%d is invalid)",sz_wz);
	_wz.Resize (sz_wz);
	for (int i=0; i<sz_wz; ++i) _wz[i] = BPy::extract<double>(Wz[i])();
	_set_wz();

	// Check
	if (len(Edges)!=24) throw new Fatal("PySortBlock3D:: List with edges must have 24 items. Ex.: Edges = {(v1,v2):tag1, (v1,v2):tag2, ... 24 edges}. (%d is invalid)",len(Edges));

	// Edges
	Array<long> edges(48); // 48 => we have to serialize for Util::Tree
	int k = 0;
	for (int i=0; i<24; ++i)
	{
		BPy::tuple ed = BPy::extract<BPy::tuple>(Edges[i])();
		edges[k  ] = BPy::extract<long>(ed[0])();
		edges[k+1] = BPy::extract<long>(ed[1])();
		k += 2;
	}

	// Tree
	Util::Tree tree(edges);
	Array<long> eds;   eds.Resize(24);  eds.SetNS(Util::_4);
	Array<long> vg2l; vg2l.Resize(20); vg2l.SetNS(Util::_4); // map: global vertex to local vertex IDs
	Array<long> vl2g; vl2g.Resize(20); vl2g.SetNS(Util::_4); // map: local vertex to global vertex IDs
	vl2g[0] = OrigID;

	// Bottom nodes
	Array<long> path; path.SetNS(Util::_4);
	tree.DelEdge   (OrigID, XPlusID);
	tree.ShortPath (XPlusID, YPlusID, path);
	vl2g[ 1] = path[1];   vg2l[path[1]] =  1;
	vl2g[ 2] = path[3];   vg2l[path[3]] =  2;
	vl2g[ 3] = path[5];   vg2l[path[5]] =  3;
	vl2g[ 8] = path[0];   vg2l[path[0]] =  8;
	vl2g[ 9] = path[2];   vg2l[path[2]] =  9;
	vl2g[10] = path[4];   vg2l[path[4]] = 10;
	vl2g[11] = path[6];   vg2l[path[6]] = 11;

	// Behind nodes
	tree.DelEdge   (OrigID, YPlusID);
	tree.ShortPath (YPlusID, ZPlusID, path);
	vl2g[19] = path[2];   vg2l[path[2]] = 19;
	vl2g[ 7] = path[3];   vg2l[path[3]] =  7;
	vl2g[15] = path[4];   vg2l[path[4]] = 15;
	vl2g[ 4] = path[5];   vg2l[path[5]] =  4;
	vl2g[16] = path[6];   vg2l[path[6]] = 16;

	// Left nodes
	tree.DelEdge   (OrigID, ZPlusID);
	tree.ShortPath (ZPlusID, XPlusID, path);
	vl2g[12] = path[2];   vg2l[path[2]] = 12;
	vl2g[ 5] = path[3];   vg2l[path[3]] =  5;
	vl2g[17] = path[4];   vg2l[path[4]] = 17;

	// Corner-front nodes
	tree.DelEdge   (vl2g[7], vl2g[15]);
	tree.DelEdge   (vl2g[7], vl2g[19]);
	tree.ShortPath (vl2g[7], vl2g[ 5], path);
	vl2g[14] = path[1];   vg2l[path[1]] = 14;
	vl2g[ 6] = path[2];   vg2l[path[2]] =  6;
	vl2g[13] = path[3];   vg2l[path[3]] = 13;
	tree.ShortPath (vl2g[7], vl2g[2], path);
	vl2g[18] = path[3];   vg2l[path[3]] = 18;

	// 
	int net = BPy::len(ETags); // number of edges with tags
	int nft = BPy::len(FTags); // number of faces with tags
	Array<long> eg2l; eg2l.SetNS(Util::_4); // map: global edge to local edge IDs
	if (net>0 || nft>0)
	{
		eg2l.Resize(24);
		tree.Reset (edges);
		eg2l[tree.GetEdge(vl2g[0],vl2g[11])] =  0;
		eg2l[tree.GetEdge(vl2g[1],vl2g[ 9])] =  1;
		eg2l[tree.GetEdge(vl2g[0],vl2g[ 8])] =  2;
		eg2l[tree.GetEdge(vl2g[3],vl2g[10])] =  3;
		eg2l[tree.GetEdge(vl2g[4],vl2g[15])] =  4;
		eg2l[tree.GetEdge(vl2g[5],vl2g[13])] =  5;
		eg2l[tree.GetEdge(vl2g[4],vl2g[12])] =  6;
		eg2l[tree.GetEdge(vl2g[7],vl2g[14])] =  7;
		eg2l[tree.GetEdge(vl2g[0],vl2g[16])] =  8;
		eg2l[tree.GetEdge(vl2g[1],vl2g[17])] =  9;
		eg2l[tree.GetEdge(vl2g[2],vl2g[18])] = 10;
		eg2l[tree.GetEdge(vl2g[3],vl2g[19])] = 11;
		eg2l[tree.GetEdge(vl2g[3],vl2g[11])] = 12;
		eg2l[tree.GetEdge(vl2g[2],vl2g[ 9])] = 13;
		eg2l[tree.GetEdge(vl2g[1],vl2g[ 8])] = 14;
		eg2l[tree.GetEdge(vl2g[2],vl2g[10])] = 15;
		eg2l[tree.GetEdge(vl2g[7],vl2g[15])] = 16;
		eg2l[tree.GetEdge(vl2g[6],vl2g[13])] = 17;
		eg2l[tree.GetEdge(vl2g[5],vl2g[12])] = 18;
		eg2l[tree.GetEdge(vl2g[6],vl2g[14])] = 19;
		eg2l[tree.GetEdge(vl2g[4],vl2g[16])] = 20;
		eg2l[tree.GetEdge(vl2g[5],vl2g[17])] = 21;
		eg2l[tree.GetEdge(vl2g[6],vl2g[18])] = 22;
		eg2l[tree.GetEdge(vl2g[7],vl2g[19])] = 23;
	}

	// Read coordinates
	int nverts = BPy::len(Verts);
	if (nverts!=20) throw new Fatal("Block::PySet3D: Number of Verts must be either 8 or 20 (%d is invalid)",nverts);
	for (int i=0; i<20; ++i)
	{
		BPy::tuple xyz = BPy::extract<BPy::tuple>(Verts[vl2g[i]])();
		_c(0,i) = BPy::extract<double>(xyz[0])();
		_c(1,i) = BPy::extract<double>(xyz[1])();
		_c(2,i) = BPy::extract<double>(xyz[2])();
	}

	// Edge tags
	_has_etags = (net>0 ? true : false);
	for (int i=0; i<net; ++i)
	{
		BPy::tuple ed = BPy::extract<BPy::tuple>(ETags.keys()[i])();
		size_t eid = tree.GetEdge (BPy::extract<long>(ed[0])(), BPy::extract<long>(ed[1])());
		_etags(eg2l[eid]) = BPy::extract<int>(ETags.values()[i])();
	}

	// Face tags
	_has_ftags = (nft>0 ? true : false);
	for (int i=0; i<nft; ++i)
	{
		BPy::tuple fa = BPy::extract<BPy::tuple>(FTags.keys()[i])();
		long ed0 = eg2l[BPy::extract<long>(fa[0])()];
		long ed1 = eg2l[BPy::extract<long>(fa[1])()];
		long ed2 = eg2l[BPy::extract<long>(fa[2])()];
		int  sum[6] = { 0,0,0,0,0,0 };
		sum[Edge2Face[ed0].L] += 1;   sum[Edge2Face[ed0].R] += 1;
		sum[Edge2Face[ed1].L] += 1;   sum[Edge2Face[ed1].R] += 1;
		sum[Edge2Face[ed2].L] += 1;   sum[Edge2Face[ed2].R] += 1;
		for (int j=0; j<6; ++j)
		{
			if (sum[j]==3) // found
			{
				_ftags(j) = BPy::extract<int>(FTags.values()[i])();
				break;
			}
		}
	}
}

inline void Block::PySetETags(BPy::list const & Tags)
{
	int n_edges = BPy::len(Tags); // 2D => 4,  3D => 12
	if (_is_3d)
	{
		if (n_edges==12)
		{
			for (int i=0; i<n_edges; ++i)
				_etags(i) = BPy::extract<int>(Tags[i])();
			_has_etags = true;
		}
		else
		{
			if (n_edges!=24) throw new Fatal("Block::PySetETags: For 3D meshes, the number of edges must be 24 (n_edges==%d is invalid)",n_edges);
			for (int i=0; i<12; ++i)
				_etags(i) = BPy::extract<int>(Tags[i*2])();
			_has_etags = true;
		}
		/*
		if (_is_o2)
		{
			if (n_edges==24) throw new Fatal("Block::PySetETags: For 3D meshes, the number of edges must be 24 (n_edges==%d is invalid)",n_edges);
		}
		else
		{
			if (n_edges!=12) throw new Fatal("Block::PySetETags: For 3D meshes, the number of edges must be 12 (n_edges==%d is invalid)",n_edges);
		}
		*/
	}
	else
	{
		if (n_edges!= 4) throw new Fatal("Block::PySetETags: For 2D meshes, the number of edges must be 4 (n_edges==%d is invalid)",n_edges);
		for (int i=0; i<n_edges; ++i)
			_etags(i) = BPy::extract<int>(Tags[i])();
		_has_etags = true;
	}
}

inline void Block::PySetFTags(BPy::list const & Tags)
{
	if (_is_3d==false) throw new Fatal("Block::PySetFTags: This block is not 3D");
	int n_faces = BPy::len(Tags); // 3D => 6
	if (n_faces!=6) throw new Fatal("Block::PySetFTags: For 3D meshes, the number of faces must be 6 (n_faces==%d is invalid)",n_faces);
	for (int i=0; i<n_faces; ++i)
		_ftags(i) = BPy::extract<int>(Tags[i])();
	_has_ftags = true;
}

#endif // USE_BOOST_PYTHON


// Methods ------------------------------------------------------------------------------------------- Structured

inline size_t Structured::Generate(Array<Block*> const & Blocks)
{
	// Check
	if (Blocks.Size()<1) throw new Fatal("Structured::Generate: Number of blocks must be greater than 0 (%d is invalid)",Blocks.Size());

	// Erase previous mesh
	_erase();

	// Check if the first block is 3D
	_is_3d = Blocks[0]->Is3D();
	_s.Resize((_is_3d ? 20 : 8)); // resize the shape values vector

	// Generate vertices and elements (with duplicates)
	for (size_t b=0; b<Blocks.Size(); ++b)
	{
		// Check if all blocks have the same space dimension
		if (Blocks[b]->Is3D()!=_is_3d) throw new Fatal("Structured::Generate: All blocks must have the same space dimension");
		
		// Check if block has all arrays/vectors/matrices with the right sizes
		Blocks[b]->Alright();

		// Generate
		double t = -1.0; // initial Z natural coordinate
		for (int k=0; k<(_is_3d ? Blocks[b]->nDivZ()+1 : 1); ++k)
		{
			double s = -1.0; // initial Y natural coordinate
			for (int j=0; j<Blocks[b]->nDivY()+1; ++j)
			{
				double r = -1.0; // initial X natural coordinate
				for (int i=0; i<Blocks[b]->nDivX()+1; ++i)
				{
					// Compute shape (interpolation) functions
					if (_is_3d) _shape_3d (r,s,t);
					else        _shape_2d (r,s);

					// New vertex
					Vertex * v = new Vertex;
					v->MyID = _verts_d.Size();                  // id
					v->C    = Blocks[b]->C() * _s;              // new x-y-z coordinates
					v->Dupl = false;                            // is this a duplicated node?
					Blocks[b]->FindLocalEdgesFacesID(i,j,k, v); // check if it is on boundary and find EdgesID where this vertex is located on
					_verts_d.Push(v);
					if (v->OnBry) _verts_d_bry.Push(v); // array with vertices on boundary

					// New element
					if (i!=0 && j!=0 && (_is_3d ? k!=0 : true))
					{
						Elem * e = new Elem;
						e->MyID = _elems.Size();    // id
						e->Tag  = Blocks[b]->Tag(); // tag
						if (_is_3d)
						{
							// connectivity
							e->V.Resize(8);
							e->V[0] = _verts_d[v->MyID - 1 - (Blocks[b]->nDivX()+1) - (Blocks[b]->nDivX()+1)*(Blocks[b]->nDivY()+1)];
							e->V[1] = _verts_d[v->MyID     - (Blocks[b]->nDivX()+1) - (Blocks[b]->nDivX()+1)*(Blocks[b]->nDivY()+1)];
							e->V[2] = _verts_d[v->MyID                              - (Blocks[b]->nDivX()+1)*(Blocks[b]->nDivY()+1)];
							e->V[3] = _verts_d[v->MyID - 1                          - (Blocks[b]->nDivX()+1)*(Blocks[b]->nDivY()+1)];
							e->V[4] = _verts_d[v->MyID - 1 - (Blocks[b]->nDivX()+1)];
							e->V[5] = _verts_d[v->MyID -     (Blocks[b]->nDivX()+1)];
							e->V[6] = _verts_d[v->MyID];
							e->V[7] = _verts_d[v->MyID - 1];
							// shares information
							for (size_t m=0; m<8; ++m)
							{
								Share s = {e,m}; // The node shares information will have this element and the local node index
								e->V[m]->Shares.Push(s);
							}
							// is on boundary ?
							e->OnBry = (e->V[0]->OnBry || e->V[1]->OnBry || e->V[2]->OnBry || e->V[3]->OnBry || e->V[4]->OnBry || e->V[5]->OnBry || e->V[6]->OnBry || e->V[7]->OnBry); // if any node is on Bry, yes
							// VTK cell type
							e->VTKCellType = VTK_HEXAHEDRON;
						}
						else
						{
							// connectivity
							e->V.Resize(4);
							e->V[0] = _verts_d[v->MyID - 1 - (Blocks[b]->nDivX()+1)];
							e->V[1] = _verts_d[v->MyID     - (Blocks[b]->nDivX()+1)];
							e->V[2] = _verts_d[v->MyID];
							e->V[3] = _verts_d[v->MyID - 1];
							// shares information
							for (size_t m=0; m<4; ++m)
							{
								Share s = {e,m};
								e->V[m]->Shares.Push(s);
							}
							// is on boundary ?
							e->OnBry = (e->V[0]->OnBry || e->V[1]->OnBry || e->V[2]->OnBry || e->V[3]->OnBry); // if any node is on Bry, yes
							// VTK cell type
							e->VTKCellType = VTK_QUAD;
						}
						if (e->OnBry)
						{
							_elems_bry.Push(e);      // array with elements on boundary
							Blocks[b]->ApplyTags(e); // apply tags to edges and faces of this element with boundary tags
						}
						_elems.Push(e); // array with all elements
					}
					// Next r
					r += (2.0/Blocks[b]->SumWeightX()) * Blocks[b]->Wx(i);
				}
				// Next s
				s += (2.0/Blocks[b]->SumWeightY()) * Blocks[b]->Wy(j);
			}
			// Next t
			t += (_is_3d ? (2.0/Blocks[b]->SumWeightZ()) * Blocks[b]->Wz(k) : 0.0);
		}
	}

	// Remove duplicates
	long ncomp = 0; // number of comparisons
	long ndupl = 0; // number of duplicates
	if (Blocks.Size()>1)
	{
		for (size_t i=0; i<_verts_d_bry.Size(); ++i)
		{
			for (size_t j=i+1; j<_verts_d_bry.Size(); ++j)
			{
				// check distance
				double dist = sqrt(          pow(_verts_d_bry[i]->C(0)-_verts_d_bry[j]->C(0),2.0)+
											 pow(_verts_d_bry[i]->C(1)-_verts_d_bry[j]->C(1),2.0)+
								   (_is_3d ? pow(_verts_d_bry[i]->C(2)-_verts_d_bry[j]->C(2),2.0) : 0.0));
				if (dist<_tol)
				{
					/* TODO: this is wrong, since corner nodes can be duplicated and are still on boundary
					// If this node is duplicated, than it is not on-boundary any longer
					_verts_d_bry[i]->OnBry = false;
					*/
					// Mark duplicated
					_verts_d_bry[j]->Dupl = true;
					// Chage elements' connectivities
					for (size_t k=0; k<_verts_d_bry[j]->Shares.Size(); ++k)
					{
						Elem * e = _verts_d_bry[j]->Shares[k].E;
						int    n = _verts_d_bry[j]->Shares[k].N;
						e->V[n] = _verts_d_bry[i];
					}
					ndupl++;
				}
				ncomp++;
			}
		}
	}

	// Set new array with non-duplicated vertices
	size_t k = 0;
	size_t m = 0;
	_verts    .Resize (_verts_d    .Size()-ndupl);
	_verts_bry.Resize (_verts_d_bry.Size()-ndupl);
	for (size_t i=0; i<_verts_d.Size(); ++i)
	{
		if (_verts_d[i]->Dupl==false)
		{
			_verts[k]       = _verts_d[i]; // copy pointer
			_verts[k]->MyID = k;           // new ID
			if (_verts[k]->OnBry)
			{
				_verts_bry[m] = _verts[k]; // copy pointer
				m++;
			}
			k++;
		}
	}

	std::cout << "number of comparisons = " << ncomp << ", number of duplicates = " << ndupl << std::endl;

	return _elems.Size();
}


#ifdef USE_BOOST_PYTHON

inline size_t Structured::PyGenerate(BPy::list const & ListOfMeshBlock)
{
	int nb = BPy::len(ListOfMeshBlock);
	if (nb<1) throw new Fatal("Structured::PyGenerate: Number of blocks must be greater than 0 (%d is invalid)",nb);
	Array<Mesh::Block*> blocks;
	blocks.Resize (nb);
	for (int i=0; i<nb; ++i)
		blocks[i] = BPy::extract<Mesh::Block*>(ListOfMeshBlock[i])();
	return Mesh::Structured::Generate (blocks);
}

#endif // USE_BOOST_PYTHON


// Private methods -- Structured

inline void Structured::_shape_2d(double r, double s)
{
	/*      3           6            2
	 *        @---------@----------@
	 *        |               (1,1)|
	 *        |       s ^          |
	 *        |         |          |
	 *        |         |          |
	 *      7 @         +----> r   @ 5
	 *        |       (0,0)        |
	 *        |                    |
	 *        |                    |
	 *        |(-1,-1)             |
	 *        @---------@----------@
	 *      0           4            1
	 */
	double rp1=1.0+r; double rm1=1.0-r;
	double sp1=1.0+s; double sm1=1.0-s;

	_s(0) = 0.25*rm1*sm1*(rm1+sm1-3.0);
	_s(1) = 0.25*rp1*sm1*(rp1+sm1-3.0);
	_s(2) = 0.25*rp1*sp1*(rp1+sp1-3.0);
	_s(3) = 0.25*rm1*sp1*(rm1+sp1-3.0);
	_s(4) = 0.50*sm1*(1.0-r*r);
	_s(5) = 0.50*rp1*(1.0-s*s);
	_s(6) = 0.50*sp1*(1.0-r*r);
	_s(7) = 0.50*rm1*(1.0-s*s);
}

inline void Structured::_shape_3d(double r, double s, double t)
{
	/*     t            
	 *     |             
	 *     +---s       4         15         7
	 *   ,'              @_______@________@
	 *  r              ,'|              ,'|     The origin is in the
	 *            12 @'  |         14 ,'  |     centre of the cube:
	 *             ,'    |16        ,@    |19      0 => (-1,-1)
	 *       5   ,'      @      6 ,'      @        6 => ( 1, 1)
	 *         @'=======@=======@'        |        
	 *         |      13 |      |         |        
	 *         |         |      |  11     |        
	 *      17 |       0 @______|_@_______@
	 *         @       ,'       @       ,' 3      
	 *         |   8 @'      18 |     ,'           
	 *         |   ,'           |   ,@ 10          
	 *         | ,'             | ,'               
	 *         @_______@________@'                 
	 *        1        9         2                 
	 */
	_s( 0) = 0.125*(1.0-r)  *(1.0-s)  *(1.0-t)  *(-r-s-t-2.0);
	_s( 1) = 0.125*(1.0+r)  *(1.0-s)  *(1.0-t)  *( r-s-t-2.0);
	_s( 2) = 0.125*(1.0+r)  *(1.0+s)  *(1.0-t)  *( r+s-t-2.0);
	_s( 3) = 0.125*(1.0-r)  *(1.0+s)  *(1.0-t)  *(-r+s-t-2.0);
	_s( 4) = 0.125*(1.0-r)  *(1.0-s)  *(1.0+t)  *(-r-s+t-2.0);
	_s( 5) = 0.125*(1.0+r)  *(1.0-s)  *(1.0+t)  *( r-s+t-2.0);
	_s( 6) = 0.125*(1.0+r)  *(1.0+s)  *(1.0+t)  *( r+s+t-2.0);
	_s( 7) = 0.125*(1.0-r)  *(1.0+s)  *(1.0+t)  *(-r+s+t-2.0);

	_s( 8) = 0.25 *(1.0-r*r)*(1.0-s)  *(1.0-t);
	_s( 9) = 0.25 *(1.0+r)  *(1.0-s*s)*(1.0-t);
	_s(10) = 0.25 *(1.0-r*r)*(1.0+s)  *(1.0-t);
	_s(11) = 0.25 *(1.0-r)  *(1.0-s*s)*(1.0-t);

	_s(12) = 0.25 *(1.0-r*r)*(1.0-s)  *(1.0+t);
	_s(13) = 0.25 *(1.0+r)  *(1.0-s*s)*(1.0+t);
	_s(14) = 0.25 *(1.0-r*r)*(1.0+s)  *(1.0+t);
	_s(15) = 0.25 *(1.0-r)  *(1.0-s*s)*(1.0+t);

	_s(16) = 0.25 *(1.0-r)  *(1.0-s)  *(1.0-t*t);
	_s(17) = 0.25 *(1.0+r)  *(1.0-s)  *(1.0-t*t);
	_s(18) = 0.25 *(1.0+r)  *(1.0+s)  *(1.0-t*t);
	_s(19) = 0.25 *(1.0-r)  *(1.0+s)  *(1.0-t*t);
}

inline void Structured::_vtk_con(Elem const * E, String & Connect) const
{
	if (_is_3d)
	{
		Connect.Printf("%d %d %d %d %d %d %d %d",E->V[1]->MyID,
		                                         E->V[2]->MyID,
		                                         E->V[3]->MyID,
		                                         E->V[0]->MyID,
		                                         E->V[5]->MyID,
		                                         E->V[6]->MyID,
		                                         E->V[7]->MyID,
		                                         E->V[4]->MyID);
	}
	else
	{
		Connect.Printf("%d %d %d %d",E->V[0]->MyID,
		                             E->V[1]->MyID,
		                             E->V[2]->MyID,
		                             E->V[3]->MyID);
	}
}

inline void Structured::_erase()
{
	for (size_t i=0; i<_verts_d.Size(); ++i) if (_verts_d[i]!=NULL) delete _verts_d[i]; // it is only necessary to delete nodes in _verts_d array
	for (size_t i=0; i<_elems.  Size(); ++i) if (_elems  [i]!=NULL) delete _elems  [i]; // it is only necessary to delete elems in _elems array
	_is_3d = false;
	_verts_d    .Resize(0);
	_verts_d_bry.Resize(0);
	_verts      .Resize(0);
	_elems      .Resize(0);
	_elems_bry  .Resize(0);
	_verts_bry  .Resize(0);
}

inline void Structured::_face_to_verts(size_t FaceLocalID, Array<size_t> & Verts) const
{
	if (_is_o2)
	{
		Verts.Resize(8);
		Verts[0] = Face2Vert[FaceLocalID].I0;
		Verts[1] = Face2Vert[FaceLocalID].I1;
		Verts[2] = Face2Vert[FaceLocalID].I2;
		Verts[3] = Face2Vert[FaceLocalID].I3;
		Verts[4] = Face2Vert[FaceLocalID].I4;
		Verts[5] = Face2Vert[FaceLocalID].I5;
		Verts[6] = Face2Vert[FaceLocalID].I6;
		Verts[7] = Face2Vert[FaceLocalID].I7;
	}
	else
	{
		Verts.Resize(4);
		Verts[0] = Face2Vert[FaceLocalID].I0;
		Verts[1] = Face2Vert[FaceLocalID].I1;
		Verts[2] = Face2Vert[FaceLocalID].I2;
		Verts[3] = Face2Vert[FaceLocalID].I3;
	}
}

inline void Structured::_face_to_edges(size_t FaceLocalID, Array<size_t> & Edges) const
{
	if (_is_o2)
	{
		Edges.Resize(8);
		Edges[0] = Face2Edge[FaceLocalID].I0;
		Edges[1] = Face2Edge[FaceLocalID].I1;
		Edges[2] = Face2Edge[FaceLocalID].I2;
		Edges[3] = Face2Edge[FaceLocalID].I3;
		Edges[4] = Face2Edge[FaceLocalID].I4;
		Edges[5] = Face2Edge[FaceLocalID].I5;
		Edges[6] = Face2Edge[FaceLocalID].I6;
		Edges[7] = Face2Edge[FaceLocalID].I7;
	}
	else
	{
		Edges.Resize(4);
		Edges[0] = Face2Edge[FaceLocalID].I0;
		Edges[1] = Face2Edge[FaceLocalID].I1;
		Edges[2] = Face2Edge[FaceLocalID].I2;
		Edges[3] = Face2Edge[FaceLocalID].I3;
	}
}


}; // namespace Mesh

#endif // MECHSYS_MESH_STRUCTURED_H
