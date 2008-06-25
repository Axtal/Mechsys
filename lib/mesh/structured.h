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
#include <cfloat>

// Blitz++
#include <blitz/tinyvec-et.h>

// MechSys
#include "util/array.h"
#include "util/exception.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/laexpr.h"

#ifndef VTU_NEWLINE_DEFINED
  #define VTU_NEWLINE_DEFINED
  #define VTU_NEWLINE(I,K,N,KMAX,OF) if (K>KMAX) { OF<<(I<N-1?"\n        ":"\n"); K=0; } else if (I==N-1) { OF<<"\n"; }
#endif

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


struct Edge
{
	int L; // Left vertex id
	int R; // Right vertex id
};

const Edge EDGE2VERT[]= // Map between local edge ID to local vertex ID
{{ 0, 3 },
 { 1, 2 },
 { 0, 1 },
 { 2, 3 },
 { 4, 7 },
 { 5, 6 },
 { 4, 5 },
 { 6, 7 },
 { 0, 4 },
 { 1, 5 },
 { 2, 6 },
 { 3, 7 }};

struct Elem;

struct Share
{
	Elem* E; ///< The element
	int   N; ///< Local node index: 2D=>0,1,2,3, 3D=>0,1,2,3,4,5,6,7
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
	long           MyID;  ///< ID
	bool           OnBry; ///< On boundary?
	Array<Vertex*> V;     ///< Connectivity
	Vector<int>    ETags; ///< Edge tags (size==nLocalEdges)
	Vector<int>    FTags; ///< Face tags (size==nLocalFaces)
};

class Block
{
public:
	// Constructor
	Block () : _n_div_x(0), _n_div_y(0), _n_div_z(0), _e_tags(NULL), _f_tags(NULL) {}

	// Methods
	void Set (Matrix<double> * C, Array<double> * Wx, Array<double> * Wy, Array<double> * Wz=NULL); ///< C=coordinates, W=weights
	/* 
	 * 2D: 8 nodes => C. Resize(2, 8)
	 *                Wx.Resize(nDivX)
	 *                Wy.Resize(nDivY)
	 *          _                         _
	 *    C  = |  x0 x1 x2 x3 x4 x5 x6 x7  |
	 *         |_ y0 y1 y2 y3 y4 y5 y6 y7 _|
	 *
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
	void SetETags (Vector<int> * Tags) { _e_tags = Tags; } ///< Set edges tags: size == 2D:4, 3D:12
	void SetFTags (Vector<int> * Tags) { _f_tags = Tags; } ///< Set faces tags: size == 2D:0, 3D: 6

	// Access methods
	bool   Is3D       ()         const { return _is_3d;        }
	int    nDivX      ()         const { return _n_div_x;      }
	int    nDivY      ()         const { return _n_div_y;      }
	int    nDivZ      ()         const { return _n_div_z;      }
	double SumWeightX ()         const { return _sum_weight_x; }
	double SumWeightY ()         const { return _sum_weight_y; }
	double SumWeightZ ()         const { return _sum_weight_z; }
	double Wx         (int iDiv) const { return (*_wx)[iDiv];  }
	double Wy         (int iDiv) const { return (*_wy)[iDiv];  }
	double Wz         (int iDiv) const { return (*_wz)[iDiv];  }

	/* Find in which edge a vertex will be located, for given the indexes
	 * of the natural coordinates corresponding to each division.
	 *
	 * Ex.:  i=0 => r=-1    i=nDivX() => r=+1
	 *       j=0 => s=-1    j=nDivY() => s=+1
	 *       k=0 => t=-1    k=nDivZ() => t=+1
	 */
	void FindLocalEdgesFacesID (int i, int j, int k, Vertex * V) const;

	// Apply tags to Edges or Faces on boundary
	void Tag (Elem * E) const;

	// Access the coordinates of all 8 or 20 nodes
	Matrix<double> const & C() const { return (*_c); }

	// Check if all arrays/vectors/matrices have the right sizes
	void Alright() const;

private:
	Matrix<double> * _c;            ///< X, Y, and Z coordinates of all 8 or 20 nodes of this block
	Array <double> * _wx;           ///< X weights
	Array <double> * _wy;           ///< Y weights
	Array <double> * _wz;           ///< Z weights
	int              _n_div_x;      ///< number of divisions along X
	int              _n_div_y;      ///< number of divisions along Y
	int              _n_div_z;      ///< number of divisions along Z
	double           _sum_weight_x; ///< sum of weights along X
	double           _sum_weight_y; ///< sum of weights along Y
	double           _sum_weight_z; ///< sum of weights along Z
	bool             _is_3d;        ///< Is 3D block?
	Vector<int>    * _e_tags;       ///< Edges tags: size = 2D:4, 3D:12
	Vector<int>    * _f_tags;       ///< Faces tags: size = 2D:0, 3D: 6

}; // class Block

class Structured
{
public:
	// Constructor
	Structured (double Tol=sqrt(DBL_EPSILON)) : _tol(Tol), _is_3d(false) {}

	// Destructor
	~Structured ();

	// Methods
	size_t Generate (Array<Block*> const & Blocks); ///< Returns the number of elements. Boundary marks are set first for Faces, then Edges, then Vertices (if any)
	void   WriteVTU (char const * FileName) const;

	// Access methods
	bool                   Is3D     () const { return _is_3d;     }
	Array<Vertex*> const & Verts    () const { return _verts;     }
	Array<Elem*>   const & Elems    () const { return _elems;     }
	Array<Elem*>   const & ElemsBry () const { return _elems_bry; } ///< Elements on boundary

private:
	// Data
	double         _tol;         ///< Tolerance to remove duplicate nodes
	bool           _is_3d;       ///< Is 3D mesh?
	Array<Vertex*> _verts_d;     ///< Vertices (with duplicates)
	Array<Vertex*> _verts_d_bry; ///< Vertices on boundary (with duplicates)
	Array<Vertex*> _verts;       ///< Vertices
	Array<Elem*>   _elems;       ///< Elements
	Array<Elem*>   _elems_bry;   ///< Elements on boundary
	Vector<double> _s;           ///< Current shape (interpolation) values, computed just after _shape(r,s,t)

	// Private methods
	void _shape_2d (double r, double s);
	void _shape_3d (double r, double s, double t);
	void _vtk_con  (Elem const * E, String & Connect) const;

}; // class Structured


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Methods -- Block

inline void Block::Set(Matrix<double> * C, Array<double> * Wx, Array<double> * Wy, Array<double> * Wz)
{
	// Coordinates
	_c     = C;
	_is_3d = (C->Rows()>2 ? true : false);

	// 2D data
	_wx           = Wx;
	_wy           = Wy;
	_n_div_x      = _wx->Size();
	_n_div_y      = _wy->Size();
	_sum_weight_x = 0.0;
	_sum_weight_y = 0.0;
	for (int i=0; i<_n_div_x; ++i) _sum_weight_x += (*_wx)[i];
	for (int i=0; i<_n_div_y; ++i) _sum_weight_y += (*_wy)[i];
	_wx->Push(0.0); // extra value just to allow looping over weights
	_wy->Push(0.0);

	// 3D data
	if (_is_3d)
	{
		if (Wz==NULL) throw new Fatal("Block::Set: For 3D blocks, Wz must be given (not NULL)");
		_wz           = Wz;
		_n_div_z      = _wz->Size();
		_sum_weight_z = 0.0;
		for (int i=0; i<_n_div_z; ++i) _sum_weight_z += (*_wz)[i];
		_wz->Push(0.0);
	}
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

inline void Block::Tag(Elem * E) const
{
	if (_is_3d)
	{
		// ETags
		E->ETags.Resize(12);
		if (_e_tags!=NULL)
		{
			for (size_t i=0; i<E->V.Size(); ++i)
			{
				if (E->V[i]->EdgesID(0)>-1) E->ETags(E->V[i]->EdgesID(0)) = (*_e_tags)(E->V[i]->EdgesID(0));
				if (E->V[i]->EdgesID(1)>-1) E->ETags(E->V[i]->EdgesID(1)) = (*_e_tags)(E->V[i]->EdgesID(1));
				if (E->V[i]->EdgesID(2)>-1) E->ETags(E->V[i]->EdgesID(2)) = (*_e_tags)(E->V[i]->EdgesID(2));
			}
		}
		else E->ETags = 0;

		// FTags
		E->FTags.Resize(6);
		if (_f_tags==NULL) { E->FTags = 0; return; }
		for (size_t i=0; i<E->V.Size(); ++i)
		{
			if (E->V[i]->FacesID(0)>-1) E->FTags(E->V[i]->FacesID(0)) = (*_f_tags)(E->V[i]->FacesID(0));
			if (E->V[i]->FacesID(1)>-1) E->FTags(E->V[i]->FacesID(1)) = (*_f_tags)(E->V[i]->FacesID(1));
			if (E->V[i]->FacesID(2)>-1) E->FTags(E->V[i]->FacesID(2)) = (*_f_tags)(E->V[i]->FacesID(2));
		}
	}
	else
	{
		// ETags
		E->ETags.Resize(4);
		if (_e_tags==NULL) { E->ETags = 0; return; }
		for (size_t i=0; i<E->V.Size(); ++i)
		{
			if (E->V[i]->EdgesID(0)>-1) E->ETags(E->V[i]->EdgesID(0)) = (*_e_tags)(E->V[i]->EdgesID(0));
			if (E->V[i]->EdgesID(1)>-1) E->ETags(E->V[i]->EdgesID(1)) = (*_e_tags)(E->V[i]->EdgesID(1));
			if (E->V[i]->EdgesID(2)>-1) E->ETags(E->V[i]->EdgesID(2)) = (*_e_tags)(E->V[i]->EdgesID(2));
		}
	}
}

inline void Block::Alright() const
{
	if (_is_3d)
	{
		if (_c->Rows()!= 3) throw new Fatal("Block::Alright failed: C matrix of 3D blocks must have 3 rows (C.Rows==%d is invalid)",    _c->Rows());
		if (_c->Cols()!=20) throw new Fatal("Block::Alright failed: C matrix of 3D blocks must have 20 columns (C.Cols==%d is invalid)",_c->Cols());
		if (_e_tags!=NULL)
			if (_e_tags->Size()!=12) throw new Fatal("Block::Alright failed: ETags array of 3D blocks must have 12 elements (ETags.Size==%d is invalid)",_e_tags->Size());
		if (_f_tags!=NULL)
			if (_f_tags->Size()!= 6) throw new Fatal("Block::Alright failed: ETags array of 3D blocks must have 6 elements (ETags.Size==%d is invalid)",_f_tags->Size());
	}
	else
	{
		if (_c->Rows()!=2) throw new Fatal("Block::Alright failed: C matrix of 2D blocks must have 2 rows (C.Rows==%d is invalid)",   _c->Rows());
		if (_c->Cols()!=8) throw new Fatal("Block::Alright failed: C matrix of 2D blocks must have 8 columns (C.Cols==%d is invalid)",_c->Cols());
		if (_e_tags!=NULL)
			if (_e_tags->Size()!=4) throw new Fatal("Block::Alright failed: ETags array of 2D blocks must have 4 elements (ETags.Size==%d is invalid)",_e_tags->Size());
	}
}

// Destructor -- Structured

inline Structured::~Structured()
{
	for (size_t i=0; i<_verts_d.Size(); ++i) if (_verts_d[i]!=NULL) delete _verts_d[i]; // it is only necessary to delete nodes in _verts_d array
	for (size_t i=0; i<_elems.  Size(); ++i) if (_elems  [i]!=NULL) delete _elems  [i];
}

// Methods -- Structured

inline size_t Structured::Generate(Array<Block*> const & Blocks)
{
	// Check
	if (Blocks.Size()<1) throw new Fatal("Structured::Generate: Number of blocks must be greater than 0 (%d is invalid)",Blocks.Size());

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
						e->MyID = _elems.Size(); // id
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
						}
						if (e->OnBry)
						{
							_elems_bry.Push(e); // array with elements on boundary
							Blocks[b]->Tag(e);  // tag edges and faces of this element with boundary tags
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

	// Set new array with non-duplicated vertices
	size_t k = 0;
	_verts.Resize(_verts_d.Size()-ndupl);
	for (size_t i=0; i<_verts_d.Size(); ++i)
	{
		if (_verts_d[i]->Dupl==false)
		{
			_verts[k]       = _verts_d[i]; // copy pointer
			_verts[k]->MyID = k;           // new ID
			k++;
		}
	}

	std::cout << "number of comparisons = " << ncomp << ", number of duplicates = " << ndupl << std::endl;

	return _elems.Size();
}

inline void Structured::WriteVTU(char const * FileName) const
{
	// Open File
	std::ofstream      of(FileName, std::ios::out);
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
		//VTU_NEWLINE (i,k,ne,nimax/(_is_3d?20:8),oss);
		VTU_NEWLINE (i,k,ne,nimax/(_is_3d?8:4),oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	size_t ossfset = 0;
	for (size_t i=0; i<ne; ++i)
	{
		//ossfset += (_is_3d?20:8);
		ossfset += (_is_3d?8:4);
		oss << (k==0?"  ":" ") << ossfset;
		k++;
		VTU_NEWLINE (i,k,ne,nimax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		//oss << (k==0?"  ":" ") << (_is_3d?25:23); // VTK_QUADRATIC_HEXAHEDRON or VTK_QUADRATIC_QUAD
		oss << (k==0?"  ":" ") << (_is_3d?12:9); // VTK_HEXAHEDRON or VTK_QUAD
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
		VTU_NEWLINE (i,k,ne,nimax,oss);
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

	/* TODO:
	 * remove this, because, after the removal of duplicates,
	 * the shares information is not valid any longer
	 * I kept this here because it may be useful for debugging purposes
	 * (before the removal of duplicates)*/
	/*
	oss << "        <DataArray type=\"Float32\" Name=\"" << "shares" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<nn; ++i)
	{
		oss << (k==0?"  ":" ") << _verts[i]->Shares.Size();
		k++;
		VTU_NEWLINE (i,k,nn,nfmax,oss);
	}
	oss << "        </DataArray>\n";
	*/

	oss << "      </PointData>\n";

	// Data -- elements
	oss << "      <CellData Scalars=\"TheScalars\">\n";
	oss << "        <DataArray type=\"Float32\" Name=\"" << "onbry" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	k = 0; oss << "        ";
	for (size_t i=0; i<ne; ++i)
	{
		oss << (k==0?"  ":" ") << _elems[i]->OnBry;
		k++;
		VTU_NEWLINE (i,k,ne,nfmax,oss);
	}
	oss << "        </DataArray>\n";
	oss << "      </CellData>\n";

	// Bottom
	oss << "    </Piece>\n";
	oss << "  </UnstructuredGrid>\n";
	oss << "</VTKFile>" << std::endl;

	// Write to file
	of << oss.str();
	of.close();
}

// Private methods

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


}; // namespace Mesh


#ifdef USE_BOOST_PYTHON

namespace boopy = boost::python;

class PyMeshBlock
{
public:
	void Set (boopy::list const & C, boopy::list const & Wx, boopy::list const & Wy)
	{
		// Read C
		int nrow = boopy::len(C); if (nrow<1) throw new Fatal("PyMeshBlock: Number of rows of C matrix must be greater than 0 (%d is invalid)",nrow);
		int ncol = boopy::len(C[0]);
		_c.Resize (nrow,ncol);
		for (int i=0; i<nrow; ++i)
		{
			boopy::list row = boopy::extract<boopy::list>(C[i])();
			if (boopy::len(row)!=ncol) throw new Fatal("PyMeshBlock: All rows of C matrix must have the same number of columns (%d is invalid)",boopy::len(row));
			for (int j=0; j<ncol; ++j) _c(i,j) = boopy::extract<double>(row[j])();
		}

		// Read Wx
		int sz_wx = boopy::len(Wx); if (sz_wx<1) throw new Fatal("PyMeshBlock: Number of elements in Wx list must be greater than 0 (%d is invalid)",sz_wx);
		_wx.Resize (sz_wx);
		for (int i=0; i<sz_wx; ++i) _wx[i] = boopy::extract<double>(Wx[i])();

		// Read Wy
		int sz_wy = boopy::len(Wy); if (sz_wy<1) throw new Fatal("PyMeshBlock: Number of elements in Wy list must be greater than 0 (%d is invalid)",sz_wy);
		_wy.Resize (sz_wy);
		for (int i=0; i<sz_wy; ++i) _wy[i] = boopy::extract<double>(Wy[i])();

		// Set _block
		_block.Set (&_c, &_wx, &_wy);
	}

	void Set (boopy::list const & C, boopy::list const & Wx, boopy::list const & Wy, boopy::list const & Wz)
	{
		// Read C
		int nrow = boopy::len(C); if (nrow<1) throw new Fatal("PyMeshBlock: Number of rows of C matrix must be greater than 0 (%d is invalid)",nrow);
		int ncol = boopy::len(C[0]);
		_c.Resize (nrow,ncol);
		for (int i=0; i<nrow; ++i)
		{
			boopy::list row = boopy::extract<boopy::list>(C[i])();
			if (boopy::len(row)!=ncol) throw new Fatal("PyMeshBlock: All rows of C matrix must have the same number of columns (%d is invalid)",boopy::len(row));
			for (int j=0; j<ncol; ++j) _c(i,j) = boopy::extract<double>(row[j])();
		}

		// Read Wx
		int sz_wx = boopy::len(Wx); if (sz_wx<1) throw new Fatal("PyMeshBlock: Number of elements in Wx list must be greater than 0 (%d is invalid)",sz_wx);
		_wx.Resize (sz_wx);
		for (int i=0; i<sz_wx; ++i) _wx[i] = boopy::extract<double>(Wx[i])();

		// Read Wy
		int sz_wy = boopy::len(Wy); if (sz_wy<1) throw new Fatal("PyMeshBlock: Number of elements in Wy list must be greater than 0 (%d is invalid)",sz_wy);
		_wy.Resize (sz_wy);
		for (int i=0; i<sz_wy; ++i) _wy[i] = boopy::extract<double>(Wy[i])();

		// Read Wz
		int sz_wz = boopy::len(Wz); if (sz_wz<1) throw new Fatal("PyMeshBlock: Number of elements in Wz list must be greater than 0 (%d is invalid)",sz_wz);
		_wz.Resize (sz_wz);
		for (int i=0; i<sz_wz; ++i) _wz[i] = boopy::extract<double>(Wz[i])();

		// Set _block
		_block.Set (&_c, &_wx, &_wy, &_wz);
	}

	void SetETags (boopy::list const & Tags)
	{
		int n_edges = boopy::len(Tags); // 2D => 4,  3D => 12
		if (!(n_edges==4 || n_edges==12)) throw new Fatal("PyMeshBlock::SetETags: Number of edges must be 4 (2D) or 12 (3D) (n_edges==%d is invalid)",n_edges);
		_e_tags.Resize(n_edges);
		for (int i=0; i<n_edges; ++i) _e_tags(i) = boopy::extract<int>(Tags[i])();
		_block.SetETags (&_e_tags);
	}

	void SetFTags (boopy::list const & Tags)
	{
		int n_faces = boopy::len(Tags); // 3D => 6
		if (n_faces!=6) throw new Fatal("PyMeshBlock::SetETags: Number of faces must be 6 (n_faces==%d is invalid)",n_faces);
		_f_tags.Resize(n_faces);
		for (int i=0; i<n_faces; ++i) _f_tags(i) = boopy::extract<int>(Tags[i])();
		_block.SetFTags (&_f_tags);
	}

	Mesh::Block * GetBlock () { return &_block; }

private:
	Matrix<double> _c;
	Array<double> _wx;
	Array<double> _wy;
	Array<double> _wz;
	Vector<int>   _v_tags;
	Vector<int>   _e_tags;
	Vector<int>   _f_tags;
	Mesh::Block   _block;

}; // class PyMeshBlock

void (PyMeshBlock::*PMBSet1)(boopy::list const & C, boopy::list const & Wx, boopy::list const & Wy)                         = &PyMeshBlock::Set;
void (PyMeshBlock::*PMBSet2)(boopy::list const & C, boopy::list const & Wx, boopy::list const & Wy, boopy::list const & Wz) = &PyMeshBlock::Set;

class PyMeshStruct
{
public:
	PyMeshStruct ()           : _ms(sqrt(DBL_EPSILON)) {}
	PyMeshStruct (double Tol) : _ms(Tol)               {}

	size_t Generate (boopy::list & ListOfPyMeshBlock)
	{
		int nb = boopy::len(ListOfPyMeshBlock); if (nb<1) throw new Fatal("PyMeshStruct: Number of blocks must be greater than 0 (%d is invalid)",nb);
		Array<Mesh::Block*> blocks;
		blocks.Resize (nb);
		for (int i=0; i<nb; ++i)
		{
			PyMeshBlock * blk = boopy::extract<PyMeshBlock*>(ListOfPyMeshBlock[i])();
			blocks[i] = blk->GetBlock();
		}
		return _ms.Generate (blocks);
	}

	void WriteVTU (boopy::str FileName) { _ms.WriteVTU(boopy::extract<char const *>(FileName)()); }

	void GetVerts (boopy::list & V)
	{
		if (_ms.Is3D())
		{
			for (size_t i=0; i<_ms.Verts().Size(); ++i)
				V.append (boopy::make_tuple(_ms.Verts()[i]->C(0), _ms.Verts()[i]->C(1), _ms.Verts()[i]->C(2)));
		}
		else
		{
			for (size_t i=0; i<_ms.Verts().Size(); ++i)
				V.append (boopy::make_tuple(_ms.Verts()[i]->C(0), _ms.Verts()[i]->C(1), 0.0));
		}
	}

	void GetElems (boopy::list & E)
	{
		for (size_t i=0; i<_ms.Elems().Size(); ++i)
		{
			boopy::list conn;
			for (size_t j=0; j<_ms.Elems()[i]->V.Size(); ++j)
				conn.append (_ms.Elems()[i]->V[j]->MyID);
			E.append (conn);
		}
	}

	void GetETags (boopy::list & Tags)
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
		for (size_t i=0; i<_ms.ElemsBry().Size(); ++i) // elements on boundary
		{
			for (int j=0; j<_ms.ElemsBry()[i]->ETags.Size(); ++j) // j is the local_edge_id
			{
				int tag = _ms.ElemsBry()[i]->ETags(j);
				if (tag<0)
				{
					int eid = _ms.ElemsBry()[i]->MyID;
					int L   = _ms.ElemsBry()[i]->V[Mesh::EDGE2VERT[j].L]->MyID;
					int R   = _ms.ElemsBry()[i]->V[Mesh::EDGE2VERT[j].R]->MyID;
					Tags.append (boopy::make_tuple(eid, L, R, tag));
				}
			}
		}
	}

private:
	Mesh::Structured _ms;

}; // class PyMeshStruct

#endif


#endif // MECHSYS_MESH_STRUCTURED_H
