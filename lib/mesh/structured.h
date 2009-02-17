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

/* LOCAL indexes of VERTICES, EDGES, and FACES

  2D:
                Vertices                              Edges
   y                                                                                         
   |        3      6      2                         3       7
   +--x      @-----@-----@                        +-----+-----+
             |           |                        |           |
             |           |                      4 |           |5
           7 @           @ 5                      +           +
             |           |                        |           |
             |           |                       0|           |1
             @-----@-----@                        +-----+-----+
            0      4      1                         2       6

  3D:
                  Vertices                            Edges                              Faces
    z
    |           4        15        7                   4      16
   ,+--y         @-------@--------@                 +-------+--------+                 +----------------+ 
 x'            ,'|              ,'|              6,'|            7 ,'|               ,'|              ,'| 
          12 @'  |         14 ,'  |             +'  |            ,'  |23           ,'  |  ___       ,'  | 
           ,'    |16        ,@    |19      18 ,'  20|          ,+    |           ,'    |,'5,'  [0],'    | 
     5   ,'      @      6 ,'      @         ,'      +  17    ,19     +         ,'      |~~~     ,'      | 
       @'=======@=======@'        |       +'=======+=======+'        |       +'===============+'  ,'|   | 
       |      13 |      |         |       |    5    |      |         |11     |   ,'|   |      |   |3|   | 
       |         |      |  11     |     21|        8|  0   |22  12   |       |   |2|   |      |   |,'   | 
    17 |       0 @- - - | @- - - -@       |         +- - - | +- - - -+       |   |,'   +- - - | +- - - -+ 
       @       ,'       @       ,' 3      +      2,'       +       ,'        |       ,'       |       ,'  
       |   8 @'      18 |     ,'          |     +'         |     ,'3         |     ,' [1]  ___|     ,'    
       |   ,'           |   ,@ 10        9|   14         10|   ,+            |   ,'      ,'4,'|   ,'      
       | ,'             | ,'              | ,'             | ,15             | ,'        ~~~  | ,'        
       @-------@--------@'                +-------+--------+'                +----------------+'          
     1         9         2                    1       13
*/

// Std Lib
#include <iostream>
#include <fstream>
#include <cfloat>   // for DBL_EPSILON
#include <cstdarg>  // for va_list, va_start, va_end
#include <ctime>    /// for std::clock

// Blitz++
#include <blitz/tinyvec-et.h>

// Boost
#include <boost/tuple/tuple_io.hpp>

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

// typedefs
typedef Array< boost::tuple<int,double,double,double> > Verts_T; // ID,x,y,z
typedef Array< boost::tuple<int,int> >                  Edges_T; // v1,v2
typedef Array< boost::tuple<int,int,int> >              ETags_T; // v1,v2,tag
typedef Array< boost::tuple<int,int,int,int,int> >      FTags_T; // v1,v2,v3,v4,tag

/* Structured Block. */
class Block
{
public:
	/* 2D:           _             _
	 *     Coords = |  x0 x1 x2 x3  |
	 *              |_ y0 y1 y2 y3 _|
	 *     or       
	 *               _                          _
	 *     Coords = |  x0 x1 x2 x3  x4 x5 x6 x7  |
	 *              |_ y0 y1 y2 y3  y4 y5 y6 y7 _|
	 *
	 * 3D:           _                          _
	 *              |  x0 x1 x2 x3 x4 x5 x6 x7   |
	 *     Coords = |  y0 y1 y2 y3 y4 y5 y6 y7   |
	 *              |_ z0 z1 z2 z3 z4 z5 z6 z7  _|
	 *     
	 *     or        _                                         _
	 *              |  x0 x1 x2 x3 x4 x5 x6 x7 ... x17 x18 x19  |
	 *     Coords = |  y0 y1 y2 y3 y4 y5 y6 y7 ... y17 y18 y19  |
	 *              |_ z0 z1 z2 z3 z4 z5 z6 z7 ... z17 z18 z19 _|
	 */
	// Constants
	static Pair Edge2Face[]; ///< Map from local edge ID to a pair of faces (local ID) that share this edge

	// Constructor
	Block () : _n_div_x(0), _n_div_y(0), _n_div_z(0), _tag(-1), _has_etags(false), _has_ftags(false) {}

	// Set methods
	void Set (int             Tag,         ///< Tag to be replicated to elements
	          Verts_T const & Verts,       ///< Vertices
	          Edges_T const & Edges,       ///< Edges
	          ETags_T const * ETags,       ///< Edges Tags (can be NULL)
	          FTags_T const * FTags,       ///< Face Tags (can be NULL)
	          long            OrigID,      ///< ID of the vertex at origin
	          long            XPlusID,     ///< ID of the vertex at the positive corner of the local x-axis
	          long            YPlusID,     ///< ID of the vertex at the positive corner of the local y-axis
	          long            ZPlusID=-1); ///< ID of the vertex at the positive corner of the local z-axis

	void SetNx (size_t Nx, double Ax=0.0, bool NonLin=false); ///< Set number of x divisions with all x weights equal to 1.0
	void SetNy (size_t Ny, double Ay=0.0, bool NonLin=false); ///< Set number of y divisions with all y weights equal to 1.0
	void SetNz (size_t Nz, double Az=0.0, bool NonLin=false); ///< Set number of z divisions with all z weights equal to 1.0

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
	double Diagonal   ()         const;

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
	void PySet (int               Tag,      ///< Tag to be replicated to elements
	            BPy::list const & Verts,    ///< Vertices: [(ID,x,y), (ID,x,y), (ID,x,y) ...] or [(ID,x,y,z)...]
	            BPy::list const & Edges,    ///< Edges: [(v1,v2), (v1,v2), (v1,v2) ...]
	            BPy::list const & ETags,    ///< Edges Tags: [(v1,v2,tag), (v1,v2,tag) ...] (can be [])
	            BPy::list const & FTags,    ///< Face Tags: [(v1,v2,v3,v4,tag) ...] (can be [])
	            long              OrigID,   ///< ID of the vertex at origin
	            long              XPlusID,  ///< ID of the vertex at the positive corner of the local x-axis
	            long              YPlusID,  ///< ID of the vertex at the positive corner of the local y-axis
	            long              ZPlusID); ///< ID of the vertex at the positive corner of the local z-axis
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
	Vector<int>    _etags;        ///< Edges tags: size = 2D:8, 3D:24
	Vector<int>    _ftags;        ///< Faces tags: size = 2D:0, 3D:6

	void _set_2d();
	void _set_3d();
	void _set_wx();
	void _set_wy();
	void _set_wz();
	void _gen_mid_nodes();

}; // class Block

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


class Structured : public virtual Mesh::Generic
{
public:
	// Constructor
	Structured (bool Is3D) : Mesh::Generic(Is3D) {}

	// Destructor
	~Structured () { Erase(); }

	// Methods
	void   Erase     ();
	void   SetBlocks (Array<Block> const & Blocks) { _bls = Blocks; }
	size_t Generate  (bool WithInfo=false); ///< Returns the number of elements. Boundary marks are set first for Faces, then Edges, then Vertices (if any)

	void GenBox (int Nx, int Ny, int Nz, double Lx=1.0, double Ly=1.0, double Lz=1.0); ///< Generate a cube with dimensions Lx,Ly,Lz and with tags on faces

#ifdef USE_BOOST_PYTHON
	void   PySetBlocks (BPy::list const & ListOfMeshBlock);
#endif

private:
	// Data
	Array<Block>   _bls;         ///< Blocks defining the input geometry
	Array<Vertex*> _verts_d;     ///< Vertices (with duplicates)
	Array<Vertex*> _verts_d_bry; ///< Vertices on boundary (with duplicates)
	Array<Vertex*> _verts_m1;    ///< X O2 Vertices (with duplicates)
	Array<Vertex*> _verts_m2;    ///< Y O2 Vertices (with duplicates)
	Array<Vertex*> _verts_m3;    ///< Z O2 Vertices (with duplicates)
	Vector<double> _s;           ///< Current shape (interpolation) values, computed just after _shape(r,s,t)

	// Private methods
	void _shape_2d (double r, double s);
	void _shape_3d (double r, double s, double t);

}; // class Structured


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// ------------------------------------------------------------------------------------------------------- Block

/* public */

inline void Block::Set(int Tag, Verts_T const & Verts, Edges_T const & Edges, ETags_T const * ETags, FTags_T const * FTags, long OrigID, long XPlusID, long YPlusID, long ZPlusID)
{
	// Check
	bool   gen_mid = false; // generate mid vertices ?
	size_t nedges  = Edges.Size();
	size_t nverts  = Verts.Size();
	if (nedges==24)
	{
		if (nverts!=20) throw new Fatal("Block::Set:: For 3D blocks with 24 edges, the number of vertices must be 20. (nverts==%d is invalid)",nverts);
		_is_3d = true;
	}
	else if (nedges==12)
	{
		if (nverts!=8) throw new Fatal("Block::Set:: For 3D blocks with 12 edges, the number of vertices must be 8. (nverts==%d is invalid)",nverts);
		_is_3d  = true;
		gen_mid = true;
	}
	else if (nedges==8)
	{
		if (nverts!=8) throw new Fatal("Block::Set:: For 2D blocks with 8 edges, the number of vertices must be 8. (nverts==%d is invalid)",nverts);
		_is_3d = false;
	}
	else if (nedges==4)
	{
		if (nverts!=4) throw new Fatal("Block::Set:: For 2D blocks with 4 edges, the number of vertices must be 4 (nverts==%d is invalid)",nverts);
		_is_3d  = false;
		gen_mid = true;
	}
	else throw new Fatal("Block::Set:: The list with edges must have 4(2D), 8(2D), 12(3D), or 24(3D) items. (nedges==%d is invalid)",nedges);

	// Set
	_tag = Tag;
	if (_is_3d) _set_3d(); else _set_2d();

	// Edges
	Array<long> edges(nedges*2); // *2 => we have to serialize for Util::Tree
	int k = 0;
	for (size_t i=0; i<nedges; ++i)
	{
		edges[k  ] = Edges[i].get<0>();
		edges[k+1] = Edges[i].get<1>();
		k += 2;
	}

	// Tree
	Util::Tree tree(edges);
	Array<long> eds;   eds.Resize(nedges);  eds.SetNS(Util::_4);
	Array<long> vl2g; vl2g.Resize(nverts); vl2g.SetNS(Util::_4); // map: local vertex to global vertex IDs

	// Sort
	vl2g[0] = OrigID;
	Array<long> path;
	tree.DelEdge   (OrigID, XPlusID);
	tree.ShortPath (XPlusID, YPlusID, path);
	if (_is_3d)
	{
		if (nedges==24)
		{
			// Bottom nodes
			vl2g[ 8] = path[0];
			vl2g[ 1] = path[1];
			vl2g[ 9] = path[2];
			vl2g[ 2] = path[3];
			vl2g[10] = path[4];
			vl2g[ 3] = path[5];
			vl2g[11] = path[6];

			// Behind nodes
			tree.DelEdge   (OrigID, YPlusID);
			tree.ShortPath (YPlusID, ZPlusID, path);
			vl2g[19] = path[2];
			vl2g[ 7] = path[3];
			vl2g[15] = path[4];
			vl2g[ 4] = path[5];
			vl2g[16] = path[6];

			// Left nodes
			tree.DelEdge   (OrigID, ZPlusID);
			tree.ShortPath (ZPlusID, XPlusID, path);
			vl2g[12] = path[2];
			vl2g[ 5] = path[3];
			vl2g[17] = path[4];

			// Corner-front nodes
			tree.DelEdge   (vl2g[7], vl2g[15]);
			tree.DelEdge   (vl2g[7], vl2g[19]);
			tree.ShortPath (vl2g[7], vl2g[ 5], path);
			vl2g[14] = path[1];
			vl2g[ 6] = path[2];
			vl2g[13] = path[3];
			tree.ShortPath (vl2g[7], vl2g[2], path);
			vl2g[18] = path[3];
		}
		else // nedges==12
		{
			// Bottom nodes
			vl2g[1] = path[0];
			vl2g[2] = path[1];
			vl2g[3] = path[2];

			// Behind nodes
			tree.DelEdge   (OrigID, YPlusID);
			tree.ShortPath (YPlusID, ZPlusID, path);
			vl2g[7] = path[1];
			vl2g[4] = path[2];

			// Left nodes
			tree.DelEdge   (OrigID, ZPlusID);
			tree.ShortPath (ZPlusID, XPlusID, path);
			vl2g[5] = path[1];

			// Corner-front nodes
			tree.DelEdge   (vl2g[4], vl2g[7]);
			tree.DelEdge   (vl2g[7], vl2g[3]);
			tree.ShortPath (vl2g[7], vl2g[5], path);
			vl2g[6] = path[1];
		}
	}
	else
	{
		if (nedges==8)
		{
			vl2g[4] = path[0];
			vl2g[1] = path[1];
			vl2g[5] = path[2];
			vl2g[2] = path[3];
			vl2g[6] = path[4];
			vl2g[3] = path[5];
			vl2g[7] = path[6];
		}
		else // nedges==4
		{
			vl2g[1] = path[0];
			vl2g[2] = path[1];
			vl2g[3] = path[2];
		}
	}
	//std::cout << "vl2g = " << vl2g << std::endl;

	// Edge and face tags
	size_t net = (ETags==NULL ? 0 : ETags->Size()); // number of edges with tags
	size_t nft = (FTags==NULL ? 0 : FTags->Size()); // number of faces with tags
	Array<long> eg2l; eg2l.SetNS(Util::_4); // map: global edge to local edge IDs
	if (net>0 || nft>0)
	{
		eg2l.Resize(nedges);
		tree.Reset (edges);
		if (_is_3d)
		{
			if (nedges==24)
			{
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
			else // nedges==12
			{
				eg2l[tree.GetEdge(vl2g[0],vl2g[3])] =  0;
				eg2l[tree.GetEdge(vl2g[1],vl2g[2])] =  1;
				eg2l[tree.GetEdge(vl2g[0],vl2g[1])] =  2;
				eg2l[tree.GetEdge(vl2g[3],vl2g[2])] =  3;
				eg2l[tree.GetEdge(vl2g[4],vl2g[7])] =  4;
				eg2l[tree.GetEdge(vl2g[5],vl2g[6])] =  5;
				eg2l[tree.GetEdge(vl2g[4],vl2g[5])] =  6;
				eg2l[tree.GetEdge(vl2g[7],vl2g[6])] =  7;
				eg2l[tree.GetEdge(vl2g[0],vl2g[4])] =  8;
				eg2l[tree.GetEdge(vl2g[1],vl2g[5])] =  9;
				eg2l[tree.GetEdge(vl2g[2],vl2g[6])] = 10;
				eg2l[tree.GetEdge(vl2g[3],vl2g[7])] = 11;
			}
		}
		else
		{
			if (nedges==8)
			{
				eg2l[tree.GetEdge(vl2g[0],vl2g[7])] = 0;
				eg2l[tree.GetEdge(vl2g[1],vl2g[5])] = 1;
				eg2l[tree.GetEdge(vl2g[0],vl2g[4])] = 2;
				eg2l[tree.GetEdge(vl2g[3],vl2g[6])] = 3;
				eg2l[tree.GetEdge(vl2g[3],vl2g[7])] = 4;
				eg2l[tree.GetEdge(vl2g[2],vl2g[5])] = 5;
				eg2l[tree.GetEdge(vl2g[1],vl2g[4])] = 6;
				eg2l[tree.GetEdge(vl2g[2],vl2g[6])] = 7;
			}
			else // nedges==4
			{
				eg2l[tree.GetEdge(vl2g[0],vl2g[3])] = 0;
				eg2l[tree.GetEdge(vl2g[1],vl2g[2])] = 1;
				eg2l[tree.GetEdge(vl2g[0],vl2g[1])] = 2;
				eg2l[tree.GetEdge(vl2g[3],vl2g[2])] = 3;
			}
		}
	}
	//std::cout << "eg2l = " << eg2l << std::endl;

	// Read coordinates
	for (size_t i=0; i<nverts; ++i)
	{
		long id    = vl2g[i];
		bool found = false;
		for (size_t j=0; j<nverts; ++j)
		{
			if (Verts[j].get<0>()==id)
			{
				            _c(0,i) = Verts[j].get<1>();
				            _c(1,i) = Verts[j].get<2>();
				if (_is_3d) _c(2,i) = Verts[j].get<3>();
				found = true;
				break;
			}
		}
		if (found==false) throw new Fatal("Block::Set: Could not find vertex ID == %d in Array of Verts",id);
	}

	// Generate mid nodes
	if (gen_mid) _gen_mid_nodes ();

	// Edge tags
	_has_etags = (net>0 ? true : false);
	for (size_t i=0; i<net; ++i)
	{
		size_t eid = tree.GetEdge ((*ETags)[i].get<0>(), (*ETags)[i].get<1>());
		_etags(eg2l[eid]) = (*ETags)[i].get<2>();
	}

	// Face tags
	if (_is_3d)
	{
		_has_ftags = (nft>0 ? true : false);
		for (size_t i=0; i<nft; ++i)
		{
			// face tag
			int ftag = (*FTags)[i].get<4>();

			// 3 edges defining the face
			long eid0 = eg2l[ tree.GetEdge ((*FTags)[i].get<0>(), (*FTags)[i].get<1>()) ];
			long eid1 = eg2l[ tree.GetEdge ((*FTags)[i].get<1>(), (*FTags)[i].get<2>()) ];
			long eid2 = eg2l[ tree.GetEdge ((*FTags)[i].get<2>(), (*FTags)[i].get<3>()) ];

			// find face given 3 edges
			int  sum[6] = { 0,0,0,0,0,0 }; // number of edges defining each face
			sum[Edge2Face[eid0].L] += 1;   sum[Edge2Face[eid0].R] += 1;
			sum[Edge2Face[eid1].L] += 1;   sum[Edge2Face[eid1].R] += 1;
			sum[Edge2Face[eid2].L] += 1;   sum[Edge2Face[eid2].R] += 1;
			bool found = false; // found face with edges eid0,eid1,eid2 ?
			for (int j=0; j<6; ++j) // loop over face local IDs
			{
				if (sum[j]==3) // found face
				{
					_ftags(j) = ftag;
					found     = true;
					break;
				}
			}
			if (found==false) throw new Fatal("Block::Set: Could not find face with vertices (%d,%d,%d,%d) and ftag==%d in block with blktag==%d", (*FTags)[i].get<0>(), (*FTags)[i].get<1>(), (*FTags)[i].get<2>(), (*FTags)[i].get<3>(), ftag, Tag);
		}
	}
	else
	{
		if (nft>0) throw new Fatal("Block::SetCoords: For 2D blocks, FTags must be null == []");
	}
}

inline void Block::SetNx(size_t Nx, double Ax, bool NonLin)
{
	_wx.Resize(Nx);
	if (NonLin) for (size_t i=0; i<Nx; ++i) _wx[i] = pow(i+1.0,Ax);
	else        for (size_t i=0; i<Nx; ++i) _wx[i] = 1.0+Ax*i;
	_set_wx();
}

inline void Block::SetNy(size_t Ny, double Ay, bool NonLin)
{
	_wy.Resize(Ny);
	if (NonLin) for (size_t i=0; i<Ny; ++i) _wy[i] = pow(i+1.0,Ay);
	else        for (size_t i=0; i<Ny; ++i) _wy[i] = 1.0+Ay*i;
	_set_wy();
}

inline void Block::SetNz(size_t Nz, double Az, bool NonLin)
{
	_wz.Resize(Nz);
	if (NonLin) for (size_t i=0; i<Nz; ++i) _wz[i] = pow(i+1.0,Az);
	else        for (size_t i=0; i<Nz; ++i) _wz[i] = 1.0+Az*i;
	_set_wz();
}

inline double Block::Diagonal() const
{
	if (_is_3d) return sqrt(pow(_c(0,6)-_c(0,0),2.0)+pow(_c(1,6)-_c(1,0),2.0)+pow(_c(2,6)-_c(2,0),2.0));
	else        return sqrt(pow(_c(0,2)-_c(0,0),2.0)+pow(_c(1,2)-_c(1,0),2.0));
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

inline void Block::PySet(int Tag, BPy::list const & Verts, BPy::list const & Edges, BPy::list const & ETags, BPy::list const & FTags, long OrigID, long XPlusID, long YPlusID, long ZPlusID)
{
	int nverts = BPy::len(Verts);
	int nedges = BPy::len(Edges);
	int netags = BPy::len(ETags);
	int nftags = BPy::len(FTags);
	Mesh::Verts_T verts(nverts);
	Mesh::Edges_T edges(nedges);
	Mesh::ETags_T etags; if (netags>0) etags.Resize(netags);
	Mesh::FTags_T ftags; if (nftags>0) ftags.Resize(nftags);

	// Extract verts
	for (int i=0; i<nverts; ++i)
	{
		BPy::tuple const & ve = BPy::extract<BPy::tuple>(Verts[i])();
		     if (BPy::len(ve)==3) verts[i] = boost::make_tuple (BPy::extract<int>(ve[0])(), BPy::extract<double>(ve[1])(), BPy::extract<double>(ve[2])(), 0.0);
		else if (BPy::len(ve)==4) verts[i] = boost::make_tuple (BPy::extract<int>(ve[0])(), BPy::extract<double>(ve[1])(), BPy::extract<double>(ve[2])(), BPy::extract<double>(ve[3])());
		else throw new Fatal("Block::PySet: The list of vertices (Verts) must contain tuples with 3 or 4 elements. Ex.: [(ID,x,y)...] or [(ID,x,y,z)...]");
	}

	// Extract edges
	for (int i=0; i<nedges; ++i)
	{
		BPy::tuple const & ed = BPy::extract<BPy::tuple>(Edges[i])();
		if (BPy::len(ed)!=2) throw new Fatal("Block::PySet: The list of edges (Edges) must contain tuples with 2 elements. Ex.: [(v1,v2)...]");
		edges[i] = boost::make_tuple (BPy::extract<int>(ed[0])(), BPy::extract<int>(ed[1])());
	}

	// Extract etags
	for (int i=0; i<netags; ++i)
	{
		BPy::tuple const & et = BPy::extract<BPy::tuple>(ETags[i])();
		if (BPy::len(et)!=3) throw new Fatal("Block::PySet: The list of edges tags (ETags) must contain tuples with 3 elements. Ex.: [(v1,v2,tag)...]");
		etags[i] = boost::make_tuple (BPy::extract<int>(et[0])(), BPy::extract<int>(et[1])(), BPy::extract<int>(et[2])());
	}

	// Extract ftags
	for (int i=0; i<nftags; ++i)
	{
		BPy::tuple const & ft = BPy::extract<BPy::tuple>(FTags[i])();
		if (BPy::len(ft)!=5) throw new Fatal("Block::PySet: The list of faces tags (FTags) must contain tuples with 5 elements. Ex.: [(v1,v2,v3,v4,tag)...]");
		ftags[i] = boost::make_tuple (BPy::extract<int>(ft[0])(), BPy::extract<int>(ft[1])(), BPy::extract<int>(ft[2])(), BPy::extract<int>(ft[3])(), BPy::extract<int>(ft[4])());
	}

	// Set Block
	Set (Tag, verts, edges, (netags>0?&etags:NULL), (nftags>0?&ftags:NULL), OrigID, XPlusID, YPlusID, ZPlusID);
}

#endif // USE_BOOST_PYTHON


/* private */

inline void Block::_set_3d()
{
	_is_3d = true;
	_c    .Resize    (3,20);
	_etags.Resize    (24);
	_ftags.Resize    (6);
	_etags.SetValues (0);
	_ftags.SetValues (0);
}

inline void Block::_set_2d()
{
	_is_3d = false;
	_c.Resize        (2,8);
	_etags.Resize    (8);
	_etags.SetValues (0);
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

inline void Block::_gen_mid_nodes()
{
	if (_is_3d)
	{
		for (size_t i=0; i<3; ++i)
		{
			_c(i, 8) = (_c(i,0) + _c(i,1))/2.0;
			_c(i, 9) = (_c(i,1) + _c(i,2))/2.0;
			_c(i,10) = (_c(i,2) + _c(i,3))/2.0;
			_c(i,11) = (_c(i,3) + _c(i,0))/2.0;

			_c(i,12) = (_c(i,4) + _c(i,5))/2.0;
			_c(i,13) = (_c(i,5) + _c(i,6))/2.0;
			_c(i,14) = (_c(i,6) + _c(i,7))/2.0;
			_c(i,15) = (_c(i,7) + _c(i,4))/2.0;

			_c(i,16) = (_c(i,0) + _c(i,4))/2.0;
			_c(i,17) = (_c(i,1) + _c(i,5))/2.0;
			_c(i,18) = (_c(i,2) + _c(i,6))/2.0;
			_c(i,19) = (_c(i,3) + _c(i,7))/2.0;
		}
	}
	else
	{
		for (size_t i=0; i<2; ++i)
		{
			_c(i,4) = (_c(i,0) + _c(i,1))/2.0;
			_c(i,5) = (_c(i,1) + _c(i,2))/2.0;
			_c(i,6) = (_c(i,2) + _c(i,3))/2.0;
			_c(i,7) = (_c(i,3) + _c(i,0))/2.0;
		}
	}
}

// --------------------------------------------------------------------------------------------------- Structured

/* public */

inline size_t Structured::Generate(bool WithInfo)
{
	// Info
	double start = std::clock();

	// Check
	if (_bls.Size()<1) throw new Fatal("Structured::Generate: Number of blocks must be greater than 0 (%d is invalid)",_bls.Size());

	// Erase previous mesh
	Erase();

	// Check if the first block is 3D
	_is_3d = _bls[0].Is3D();
	_s.Resize((_is_3d ? 20 : 8)); // resize the shape values vector

	// Generate vertices and elements (with duplicates)
	double min_diagonal = _bls[0].Diagonal(); // minimum diagonal among all elements
	for (size_t b=0; b<_bls.Size(); ++b)
	{
		// Check if all blocks have the same space dimension
		if (_bls[b].Is3D()!=_is_3d) throw new Fatal("Structured::Generate: All blocks must have the same space dimension");
		
		// Check if block has all arrays/vectors/matrices with the right sizes
		_bls[b].Alright();

		// Generate
		double t_before = -2.0;
		double t = -1.0; // initial Z natural coordinate
		for (int k=0; k<(_is_3d ? _bls[b].nDivZ()+1 : 1); ++k)
		{
			double s_before = -2.0;
			double s = -1.0; // initial Y natural coordinate
			for (int j=0; j<_bls[b].nDivY()+1; ++j)
			{
				double r_before = -2.0;
				double r = -1.0; // initial X natural coordinate
				for (int i=0; i<_bls[b].nDivX()+1; ++i)
				{
					// Compute shape (interpolation) functions
					if (_is_3d) _shape_3d (r,s,t);
					else        _shape_2d (r,s);

					// New vertex
					Vertex  * v = new Vertex;
					v->MyID     = _verts_d.Size();           // id
					v->BlockNum = b;                         // block number (used when removing duplicates)
					v->C        = _bls[b].C() * _s;          // new x-y-z coordinates
					v->Dupl     = false;                     // is this a duplicated node?
					_bls[b].FindLocalEdgesFacesID(i,j,k, v); // check if it is on boundary and find EdgesID where this vertex is located on
					_verts_d.Push(v);
					if (v->OnBry) _verts_d_bry.Push(v); // array with vertices on boundary

					// New O2 vertices
					Vertex * v1 = NULL;
					Vertex * v2 = NULL;
					Vertex * v3 = NULL;
					if (_is_o2)
					{
						if (i!=0)
						{
							if (_is_3d) _shape_3d ((r_before+r)/2.0,s,t);
							else        _shape_2d ((r_before+r)/2.0,s);
							v1           = new Vertex;
							v1->MyID     = _verts_m1.Size();
							v1->BlockNum = b;
							v1->C        = _bls[b].C() * _s;
							v1->Dupl     = false;
							if (i==_bls[b].nDivX()) _bls[b].FindLocalEdgesFacesID(i-1,j,k, v1);
							else                    _bls[b].FindLocalEdgesFacesID(i,  j,k, v1);
							_verts_m1.Push(v1);
							if (v1->OnBry) _verts_d_bry.Push(v1);
						}
						if (j!=0)
						{
							if (_is_3d) _shape_3d (r,(s_before+s)/2.0,t);
							else        _shape_2d (r,(s_before+s)/2.0);
							v2           = new Vertex;
							v2->MyID     = _verts_m2.Size();
							v2->BlockNum = b;
							v2->C        = _bls[b].C() * _s;
							v2->Dupl     = false;
							if (j==_bls[b].nDivY()) _bls[b].FindLocalEdgesFacesID(i,j-1,k, v2);
							else                    _bls[b].FindLocalEdgesFacesID(i,j,  k, v2);
							_verts_m2.Push(v2);
							if (v2->OnBry) _verts_d_bry.Push(v2);
						}
						if (k!=0)
						{
							if (_is_3d) _shape_3d (r,s,(t_before+t)/2.0);
							v3           = new Vertex;
							v3->MyID     = _verts_m3.Size();
							v3->BlockNum = b;
							v3->C        = _bls[b].C() * _s;
							v3->Dupl     = false;
							if (k==_bls[b].nDivZ()) _bls[b].FindLocalEdgesFacesID(i,j,k-1, v3);
							else                    _bls[b].FindLocalEdgesFacesID(i,j,k,   v3);
							_verts_m3.Push(v3);
							if (v3->OnBry) _verts_d_bry.Push(v3);
						}
					}

					// New element
					double diagonal = 0.0; // element diagonal
					if (i!=0 && j!=0 && (_is_3d ? k!=0 : true))
					{
						Elem * e = new Elem;
						e->MyID     = _elems.Size(); // id
						e->Tag      = _bls[b].Tag(); // tag
						if (_is_3d)
						{
							// connectivity
							e->V.Resize((_is_o2?20:8));
							e->V[0] = _verts_d[v->MyID - 1 - (_bls[b].nDivX()+1) - (_bls[b].nDivX()+1)*(_bls[b].nDivY()+1)];
							e->V[1] = _verts_d[v->MyID     - (_bls[b].nDivX()+1) - (_bls[b].nDivX()+1)*(_bls[b].nDivY()+1)];
							e->V[2] = _verts_d[v->MyID                           - (_bls[b].nDivX()+1)*(_bls[b].nDivY()+1)];
							e->V[3] = _verts_d[v->MyID - 1                       - (_bls[b].nDivX()+1)*(_bls[b].nDivY()+1)];
							e->V[4] = _verts_d[v->MyID - 1 - (_bls[b].nDivX()+1)];
							e->V[5] = _verts_d[v->MyID -     (_bls[b].nDivX()+1)];
							e->V[6] = _verts_d[v->MyID];
							e->V[7] = _verts_d[v->MyID - 1];
							if (_is_o2)
							{
								e->V[ 8] = _verts_m1[v1->MyID - _bls[b].nDivX() - _bls[b].nDivX()*(_bls[b].nDivY()+1)];
								e->V[ 9] = _verts_m2[v2->MyID                   - _bls[b].nDivY()*(_bls[b].nDivX()+1)];
								e->V[10] = _verts_m1[v1->MyID                   - _bls[b].nDivX()*(_bls[b].nDivY()+1)];
								e->V[11] = _verts_m2[v2->MyID -1                - _bls[b].nDivY()*(_bls[b].nDivX()+1)];

								e->V[12] = _verts_m1[v1->MyID - _bls[b].nDivX()];
								e->V[13] = _verts_m2[v2->MyID];
								e->V[14] = _verts_m1[v1->MyID];
								e->V[15] = _verts_m2[v2->MyID - 1];

								e->V[16] = _verts_m3[v3->MyID - 1 - (_bls[b].nDivX()+1)];
								e->V[17] = _verts_m3[v3->MyID -     (_bls[b].nDivX()+1)];
								e->V[18] = _verts_m3[v3->MyID];
								e->V[19] = _verts_m3[v3->MyID - 1];
							}
							// shares information
							for (size_t m=0; m<(_is_o2?20:8); ++m)
							{
								Share s = {e,m}; // The node shares information will have this element and the local node index
								e->V[m]->Shares.Push(s);
							}
							// is on boundary ?
							e->OnBry = (e->V[0]->OnBry || e->V[1]->OnBry || e->V[2]->OnBry || e->V[3]->OnBry || e->V[4]->OnBry || e->V[5]->OnBry || e->V[6]->OnBry || e->V[7]->OnBry); // if any node is on Bry, yes
							// VTK cell type
							e->VTKCellType = (_is_o2 ? VTK_QUADRATIC_HEXAHEDRON : VTK_HEXAHEDRON);
							// Diagonal
							diagonal = sqrt(pow(e->V[6]->C(0) - e->V[0]->C(0),2.0)+
							                pow(e->V[6]->C(1) - e->V[0]->C(1),2.0)+
							                pow(e->V[6]->C(2) - e->V[0]->C(2),2.0));
						}
						else
						{
							// connectivity
							e->V.Resize((_is_o2?8:4));
							e->V[0] = _verts_d[v->MyID - 1 - (_bls[b].nDivX()+1)];
							e->V[1] = _verts_d[v->MyID     - (_bls[b].nDivX()+1)];
							e->V[2] = _verts_d[v->MyID];
							e->V[3] = _verts_d[v->MyID - 1];
							if (_is_o2)
							{
								e->V[4] = _verts_m1[v1->MyID - _bls[b].nDivX()];
								e->V[5] = _verts_m2[v2->MyID];
								e->V[6] = _verts_m1[v1->MyID];
								e->V[7] = _verts_m2[v2->MyID - 1];
							}
							// shares information
							for (size_t m=0; m<(_is_o2?8:4); ++m)
							{
								Share s = {e,m};
								e->V[m]->Shares.Push(s);
							}
							// is on boundary ?
							e->OnBry = (e->V[0]->OnBry || e->V[1]->OnBry || e->V[2]->OnBry || e->V[3]->OnBry); // if any node is on Bry, yes
							// VTK cell type
							e->VTKCellType = (_is_o2 ? VTK_QUADRATIC_QUAD : VTK_QUAD);
							// Diagonal
							diagonal = sqrt(pow(e->V[2]->C(0) - e->V[0]->C(0),2.0)+
							                pow(e->V[2]->C(1) - e->V[0]->C(1),2.0));
						}
						if (e->OnBry)
						{
							_elems_bry.Push(e);   // array with elements on boundary
							_bls[b].ApplyTags(e); // apply tags to edges and faces of this element with boundary tags
						}
						_elems.Push(e); // array with all elements
						// Minimum diagonal
						if (diagonal<min_diagonal) min_diagonal = diagonal;
					}
					// Next r
					r_before = r;
					r += (2.0/_bls[b].SumWeightX()) * _bls[b].Wx(i);
				}
				// Next s
				s_before = s;
				s += (2.0/_bls[b].SumWeightY()) * _bls[b].Wy(j);
			}
			// Next t
			t_before = t;
			t += (_is_3d ? (2.0/_bls[b].SumWeightZ()) * _bls[b].Wz(k) : 0.0);
		}
	}

	// Remove duplicates
	long ncomp = 0; // number of comparisons
	long ndupl = 0; // number of duplicates
	double tol = min_diagonal*0.001; // tolerance decide whether two vertices are coincident or not
	if (_bls.Size()>1)
	{
		for (size_t i=0; i<_verts_d_bry.Size(); ++i)
		{
			for (size_t j=i+1; j<_verts_d_bry.Size(); ++j)
			{
				if (_verts_d_bry[i]->BlockNum!=_verts_d_bry[j]->BlockNum) // Vertices are located on different blocks
				{
					// check distance
					double dist = sqrt(          pow(_verts_d_bry[i]->C(0)-_verts_d_bry[j]->C(0),2.0)+
					                             pow(_verts_d_bry[i]->C(1)-_verts_d_bry[j]->C(1),2.0)+
					                   (_is_3d ? pow(_verts_d_bry[i]->C(2)-_verts_d_bry[j]->C(2),2.0) : 0.0));
					if (dist<tol)
					{
						/* TODO: this is wrong, since corner nodes can be duplicated and are still on boundary
						// If this node is duplicated, than it is not on-boundary any longer
						_verts_d_bry[i]->OnBry = false;
						*/
						// Mark duplicated
						if (_verts_d_bry[j]->Dupl==false) // vertex not tagged as duplicated yet
						{
							ndupl++;
							_verts_d_bry[j]->Dupl = true;
							// Chage elements' connectivities
							for (size_t k=0; k<_verts_d_bry[j]->Shares.Size(); ++k)
							{
								Elem * e = _verts_d_bry[j]->Shares[k].E;
								int    n = _verts_d_bry[j]->Shares[k].N;
								e->V[n] = _verts_d_bry[i];
							}
						}
					}
					ncomp++;
				}
			}
		}
	}

	// Set new array with non-duplicated vertices
	size_t k = 0;
	size_t m = 0;
	_verts    .Resize (_verts_d    .Size() + _verts_m1.Size() + _verts_m2.Size() + _verts_m3.Size() - ndupl);
	_verts_bry.Resize (_verts_d_bry.Size() - ndupl);
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
	for (size_t i=0; i<_verts_m1.Size(); ++i)
	{
		if (_verts_m1[i]->Dupl==false)
		{
			_verts[k]       = _verts_m1[i]; // copy pointer
			_verts[k]->MyID = k;            // new ID
			if (_verts[k]->OnBry)
			{
				_verts_bry[m] = _verts[k]; // copy pointer
				m++;
			}
			k++;
		}
	}
	for (size_t i=0; i<_verts_m2.Size(); ++i)
	{
		if (_verts_m2[i]->Dupl==false)
		{
			_verts[k]       = _verts_m2[i]; // copy pointer
			_verts[k]->MyID = k;            // new ID
			if (_verts[k]->OnBry)
			{
				_verts_bry[m] = _verts[k]; // copy pointer
				m++;
			}
			k++;
		}
	}
	for (size_t i=0; i<_verts_m3.Size(); ++i)
	{
		if (_verts_m3[i]->Dupl==false)
		{
			_verts[k]       = _verts_m3[i]; // copy pointer
			_verts[k]->MyID = k;            // new ID
			if (_verts[k]->OnBry)
			{
				_verts_bry[m] = _verts[k]; // copy pointer
				m++;
			}
			k++;
		}
	}

	// Info
	if (WithInfo)
	{
		double total = std::clock() - start;
		std::cout << "[1;33m\n--- Structured Mesh Generation ---------------------------------[0m\n";
		if (_is_o2) std::cout << "[1;36m    Time elapsed (o2)     = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";
		else        std::cout << "[1;36m    Time elapsed          = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";
		std::cout << "[1;35m    Number of comparisons = " << ncomp         << "[0m\n";
		std::cout << "[1;35m    Number of duplicates  = " << ndupl         << "[0m\n";
		std::cout << "[1;35m    Minimum diagonal      = " << min_diagonal  << "[0m\n";
		std::cout << "[1;32m    Number of elements    = " << _elems.Size() << "[0m" << std::endl;
	}

	return _elems.Size();
}

inline void Structured::GenBox(int Nx, int Ny, int Nz, double Lx, double Ly, double Lz)
{
    /*
                      4----------------7  
                    ,'|              ,'| 
                  ,'  |            ,'  | 
                ,'    | -6    -1 ,'    | 
              ,'      |        ,'      | 
            5'===============6'        | 
            |         |      |    -4   | 
            |    -3   |      |         | 
            |         0- - - | -  - - -3  
            |       ,'       |       ,'  
            |     ,' -2      |     ,'    
            |   ,'        -5 |   ,'      
            | ,'             | ,'        
            1----------------2'          
    */

	// Blocks
	Array<Mesh::Block> bks(1);

	// Block # 0
	Mesh::Verts_T ve0( 8);
	Mesh::Edges_T ed0(12);
	Mesh::FTags_T ft0( 6);
	ve0 = boost::make_tuple(0, 0.,0.,0.), boost::make_tuple(1, Lx,0.,0.), boost::make_tuple(2, Lx,Ly,0.), boost::make_tuple(3, 0.,Ly,0.),
	      boost::make_tuple(4, 0.,0.,Lz), boost::make_tuple(5, Lx,0.,Lz), boost::make_tuple(6, Lx,Ly,Lz), boost::make_tuple(7, 0.,Ly,Lz);
	ed0 = boost::make_tuple(0,1), boost::make_tuple(1,2), boost::make_tuple(2,3), boost::make_tuple(3,0),
	      boost::make_tuple(4,5), boost::make_tuple(5,6), boost::make_tuple(6,7), boost::make_tuple(7,4),
	      boost::make_tuple(0,4), boost::make_tuple(1,5), boost::make_tuple(2,6), boost::make_tuple(3,7);
	ft0 = boost::make_tuple(0,3,7,4,-1), boost::make_tuple(1,2,6,5,-2), boost::make_tuple(1,0,4,5,-3),
	      boost::make_tuple(2,3,7,6,-4), boost::make_tuple(0,1,2,3,-5), boost::make_tuple(4,5,6,7,-6);
	bks[0].Set   (-1, ve0, ed0, NULL, &ft0, /*orig*/0, /*xplus*/1, /*yplus*/3, /*zplus*/4);
	bks[0].SetNx (Nx);
	bks[0].SetNy (Ny);
	bks[0].SetNz (Nz);

	// Generate
	SetBlocks (bks);
	Generate  (/*WithInfo*/false);
}

#ifdef USE_BOOST_PYTHON

inline void Structured::PySetBlocks(BPy::list const & ListOfMeshBlock)
{
	int nb = BPy::len(ListOfMeshBlock);
	if (nb<1) throw new Fatal("Structured::PyGenerate: Number of blocks must be greater than 0 (%d is invalid)",nb);
	_bls.Resize (nb);
	for (int i=0; i<nb; ++i) _bls[i] = BPy::extract<Mesh::Block const &>(ListOfMeshBlock[i])();
}

#endif // USE_BOOST_PYTHON


/* private */

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

inline void Structured::Erase()
{
	for (size_t i=0; i<_verts_d .Size(); ++i) if (_verts_d [i]!=NULL) delete _verts_d [i]; // it is only necessary to delete nodes in _verts_d array
	for (size_t i=0; i<_verts_m1.Size(); ++i) if (_verts_m1[i]!=NULL) delete _verts_m1[i]; // it is only necessary to delete nodes in _verts_m1 array
	for (size_t i=0; i<_verts_m2.Size(); ++i) if (_verts_m2[i]!=NULL) delete _verts_m2[i]; // it is only necessary to delete nodes in _verts_m2 array
	for (size_t i=0; i<_verts_m3.Size(); ++i) if (_verts_m3[i]!=NULL) delete _verts_m3[i]; // it is only necessary to delete nodes in _verts_m3 array
	for (size_t i=0; i<_elems   .Size(); ++i) if (_elems   [i]!=NULL) delete _elems   [i]; // it is only necessary to delete elems in _elems array
	_verts_d    .Resize(0);
	_verts_m1   .Resize(0);
	_verts_m2   .Resize(0);
	_verts_m3   .Resize(0);
	_verts_d_bry.Resize(0);
	_verts      .Resize(0);
	_elems      .Resize(0);
	_elems_bry  .Resize(0);
	_verts_bry  .Resize(0);
}

}; // namespace Mesh

#endif // MECHSYS_MESH_STRUCTURED_H
