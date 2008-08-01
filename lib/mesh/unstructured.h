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

#ifndef MECHSYS_MESH_UNSTRUCTURED_H
#define MECHSYS_MESH_UNSTRUCTURED_H

/* LOCAL indexes of Vertices, Edges, and Faces

  2D:
             Nodes                 Faces 
                           
   y           2                                                        
   |           @                     @                                  
   +--x       / \                   / \                                  
           5 /   \ 4               /   \                                 
            @     @             2 /     \ 1                              
           /       \             /       \                               
          /         \           /         \                              
         @-----@-----@         @-----------@                             
        0      3      1              0               
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

// Jonathan R Shewchuk' Triangle
extern "C"
{
	#define REAL double
	#define ANSI_DECLARATORS
	#define VOID int
	  #include "jrs_triangle.h"
	#undef REAL
	#undef ANSI_DECLARATORS
	#undef VOID
}

// MechSys
#include "util/array.h"
#include "util/exception.h"
#include "util/lineparser.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/laexpr.h"
#include "mesh/mesh.h"

using LinAlg::Vector;
using LinAlg::Matrix;
using blitz::TinyVector;

namespace Mesh
{

class Region
{
}; // class Region

class Hole
{
}; // class Hole

/** JRS' Triangle Input/Output structure. */
typedef triangulateio TriIO;

class Unstructured : public virtual Mesh::Generic
{
public:
	// Constants
	static Edge Edge2Vert[]; ///< Map from local edge ID to local vertex ID

	// Constructor
	Unstructured (double Tol=sqrt(DBL_EPSILON)); ///< Tol is the tolerance to regard two vertices as coincident

	// Destructor
	~Unstructured () { _tri_deallocate_all(_tin); _tri_deallocate_all(_tou); }

	// Set Methods
	void SetPolySize    (size_t NPoints, size_t NSegments, size_t NRegions=0, size_t NHoles=0);                                                                                                 ///< Erase any previous input PSLG and set the number of points and segments of the polygon. Also set the number of holes and regions.
	void SetPolyPoint   (size_t i, double X, double Y)                                { _tin.pointlist  [i*2]=X;          _tin.pointlist  [i*2+1]=Y;           }                                ///< Set the coordinates of point i of input PSLG. SetPolySize MUST be called first.
	void SetPolySegment (size_t i, size_t iPointLeft, size_t iPointRight, long Tag=0) { _tin.segmentlist[i*2]=iPointLeft; _tin.segmentlist[i*2+1]=iPointRight; _tin.segmentmarkerlist[i]=Tag; } ///< Set the left and right points of a segment i of input PSLG. SetPolySize MUST be called first.
	void SetPolyRegion  (size_t i, double X, double Y)                                { _tin.regionlist [i*2]=X;          _tin.regionlist [i*2+1]=Y;           }                                ///< Set the coordinates of a point defining a region i.
	void SetPolyHole    (size_t i, double X, double Y)                                { _tin.holelist   [i*2]=X;          _tin.holelist   [i*2+1]=Y;           }                                ///< Set the coordinates of a point defining a hole i.

	// Methods
	size_t Generate ();

#ifdef USE_BOOST_PYTHON
// {
	//size_t PyGenerate (BPy::list const & Regions, BPy::list const & Holes);
// }
#endif

private:
	// Data
	TriIO _tin; ///< Triangle IO's input structure
	TriIO _tou; ///< Triangle IO's output structure

	// Overloaded private methods
	void _vtk_con          (Elem const * E, String & Connect) const;
	int  _edge_to_lef_vert (int EdgeLocalID) const { return Edge2Vert[EdgeLocalID].L; }
	int  _edge_to_rig_vert (int EdgeLocalID) const { return Edge2Vert[EdgeLocalID].R; }

	// Private methods
	void _tri_set_all_to_null (TriIO & Tio); ///< Set all elements of Triangle's IO structure to NULL and 0
	void _tri_deallocate_all  (TriIO & Tio); ///< Deallocate all arrays inside Triangle's IO structure

}; // class Unstructured

Edge Unstructured::Edge2Vert[]= {{ 0, 1 },
                                 { 1, 2 },
                                 { 0, 2 }};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline Unstructured::Unstructured(double Tol)
	: Mesh::Generic (Tol)
{
	_tri_set_all_to_null (_tin);
	_tri_set_all_to_null (_tou);
}

inline void Unstructured::SetPolySize(size_t NPoints, size_t NSegments, size_t NRegions, size_t NHoles)
{
	// Erase previous PSLG
	_tri_deallocate_all (_tin);
	_tri_deallocate_all (_tou);

	// Points
	_tin.pointlist      = (double*)malloc(NPoints*2*sizeof(double));
	_tin.numberofpoints = NPoints;

	// Segments
	_tin.segmentlist       = (int*)malloc(NSegments*2*sizeof(int));
	_tin.segmentmarkerlist = (int*)malloc(NSegments * sizeof(int));
	_tin.numberofsegments  = NSegments;
	for (size_t i=0; i<NSegments; ++i) _tin.segmentmarkerlist[i]=0;

	// Regions
	if (NRegions>0)
	{
		_tin.regionlist      = (double*)malloc(NRegions*2*sizeof(double));
		_tin.numberofregions = NRegions;
	}

	// Holes
	if (NHoles>0)
	{
		_tin.holelist      = (double*)malloc(NHoles*2*sizeof(double));
		_tin.numberofholes = NHoles;
	}
}

inline size_t Unstructured::Generate()
{
	// Generate (via JRS' triangle)
	triangulate ("pqz", &_tin, &_tou, NULL); // Q=quiet, p=poly, q=quality, z=zero

	/* After triangulate (with -p switch), _tou.regionlist gets the content of _tin.regionlist and
	 * _tou.holelist gets the content of _tin.regionlist. Thus, these output variables must be set
	 * to NULL in order to tell to the destructor to ignore them and then to not double-free memory. */
	_tou.regionlist      = NULL;
	_tou.numberofregions = 0;
	_tou.holelist        = NULL;
	_tou.numberofholes   = 0;

	// Set Vertices
	SetNVerts (_tou.numberofpoints);
	for (int i=0; i<_tou.numberofpoints; ++i)
		SetVert2D (i, false, _tou.pointlist[i*2], _tou.pointlist[i*2+1]);

	// Set Elements
	SetNElems (_tou.numberoftriangles);
	for (int i=0; i<_tou.numberoftriangles; ++i)
	{
		SetElem (i, -1, false, VTK_TRIANGLE, _tou.numberofcorners);
		for (int j=0; j<_tou.numberofcorners; ++j)
			SetElemCon (i, j, _tou.trianglelist[i*_tou.numberofcorners+j]);
	}

	// Return number of elements==triangles
	return _tou.numberoftriangles;
}


/* private */

inline void Unstructured::_tri_set_all_to_null(TriIO & Tio)
{
	// Points
	Tio.pointlist               = NULL;
	Tio.pointattributelist      = NULL;
	Tio.pointmarkerlist         = NULL;
	Tio.numberofpoints          = 0;
	Tio.numberofpointattributes = 0;

	// Unstructureds
	Tio.trianglelist               = NULL;
	Tio.triangleattributelist      = NULL;
	Tio.trianglearealist           = NULL;
	Tio.neighborlist               = NULL;
	Tio.numberoftriangles          = 0;
	Tio.numberofcorners            = 0;
	Tio.numberoftriangleattributes = 0;

	// Segments
	Tio.segmentlist       = NULL;
	Tio.segmentmarkerlist = NULL;
	Tio.numberofsegments  = 0;

	// Holes
	Tio.holelist      = NULL;
	Tio.numberofholes = 0;

	// Regions
	Tio.regionlist      = NULL;
	Tio.numberofregions = 0;

	// Edges
	Tio.edgelist       = NULL;
	Tio.edgemarkerlist = NULL;
	Tio.normlist       = NULL;
	Tio.numberofedges  = 0;
}

inline void Unstructured::_tri_deallocate_all(TriIO & Tio)
{
	// Points
	if (Tio.pointlist          != NULL) free(Tio.pointlist);
	if (Tio.pointattributelist != NULL) free(Tio.pointattributelist);
	if (Tio.pointmarkerlist    != NULL) free(Tio.pointmarkerlist);

	// Triangles
	if (Tio.trianglelist          != NULL) free(Tio.trianglelist);
	if (Tio.triangleattributelist != NULL) free(Tio.triangleattributelist);
	if (Tio.trianglearealist      != NULL) free(Tio.trianglearealist);
	if (Tio.neighborlist          != NULL) free(Tio.neighborlist);

	// Segments
	if (Tio.segmentlist       != NULL) free(Tio.segmentlist);
	if (Tio.segmentmarkerlist != NULL) free(Tio.segmentmarkerlist);

	// Holes
	if (Tio.holelist != NULL) free(Tio.holelist);

	// Regions
	if (Tio.regionlist != NULL) free(Tio.regionlist);

	// Edges
	if (Tio.edgelist       != NULL) free(Tio.edgelist);
	if (Tio.edgemarkerlist != NULL) free(Tio.edgemarkerlist);
	if (Tio.normlist       != NULL) free(Tio.normlist);

	// Clear all
	_tri_set_all_to_null (Tio);
}


#ifdef USE_BOOST_PYTHON
// {

/*
inline size_t Unstructured::PyGenerate(BPy::list const & ListOfMeshBlock)
{
	int nb = BPy::len(ListOfMeshBlock);
	if (nb<1) throw new Fatal("Unstructured::PyGenerate: Number of blocks must be greater than 0 (%d is invalid)",nb);
	Array<Mesh::Block*> blocks;
	blocks.Resize (nb);
	for (int i=0; i<nb; ++i)
		blocks[i] = BPy::extract<Mesh::Block*>(ListOfMeshBlock[i])();
	return Mesh::Unstructured::Generate (blocks);
}
*/

// }
#endif

// Private methods -- Unstructured

inline void Unstructured::_vtk_con(Elem const * E, String & Connect) const
{
	if (_is_3d) // Tetrahedrons
	{
		throw new Fatal("Unstructured::_vtk_con: 3D elements (Tetrahedrons) are not implemented yet");
	}
	else // Triangles
	{
		if (_is_o2) // quadratic
		{
			Connect.Printf("%d %d %d %d %d %d",E->V[0]->MyID,
			                                   E->V[1]->MyID,
			                                   E->V[2]->MyID,
			                                   E->V[3]->MyID,
			                                   E->V[4]->MyID,
			                                   E->V[5]->MyID);
		}
		else // linear
		{
			Connect.Printf("%d %d %d",E->V[0]->MyID,
			                          E->V[1]->MyID,
			                          E->V[2]->MyID);
		}
	}
}

}; // namespace Mesh


#endif // MECHSYS_MESH_UNSTRUCTURED_H
