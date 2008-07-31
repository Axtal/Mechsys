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
	Unstructured (double Tol=sqrt(DBL_EPSILON)) : Mesh::Generic(Tol) {} ///< Tol is the tolerance to regard two vertices as coincident

	// Destructor
	~Unstructured () { _erase(); }

	// Methods
	void   SetNVerts (size_t NVerts) ///< Erase any previous mesh and set the number of vertices.
	{
		
	}
	size_t Generate () { return 0; }
	//size_t Generate  (Array<Region*> const & Regions, Array<Hole*> const & Holes); ///< Returns the number of elements. Boundary marks are set first for Faces, then Edges, then Vertices (if any)

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
	void _erase            ();
	int  _edge_to_lef_vert (int EdgeLocalID) const { return Edge2Vert[EdgeLocalID].L; }
	int  _edge_to_rig_vert (int EdgeLocalID) const { return Edge2Vert[EdgeLocalID].R; }

	// Private methods
	void _tri_set_all_to_null (TriIO & Tio)
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

	void _tri_allocate_point_list (size_t nPoints, TriIO & Tio) ///< Allocate only the pointlist array inside Triangle's IO structure
	{
		Tio.pointlist      = (double*)malloc(nPoints*2*sizeof(double));
		Tio.numberofpoints = nPoints;
	}

	void _tri_allocate_edge_list (size_t nEdges, TriIO & Tio) ///< Allocate only the edgelist array inside Triangle's IO structure
	{
		Tio.edgelist      = (int*)malloc(nEdges*2*sizeof(int));
		Tio.numberofedges = nEdges;
	}

	void _tri_deallocate_all (TriIO & Tio) ///< Deallocate all arrays inside Triangle's IO structure
	{
		// Points
		if (Tio.pointlist          != NULL) free(Tio.pointlist);
		if (Tio.pointattributelist != NULL) free(Tio.pointattributelist);
		if (Tio.pointmarkerlist    != NULL) free(Tio.pointmarkerlist);

		// Unstructureds
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

}; // class Unstructured

Edge Unstructured::Edge2Vert[]= {{ 0, 1 },
                                 { 1, 2 },
                                 { 0, 2 }};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


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
		Connect.Printf("%d %d %d %d %d %d %d %d",E->V[1]->MyID,
		                                         E->V[2]->MyID,
		                                         E->V[3]->MyID,
		                                         E->V[0]->MyID,
		                                         E->V[5]->MyID,
		                                         E->V[6]->MyID,
		                                         E->V[7]->MyID,
		                                         E->V[4]->MyID);
	}
	else // Triangles
	{
		if (true) // o2
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
		}
	}
}

inline void Unstructured::_erase()
{
	for (size_t i=0; i<_verts.Size(); ++i) if (_verts[i]!=NULL) delete _verts[i]; // it is only necessary to delete nodes in _verts_d array
	for (size_t i=0; i<_elems.Size(); ++i) if (_elems[i]!=NULL) delete _elems[i]; // it is only necessary to delete elems in _elems array
	_is_3d = false;
	_verts    .Resize(0);
	_verts_bry.Resize(0);
	_verts    .Resize(0);
	_elems    .Resize(0);
	_elems_bry.Resize(0);
	_verts_bry.Resize(0);
	_tri_deallocate_all (_tin);
	_tri_deallocate_all (_tou);
}

}; // namespace Mesh


#endif // MECHSYS_MESH_UNSTRUCTURED_H
