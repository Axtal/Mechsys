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

/** JRS' Triangle Input/Output structure. */
typedef triangulateio TriIO;

class Unstructured : public virtual Mesh::Generic
{
public:
	// Constants
	static Edge Edge2Vert[]; ///< Map from local edge ID to local vertex ID
	static Face Face2Vert[]; ///< Map from local face ID to local vertex ID
	static Face Face2Edge[]; ///< Map from local face ID to local edge ID

	// Constructor
	Unstructured ();

	// Destructor
	~Unstructured () { _tri_deallocate_all(_tin); _tri_deallocate_all(_tou); }

	// Set Methods
	void SetPolySize    (size_t NPoints, size_t NSegments, size_t NRegions=0, size_t NHoles=0); ///< Erase any previous input PSLG and set the number of points and segments of the polygon. Also set the number of holes and regions.
	void SetPolyPoint   (size_t i, double X, double Y, double Z=0);                             ///< Set the coordinates of point i of input PSLG. SetPolySize MUST be called first.
	void SetPolySegment (size_t i, size_t iPointLeft, size_t iPointRight, int Tag=0);           ///< Set the left and right points of a segment i of input PSLG. SetPolySize MUST be called first.
	void SetPolyRegion  (size_t i, int Tag, double MaxArea, double X, double Y, double Z=0);    ///< Set the coordinates of a point defining a region i.
	void SetPolyHole    (size_t i, double X, double Y, double Z=0);                             ///< Set the coordinates of a point defining a hole i.

	// Get methods -- derived
	size_t NVerts          ()                   const { return _tou.numberofpoints; }
	size_t NVertsBry       ()                   const { return _vbry.Size(); }
	size_t NElems          ()                   const { return _tou.numberoftriangles; }
	size_t NElemsBry       ()                   const { return NElems(); } ///< TODO
	long   VertBry         (size_t i)           const { return _vbry[i]; }
	long   ElemBry         (size_t i)           const { return i; } ///< TODO
	bool   IsVertOnBry     (size_t i)           const { return (_vbry.Find(i)<0 ? false : true); }
	double VertX           (size_t i)           const { return _tou.pointlist[i*2  ]; }
	double VertY           (size_t i)           const { return _tou.pointlist[i*2+1]; }
	double VertZ           (size_t i)           const { return 0; }
	int    ElemTag         (size_t i)           const { return (_tou.numberoftriangleattributes>0 ? _tou.triangleattributelist[i*_tou.numberoftriangleattributes] : -1); }
	bool   IsElemOnBry     (size_t i)           const { return (_tou.triedgemarks[i*3]+_tou.triedgemarks[i*3+1]+_tou.triedgemarks[i*3+2]==0 ? false : true); }
	int    ElemVTKCellType (size_t i)           const { return (_is_o2?VTK_QUADRATIC_TRIANGLE:VTK_TRIANGLE); }
	size_t ElemNVerts      (size_t i)           const { return _tou.numberofcorners; }
	size_t ElemCon         (size_t i, size_t j) const { return _tou.trianglelist[i*_tou.numberofcorners+j]; }
	size_t ElemNETags      (size_t i)           const { return 3; }
	size_t ElemNFTags      (size_t i)           const { return 0; }
	int    ElemETag        (size_t i, size_t j) const { return _tou.triedgemarks[i*3+j]; }

	// Methods
	size_t Generate (double MaxAreaGlobal=-1, double MinAngle=-1); ///< [opt] MaxAreaGlobal = Uniform maximum area constraint, [opt] MinAngle = Minium angle constraint

private:
	// Data
	TriIO        _tin;  ///< Triangle IO's input structure
	TriIO        _tou;  ///< Triangle IO's output structure
	Array<long>  _vbry; ///< Array with IDs of vertices on boundary

	// Overloaded private methods
	void   _vtk_con          (size_t i, String & Connect) const;
	size_t _edge_to_lef_vert (size_t EdgeLocalID)         const { return Edge2Vert[EdgeLocalID].L; }
	size_t _edge_to_rig_vert (size_t EdgeLocalID)         const { return Edge2Vert[EdgeLocalID].R; }
	void   _face_to_verts    (size_t FaceLocalID, Array<size_t> & Verts) const;
	void   _face_to_edges    (size_t FaceLocalID, Array<size_t> & Edges) const;

	// Private methods
	void _tri_set_all_to_null (TriIO & Tio); ///< Set all elements of Triangle's IO structure to NULL and 0
	void _tri_deallocate_all  (TriIO & Tio); ///< Deallocate all arrays inside Triangle's IO structure

}; // class Unstructured

Edge Unstructured::Edge2Vert[]= {{ 0, 1 },
                                 { 1, 2 },
                                 { 0, 2 }};

Face Unstructured::Face2Vert[]= {{  0,  2,  3,  6,  9,  7 },  // Face # 0 => Vertices 0,2,3...
                                 {  0,  1,  3,  4,  8,  7 },  // Face # 1
                                 {  0,  1,  2,  4,  5,  6 },  // Face # 2
                                 {  1,  2,  3,  5,  9,  8 }}; // Face # 3

Face Unstructured::Face2Edge[]= {{  0,  3,  5,  6,  9, 11 },  // Face # 0 => Edges 0,4,8
                                 {  1,  3,  4,  7,  9, 10 },  // Face # 1
                                 {  0,  1,  2,  6,  7,  8 },  // Face # 2
                                 {  2,  4,  5,  8, 10, 11 }}; // Face # 3


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline Unstructured::Unstructured()
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
		_tin.regionlist      = (double*)malloc(NRegions*4*sizeof(double));
		_tin.numberofregions = NRegions;
	}

	// Holes
	if (NHoles>0)
	{
		_tin.holelist      = (double*)malloc(NHoles*2*sizeof(double));
		_tin.numberofholes = NHoles;
	}
}

inline void Unstructured::SetPolyPoint(size_t i, double X, double Y, double Z)
{
	_tin.pointlist[i*2  ] = X;
	_tin.pointlist[i*2+1] = Y;
}

inline void Unstructured::SetPolySegment(size_t i, size_t iPointLeft, size_t iPointRight, int Tag)
{
	_tin.segmentlist[i*2  ]   = iPointLeft;
	_tin.segmentlist[i*2+1]   = iPointRight;
	_tin.segmentmarkerlist[i] = Tag;
}

inline void Unstructured::SetPolyRegion(size_t i, int Tag, double MaxArea, double X, double Y, double Z)
{
	_tin.regionlist[i*4  ] = X;
	_tin.regionlist[i*4+1] = Y;
	_tin.regionlist[i*4+2] = Tag;
	_tin.regionlist[i*4+3] = MaxArea;
}

inline void Unstructured::SetPolyHole(size_t i, double X, double Y, double Z)
{
	_tin.holelist[i*2  ] = X;
	_tin.holelist[i*2+1] = Y;
}

inline size_t Unstructured::Generate(double MaxAreaGlobal, double MinAngle)
{
	// Generate (via JRS' triangle)
	String prms("pzA"); // Q=quiet, p=poly, q=quality, z=zero
	if (MaxAreaGlobal>0) prms.Printf("%sa%f", prms.CStr(), MaxAreaGlobal);
	if (MinAngle     >0) prms.Printf("%sq%f", prms.CStr(), MinAngle);
	else                 prms.Printf("%sq",   prms.CStr());
	prms.Printf("%sa", prms.CStr());
	std::cout << "JRS' triangle parameters = " << prms << std::endl;
	triangulate (prms.CStr(), &_tin, &_tou, NULL);

	/* Find vertices on boundary
	 * _tou.pointmarkerlist[ipoint] will be equal to:
	 * == edgeTag (<0) => on edge with tag
	 * == 0            => internal vertex (not on boundary)
	 * == 1            => on boundary
	 */
	_vbry.Resize(0);
	for (int i=0; i<_tou.numberofpoints; ++i)
	{
		int mark = _tou.pointmarkerlist[i];
		if (mark<0 || mark==1) _vbry.Push(i);
		/*
		for (int j=0; j<_tou.numberofpointattributes; j++)
			std::cout << _tou.pointattributelist[i*_tou.numberofpointattributes+j] << "  ";
		std::cout << std::endl;
		*/
	}

	/*
	for (int i=0; i<_tou.numberoftriangles; ++i)
	{
	   for (int j=0; j<_tou.numberoftriangleattributes; j++)
		std::cout << _tou.triangleattributelist[i*_tou.numberoftriangleattributes+j] << "  ";
	   std::cout << std::endl;
	}
	*/

	/* After triangulate (with -p switch), _tou.regionlist gets the content of _tin.regionlist and
	 * _tou.holelist gets the content of _tin.regionlist. Thus, these output variables must be set
	 * to NULL in order to tell to the destructor to ignore them and then to not double-free memory. */
	_tou.regionlist      = NULL;
	_tou.numberofregions = 0;
	_tou.holelist        = NULL;
	_tou.numberofholes   = 0;

	/*
	// Set Vertices
	SetNVerts (_tou.numberofpoints);
	for (int i=0; i<_tou.numberofpoints; ++i)
		SetVert2D (i, false, _tou.pointlist[i*2], _tou.pointlist[i*2+1]);

	// Set Elements
	SetNElems (_tou.numberoftriangles);
	for (int i=0; i<_tou.numberoftriangles; ++i)
	{
		bool onbry = (_tou.triedgemarks[i*3]+_tou.triedgemarks[i*3+1]+_tou.triedgemarks[i*3+2]==0 ? false : true);
		SetElem (i, -1, onbry, VTK_TRIANGLE, _tou.numberofcorners);
		// Connectivity
		for (int j=0; j<_tou.numberofcorners; ++j)
			SetElemCon (i, j, _tou.trianglelist[i*_tou.numberofcorners+j]);
		// ETags
		SetElemETag (i, 0, _tou.triedgemarks[i*3]  );
		SetElemETag (i, 1, _tou.triedgemarks[i*3+1]);
		SetElemETag (i, 2, _tou.triedgemarks[i*3+2]);
	}
	*/

	// Return number of elements==triangles
	return _tou.numberoftriangles;
}


/* private */

inline void Unstructured::_face_to_verts(size_t FaceLocalID, Array<size_t> & Verts) const
{
	if (_is_o2)
	{
		Verts.Resize(6);
		Verts[0] = Face2Vert[FaceLocalID].I0;
		Verts[1] = Face2Vert[FaceLocalID].I1;
		Verts[2] = Face2Vert[FaceLocalID].I2;
		Verts[3] = Face2Vert[FaceLocalID].I3;
		Verts[4] = Face2Vert[FaceLocalID].I4;
		Verts[5] = Face2Vert[FaceLocalID].I5;
	}
	else
	{
		Verts.Resize(3);
		Verts[0] = Face2Vert[FaceLocalID].I0;
		Verts[1] = Face2Vert[FaceLocalID].I1;
		Verts[2] = Face2Vert[FaceLocalID].I2;
	}
}

inline void Unstructured::_face_to_edges(size_t FaceLocalID, Array<size_t> & Edges) const
{
	if (_is_o2)
	{
		Edges.Resize(6);
		Edges[0] = Face2Edge[FaceLocalID].I0;
		Edges[1] = Face2Edge[FaceLocalID].I1;
		Edges[2] = Face2Edge[FaceLocalID].I2;
		Edges[3] = Face2Edge[FaceLocalID].I3;
		Edges[4] = Face2Edge[FaceLocalID].I4;
		Edges[5] = Face2Edge[FaceLocalID].I5;
	}
	else
	{
		Edges.Resize(3);
		Edges[0] = Face2Edge[FaceLocalID].I0;
		Edges[1] = Face2Edge[FaceLocalID].I1;
		Edges[2] = Face2Edge[FaceLocalID].I2;
	}
}

inline void Unstructured::_tri_set_all_to_null(TriIO & Tio)
{
	// Points
	Tio.pointlist               = NULL;
	Tio.pointattributelist      = NULL;
	Tio.pointmarkerlist         = NULL;
	Tio.numberofpoints          = 0;
	Tio.numberofpointattributes = 0;

	// Triangles
	Tio.trianglelist               = NULL;
	Tio.triangleattributelist      = NULL;
	Tio.trianglearealist           = NULL;
	Tio.neighborlist               = NULL;
	Tio.numberoftriangles          = 0;
	Tio.numberofcorners            = 0;
	Tio.numberoftriangleattributes = 0;
	Tio.triedgemarks               = NULL;

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
	if (Tio.triedgemarks          != NULL) free(Tio.triedgemarks);

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


}; // namespace Mesh


#endif // MECHSYS_MESH_UNSTRUCTURED_H
