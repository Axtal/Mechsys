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

#ifndef MECHSYS_MESH_ALPHASHAPE_H
#define MECHSYS_MESH_ALPHASHAPE_H

// STL
#include <list>
#include <map>

// Blitz++
#include <blitz/tinyvec-et.h>

// CGAL
#include <CGAL/cartesian.h>
#include <CGAL/algorithm.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>

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

/** CGAL Structures. */
typedef CGAL::Cartesian<double>                                              CGAL_K;
typedef CGAL::Alpha_shape_vertex_base_2<CGAL_K>                              CGAL_AV;
typedef CGAL::Triangulation_face_base_2<CGAL_K>                              CGAL_TF;
typedef CGAL::Alpha_shape_face_base_2<CGAL_K,CGAL_TF>                        CGAL_AF;
typedef CGAL::Triangulation_default_data_structure_2<CGAL_K,CGAL_AV,CGAL_AF> CGAL_TDS;
typedef CGAL::Delaunay_triangulation_2<CGAL_K,CGAL_TDS>                      CGAL_DT;
typedef CGAL::Alpha_shape_2<CGAL_DT>                                         CGAL_AlphaShape2;

class AlphaShape : public virtual Mesh::Generic
{
public:
	// Set Methods
	void ResetCloud    ();                               ///< Reset cloud of points
	void AddCloudPoint (double X, double Y, double Z=0); ///< Add new point to the list of points i of input PSLG. SetCloudSize MUST be called first.

	// Methods
	size_t Generate (double Alpha=-1);

private:
	// Data
	std::list<CGAL_K::Point_2> _pts_in_2; ///< List of input points (2D)
	std::list<CGAL_K::Point_3> _pts_in_3; ///< List of input points (3D)

	// Overloaded private methods
	void   _vtk_con          (size_t i, String & Connect) const;
	size_t _edge_to_lef_vert (size_t EdgeLocalID)         const { return 0; }
	size_t _edge_to_rig_vert (size_t EdgeLocalID)         const { return 1; }
	void   _face_to_verts    (size_t FaceLocalID, Array<size_t> & Verts) const {}
	void   _face_to_edges    (size_t FaceLocalID, Array<size_t> & Edges) const {}

}; // class AlphaShape


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void AlphaShape::ResetCloud()
{
	_pts_in_2.clear();
	_pts_in_3.clear();
}

inline void AlphaShape::AddCloudPoint(double X, double Y, double Z)
{
	if (_is_3d) _pts_in_3.push_back(CGAL_K::Point_3(X,Y,Z));
	else        _pts_in_2.push_back(CGAL_K::Point_2(X,Y));
}

inline size_t AlphaShape::Generate(double Alpha)
{
	if (_is_3d)
	{
		throw new Fatal("AlphaShape::Generate: 3D: Feature not implemented yet.");
		return 0;
	}
	else
	{
		// Alpha-shape structure
		CGAL_AlphaShape2 as(_pts_in_2.begin(), _pts_in_2.end());

		// Find optimal alpha
		double alp = Alpha;
		if (alp<0) alp = (*as.find_optimal_alpha(1));

		// Generate alpha shape
		as.set_alpha (alp);

		// Erase old mesh
		_erase ();

		// Set Vertices
		size_t id = 0;
		_verts    .Resize (0);
		_verts_bry.Resize (0);
		std::map<CGAL_AlphaShape2::Vertex_handle,size_t> vs;
		for (CGAL_AlphaShape2::Alpha_shape_vertices_iterator it=as.alpha_shape_vertices_begin(); it!=as.alpha_shape_vertices_end(); ++it)
		{
			// Set _verts
			_verts.Push (new Vertex);
			_verts[id]->MyID    = id;
			_verts[id]->OnBry   = true;
			_verts[id]->EdgesID = -1;
			_verts[id]->FacesID = -1;
			_verts[id]->Dupl    = false;
			_verts[id]->C.Resize(2);
			_verts[id]->C = (*it)->point().x(), (*it)->point().y();

			// Set _verts_bry
			_verts_bry.Push (_verts[id]);

			// Set vs
			vs[(*it)] = id++;
		}

		// Set Elements
		id = 0;
		for (size_t i=0; i<_elems.Size(); ++i) if (_elems[i]!=NULL) delete _elems[i];
		_elems    .Resize (0);
		_elems_bry.Resize (0);
		for (CGAL_AlphaShape2::Alpha_shape_edges_iterator it=as.alpha_shape_edges_begin(); it!=as.alpha_shape_edges_end(); ++it)
		{
			// Set _elems
			_elems.Push (new Elem);
			_elems[id]->MyID        = id;
			_elems[id]->Tag         = 0;
			_elems[id]->OnBry       = false;
			_elems[id]->VTKCellType = VTK_LINE;
			_elems[id]->V.Resize (2);

			// Set _elems_bry
			_elems_bry.Push (_elems[id]);

			// Set connectivity
			_elems[id]->V[0] = _verts[vs[it->first->vertex(as.ccw(it->second))]];
			_elems[id]->V[1] = _verts[vs[it->first->vertex(as. cw(it->second))]];

			// Next ID
			id++;
		}
		return _elems.Size();
	}
}


/* private */

inline void AlphaShape::_vtk_con(size_t i, String & Connect) const
{
	if (_is_3d) throw new Fatal("AlphaShape::_vtk_con: 3D unstructured elements are not available (yet).");
	else
	{
		Connect.Printf("%d %d",ElemCon(i,0),
		                       ElemCon(i,1));
	}
}

}; // namespace Mesh


#endif // MECHSYS_MESH_ALPHASHAPE_H
