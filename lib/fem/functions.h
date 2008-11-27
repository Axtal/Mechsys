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

#ifndef MECHSYS_FEM_FUNCTIONS_H
#define MECHSYS_FEM_FUNCTIONS_H

// STL
#include <iostream>
#include <fstream>
#include <cfloat> // for DBL_EPSILON
#include <cstring>

// Boost
#include <boost/tuple/tuple_io.hpp>

// MechSys
#include "fem/node.h"
#include "fem/element.h"
#include "fem/geometry.h"
#include "util/array.h"
#include "util/numstreams.h"
#include "util/exception.h"
#include "mesh/mesh.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;

namespace FEM
{

typedef Array< boost::tuple<double,double,double, char const *,double> > NBrys_T; // Node: x,y,z, key, val
typedef Array< boost::tuple<                 int, char const *,double> > EBrys_T; // Edge:   tag, key, val
typedef Array< boost::tuple<                 int, char const *,double> > FBrys_T; // Face:   tag, key, val
typedef Array< boost::tuple<int, char const*, char const*, char const*,
                            char const*, char const*, bool> > EAtts_T; // Elem: tag, type, model, prms, inis, props, active

inline void SetNodesElems (Mesh::Generic const * M,          ///< In: The mesh
                           EAtts_T       const * ElemsAtts,  ///< In: Elements attributes
                           FEM::Geom           * G,          ///< Out: The FE geometry
                           double                Tol=1.0e-5,      ///< In: Tolerance for defining whether X-Y-Z coordinates are in a same plane or not
                           bool                  OnlyFrame=false) ///< Only frame (beam/truss) structure ?
{
	/* Example:
	
		// Elements attributes
		FEM::EAtts_T eatts;
		eatts.Push (make_tuple(-1, "Quad4PStrain", "LinElastic", "E=207 nu=0.3", "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0", "gam=20")); // tag, type, model, prms, inis, props
	*/

	// 3D mesh?
	bool is3d = M->Is3D();

	// Set nodes
	size_t nn = M->NVerts();
	G->SetNNodes (nn);
	double diff_x = 0.0; // used for verification (try to find which plane the problem is defined in)
	double diff_y = 0.0;
	double diff_z = 0.0;
	for (size_t i=0; i<nn; ++i) // loop over all vertices
	{
		// New node
		G->SetNode (i, M->VertX(i), M->VertY(i), (is3d ? M->VertZ(i) : 0.0));
		          diff_x += fabs(M->VertX(0)-M->VertX(i));
		          diff_y += fabs(M->VertY(0)-M->VertY(i));
		if (is3d) diff_z += fabs(M->VertZ(0)-M->VertZ(i));
	}

	// Check working plane
	if (is3d)
	{
		if (diff_x<Tol || diff_y<Tol || diff_z<Tol)
			throw new Fatal("FEM::SetNodesElems: For 3D problems, vertices cannot be all in a same plane (diff_x=%f, diff_y=%f, diff_z=%f)",diff_x,diff_y,diff_z);
	}
	else
	{
		if (diff_z>Tol)
			throw new Fatal("FEM::SetNodesElems: For 2D problems, only the X and Y coordinates must be used (diff_z=%f)",diff_z);
	}

	// Number of beams (with duplicates)
	Array< boost::tuple<size_t,size_t,size_t,int,int,int> > beams; // elem_id, eatt_id, local_edge_id, beam_tag, v0, v1
	if (OnlyFrame==false)
	{
		for (size_t k=0; k<ElemsAtts->Size(); ++k)
		{
			if (strcmp((*ElemsAtts)[k].get<1>(),"Beam")==0)
			{
				int beam_edge_tag = (*ElemsAtts)[k].get<0>();
				for (size_t i=0; i<M->NElems(); ++i)
				{
					for (size_t j=0; j<M->ElemNETags(i); ++j)
					{
						if (M->ElemETag(i,j)==beam_edge_tag)
						{
							int v0 = M->EdgeToLef(i, j);
							int v1 = M->EdgeToRig(i, j);
							bool is_new = true;
							for (size_t m=0; m<beams.Size(); ++m)
							{
								int w0 = beams[m].get<4>();
								int w1 = beams[m].get<5>();
								if ((v0==w0 || v0==w1) && (v1==w0 || v1==w1)) // coincident
								{
									is_new = false;
									break;
								}
							}
							if (is_new) beams.Push (boost::make_tuple(i, k, j, beam_edge_tag, v0, v1));
						}
					}
				}
			}
		}
		G->SetNBeams (beams.Size());
	}

	// Set elements
	G->SetNElems (M->NElems() + beams.Size());
	for (size_t i=0; i<M->NElems(); ++i)
	{
		// Set element
		bool found = false;
		for (size_t j=0; j<ElemsAtts->Size(); ++j)
		{
			if (M->ElemTag(i)==(*ElemsAtts)[j].get<0>())
			{
				// New finite element
				found = true;
				FEM::Element * fe = G->SetElem (i, (*ElemsAtts)[j].get<1>(), (*ElemsAtts)[j].get<6>(), M->ElemTag(i));

				// Set connectivity
				for (size_t k=0; k<M->ElemNVerts(i); ++k)
					fe->Connect (k, G->Nod(M->ElemCon(i,k)));

				// Set parameters and initial values
				fe->SetModel ((*ElemsAtts)[j].get<2>(), (*ElemsAtts)[j].get<3>(), (*ElemsAtts)[j].get<4>());

				// Set properties
				fe->SetProps ((*ElemsAtts)[j].get<5>());
				break;
			}
		}
		if (found==false) throw new Fatal("SetGeom: Could not find Tag==%d for Element %d in the ElemsAtts list",M->ElemTag(i),i);
	}

	// Set beams
	size_t ie = M->NElems();
	for (size_t i=0; i<beams.Size(); ++i)
	{
		// Data
		size_t elem_id       = beams[i].get<0>();
		size_t eatt_id       = beams[i].get<1>();
		size_t local_edge_id = beams[i].get<2>();
		int    beam_tag      = beams[i].get<3>();

		// New finite element
		FEM::Element * fe = G->SetElem (ie, (*ElemsAtts)[eatt_id].get<1>(), (*ElemsAtts)[eatt_id].get<6>(), beam_tag);

		// Set connectivity
		fe->Connect (0, G->Nod(M->EdgeToLef(elem_id, local_edge_id)));
		fe->Connect (1, G->Nod(M->EdgeToRig(elem_id, local_edge_id)));

		// Set parameters and initial values
		fe->SetModel ((*ElemsAtts)[eatt_id].get<2>(), (*ElemsAtts)[eatt_id].get<3>(), (*ElemsAtts)[eatt_id].get<4>());

		// Set properties
		fe->SetProps ((*ElemsAtts)[eatt_id].get<5>());
		ie++;

		// Set beam
		G->SetBeam (i, fe, beam_tag);
	}
}

inline void SetBrys (Mesh::Generic const * M,          ///< In: The mesh
                     NBrys_T       const * NodesBrys,  ///< In: Give NULL when there are no nodes boundary conditions
                     EBrys_T       const * EdgesBrys,  ///< In: Give NULL for 3D meshes without edges boundary conditions
                     FBrys_T       const * FacesBrys,  ///< In: Give NULL for 2D meshes
                     FEM::Geom           * G,          ///< Out: The FE geometry
                     double                Tol=1.0e-5) ///< In: Tolerance to be used when comparing Nodes
{
	/* Example:
	
		// Nodes brys
		FEM::NBrys_T nbrys;
		nbrys.Push (make_tuple(L/2., 0.0, 0.0, "ux", 0.0)); // x,y,z, key, val

		// Edges brys (the order matters!)
		FEM::EBrys_T ebrys;
		ebrys.Push (make_tuple(-10, "uy", 0.0)); // tag, key, val
		ebrys.Push (make_tuple(-20, "fy",  -q)); // tag, key, val

		// Faces brys (the order matters!)
		FEM::FBrys_T fbrys;
		fbrys.Push (make_tuple(-100, "uy", 0.0)); // tag, key, val
		fbrys.Push (make_tuple(-200, "fy",  -q)); // tag, key, val

	*/

	// Erase previous boundary conditions
	for (size_t i=0; i<G->NNodes(); ++i) G->Nod(i)->ClearBryValues();

	// 3D mesh?
	bool is3d = M->Is3D();

	// Set faces boundaries (the order matters)
	if (is3d && FacesBrys!=NULL)
	{
		for (size_t k=0; k<FacesBrys->Size(); ++k)
		{
			for (size_t b=0; b<M->NElemsBry(); ++b) // loop over all elements on boundary
			{
				int i = M->ElemBry(b);
				for (size_t j=0; j<M->ElemNFTags(i); ++j) // j is the local face id
				{
					int tag = M->ElemFTag(i, j);
					if (tag<0) // this element has a face tag
					{
						if (tag==(*FacesBrys)[k].get<0>())
						{
							G->Ele(i)->FaceBry ((*FacesBrys)[k].get<1>(), (*FacesBrys)[k].get<2>(), j);
							break; // go to the next element on boundary
						}
					}
				}
			}
		}
	}

	// Set edges boundaries (the order matters)
	if (EdgesBrys!=NULL)
	{
		for (size_t k=0; k<EdgesBrys->Size(); ++k)
		{
			for (size_t b=0; b<M->NElemsBry(); ++b) // loop over all elements on boundary
			{
				int i = M->ElemBry(b);
				for (size_t j=0; j<M->ElemNETags(i); ++j) // j is the local edge id
				{
					int tag = M->ElemETag(i, j);
					if (tag<0) // this element has an edge tag
					{
						if (tag==(*EdgesBrys)[k].get<0>())
						{
							G->Ele(i)->EdgeBry ((*EdgesBrys)[k].get<1>(), (*EdgesBrys)[k].get<2>(), j);
							break; // go to the next element on boundary
						}
					}
				}
			}
			for (size_t b=0; b<G->NBeams(); ++b)
			{
				if (G->BTag(b)==(*EdgesBrys)[k].get<0>())
					G->Beam(b)->EdgeBry ((*EdgesBrys)[k].get<1>(), (*EdgesBrys)[k].get<2>(), 0);
			}
		}
	}

	// Set nodes boundaries
	if (NodesBrys!=NULL)
	{
		for (size_t j=0; j<NodesBrys->Size(); ++j)
		{
			for (size_t b=0; b<M->NVertsBry(); ++b) // loop over all vertices on boundary
			{
				int i = M->VertBry(b);
				double x =         (*NodesBrys)[j].get<0>();
				double y =         (*NodesBrys)[j].get<1>();
				double z = (is3d ? (*NodesBrys)[j].get<2>() : 0.0);
				double d = sqrt(pow(x - M->VertX(i),2.0) + pow(y - M->VertY(i),2.0) + (is3d ? pow(z - M->VertZ(i),2.0) : 0.0));
				if (d<Tol) G->Nod(i)->Bry ((*NodesBrys)[j].get<3>(), (*NodesBrys)[j].get<4>());
			}
		}
	}
}

#endif // USE_BOOST || USE_BOOST_PYTHON

}; // namespace FEM

#ifdef USE_BOOST_PYTHON
// {

namespace BPy = boost::python;

void PySetNodesElems (Mesh::Generic const & M,          ///< In: The mesh
                      BPy::list     const & ElemsAtts,  ///< In: Elements attributes
                      FEM::Geom           & G,          ///< Out: The FE geometry
                      double                Tol=1.0e-5,      ///< In: Tolerance for defining whether X-Y-Z coordinates are in a same plane or not
                      bool                  OnlyFrame=false) ///< Only frame (beam/truss) structure ?
{
	/* Example:
	 *           
	 *           # Elements attributes
	 *           eatts = [[-1, 'Quad4PStrain', 'LinElastic', 'E=%f nu=%f'%(E,nu), 'Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0', 'gam=20']] # tag, type, model, prms, inis, props
	 */

	// Extract list with elements attributes
	FEM::EAtts_T eatts;
	int eatts_size = len(ElemsAtts);
	if (eatts_size>0) eatts.Resize(eatts_size);
	for (int i=0; i<eatts_size; ++i)
	{
		if (len(ElemsAtts[i])==6)
		{
			BPy::list lst = BPy::extract<BPy::list>(ElemsAtts[i])();
			eatts[i] = boost::make_tuple(BPy::extract<int        >(lst[0])(),
			                             BPy::extract<char const*>(lst[1])(),
			                             BPy::extract<char const*>(lst[2])(),
			                             BPy::extract<char const*>(lst[3])(),
			                             BPy::extract<char const*>(lst[4])(),
			                             BPy::extract<char const*>(lst[5])());
		}
		else throw new Fatal("PySetNodesElems: Each sublist in ElemsAtts must have 6 items: tag, type, model, prms, inis, props\n\tExample: ElemsAtts = [[-1, 'Quad4PStrain', 'LinElastic', 'E=207.0 nu=0.3', 'Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0', 'gam=20']]");
	}

	// Set geometry
	FEM::SetNodesElems (&M, &eatts, &G, Tol, OnlyFrame);
}

void PySetBrys (Mesh::Generic const & M,          ///< In: The mesh
                BPy::list     const & NodesBrys,  ///< In: Give [] when there are no nodes boundary conditions
                BPy::list     const & EdgesBrys,  ///< In: Give [] for 3D mesh without edge boundary conditions
                BPy::list     const & FacesBrys,  ///< In: Give [] for 2D meshes
                FEM::Geom           & G,          ///< Out: The FE geometry
                double                Tol=1.0e-5) ///< In: Tolerance to be used when comparing Nodes
{
	/* Example:
	 *           # Nodes brys
	 *           nbrys = [[L/2., 0.0, 0.0, 'ux', 0.0]] # x,y,z, key, val
	 *
	 *           # Edges brys (the order matters!)
	 *           ebrys = [[-10, 'uy', 0.0], # [tag], [key], [val]
	 *                    [-20, 'fy',  -q]] # [tag], [key], [val]
	 *           
	 *           # Faces brys (the order matters!)
	 *           fbrys = [[-100, 'uy', 0.0], # [tag], [key], [val]
	 *                    [-200, 'fy',  -q]] # [tag], [key], [val]
	 */

	// Extract list with nodes boundaries
	int nbrys_size = len(NodesBrys);
	FEM::NBrys_T * nbrys = (nbrys_size>0 ? new FEM::NBrys_T : NULL);
	if (nbrys!=NULL) nbrys->Resize(nbrys_size);
	for (int i=0; i<nbrys_size; ++i)
	{
		if (len(NodesBrys[i])==5)
		{
			BPy::list lst = BPy::extract<BPy::list>(NodesBrys[i])();
			(*nbrys)[i] = boost::make_tuple(BPy::extract<double     >(lst[0])(),
			                                BPy::extract<double     >(lst[1])(),
			                                BPy::extract<double     >(lst[2])(),
			                                BPy::extract<char const*>(lst[3])(),
			                                BPy::extract<double     >(lst[4])());
		}
		else throw new Fatal("PySetGeom: Each sublist in NodesBrys must have 5 items: x,y,z, key, val\n\tExample: NodesBrys = [[1.0, 0.0, 0.0, 'ux', 0.0]]");
	}

	// Extract list with edges boundaries
	int ebrys_size = len(EdgesBrys);
	FEM::EBrys_T * ebrys = (ebrys_size>0 ? new FEM::EBrys_T : NULL);
	if (ebrys!=NULL) ebrys->Resize(ebrys_size);
	for (int i=0; i<ebrys_size; ++i)
	{
		if (len(EdgesBrys[i])==3)
		{
			BPy::list lst = BPy::extract<BPy::list>(EdgesBrys[i])();
			(*ebrys)[i] = boost::make_tuple(BPy::extract<int        >(lst[0])(),
			                                BPy::extract<char const*>(lst[1])(),
			                                BPy::extract<double     >(lst[2])());
		}
		else throw new Fatal("PySetGeom: Each sublist in EdgesBrys must have 3 items: tag, key, val\n\tExample: EdgesBrys = [[-10, 'uy', 0.0], [-20, 'fy', -1]]");
	}

	// Extract list with faces boundaries
	int fbrys_size = len(FacesBrys);
	FEM::FBrys_T * fbrys = (fbrys_size>0 ? new FEM::FBrys_T : NULL);
	if (fbrys!=NULL) fbrys->Resize(fbrys_size);
	for (int i=0; i<fbrys_size; ++i)
	{
		if (len(FacesBrys[i])==3)
		{
			BPy::list lst = BPy::extract<BPy::list>(FacesBrys[i])();
			(*fbrys)[i] = boost::make_tuple(BPy::extract<int        >(lst[0])(),
			                                BPy::extract<char const*>(lst[1])(),
			                                BPy::extract<double     >(lst[2])());
		}
		else throw new Fatal("PySetGeom: Each sublist in FacesBrys must have 3 items: tag, key, val\n\tExample: FacesBrys = [[-10, 'uy', 0.0], [-20, 'fy', -1]]");
	}

	// Set geometry
	FEM::SetBrys (&M, nbrys, ebrys, fbrys, &G, Tol);

	// Clean up
	if (nbrys!=NULL) delete nbrys;
	if (ebrys!=NULL) delete ebrys;
	if (fbrys!=NULL) delete fbrys;
}

#endif // MECHSYS_FEM_FUNCTIONS_H
