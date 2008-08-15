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

// Boost
#if defined(USE_BOOST) || defined(USE_BOOST_PYTHON)
  #include <boost/tuple/tuple_io.hpp>
#endif

// MechSys
#include "fem/node.h"
#include "fem/element.h"
#include "fem/geometry.h"
#include "util/array.h"
#include "util/numstreams.h"
#include "util/exception.h"
#include "mesh/mesh.h"
#include "mesh/structured.h"

namespace FEM
{

#if defined(USE_BOOST) || defined(USE_BOOST_PYTHON)

typedef Array< boost::tuple<double,double,double, char const *,double> >               NBrys_T; // Node: x,y,z, key, val
typedef Array< boost::tuple<                 int, char const *,double> >               EBrys_T; // Edge:   tag, key, val
typedef Array< boost::tuple<                 int, char const *,double> >               FBrys_T; // Face:   tag, key, val
typedef Array< boost::tuple<int, char const*, char const*, char const*, char const*> > EAtts_T; // Elem:   tag, type, model, prms, inis

inline void SetGeom (Mesh::Generic const * M,          ///< In: The mesh
                     NBrys_T       const * NodesBrys,  ///< In: Give NULL when there are no nodes boundary conditions
                     EBrys_T       const * EdgesBrys,  ///< In: Give NULL for 3D meshes without edges boundary conditions
                     FBrys_T       const * FacesBrys,  ///< In: Give NULL for 2D meshes
                     EAtts_T       const * ElemsAtts,  ///< In: Elements attributes
                     FEM::Geom           * G,          ///< Out: The FE geometry
                     double                Tol=1.0e-5) ///< In: Tolerance to be used when comparing Nodes
{
	/* Example:
	
		// Nodes brys
		FEM::NBrys_T nbrys;
		nbrys.Push (make_tuple(L/2., 0.0, 0.0, "ux", 0.0)); // x,y,z, key, val

		// Edges brys
		FEM::EBrys_T ebrys;
		ebrys.Push (make_tuple(-10, "uy", 0.0)); // tag, key, val
		ebrys.Push (make_tuple(-20, "fy",  -q)); // tag, key, val

		// Faces brys
		FEM::FBrys_T fbrys;
		fbrys.Push (make_tuple(-100, "uy", 0.0)); // tag, key, val
		fbrys.Push (make_tuple(-200, "fy",  -q)); // tag, key, val

		// Elements attributes
		FEM::EAtts_T eatts;
		eatts.Push (make_tuple(-1, "Quad4PStrain", "LinElastic", "E=207 nu=0.3", "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0")); // tag, type, model, prms, inis
	*/

	// 3D mesh?
	bool is3d = M->Is3D();

	// Set nodes
	size_t nn = M->NVerts();
	G->SetNNodes (nn);
	for (size_t i=0; i<nn; ++i) // loop over all vertices
	{
		// New node
		G->SetNode (i, M->VertX(i), M->VertY(i), (is3d ? M->VertZ(i) : 0.0));
	}

	// Set elements
	size_t ne = M->NElems();
	G->SetNElems (ne);
	for (size_t i=0; i<ne; ++i)
	{
		// Set element
		bool found = false;
		for (size_t j=0; j<ElemsAtts->Size(); ++j)
		{
			if (M->ElemTag(i)==(*ElemsAtts)[j].get<0>())
			{
				// New finite element
				found = true;
				FEM::Element * fe = G->SetElem (i, (*ElemsAtts)[j].get<1>());

				// Set connectivity
				for (size_t k=0; k<M->ElemNVerts(i); ++k)
					fe->Connect (k, G->Nod(M->ElemCon(i,k)));

				// Set parameters and initial values
				fe->SetModel ((*ElemsAtts)[j].get<2>(), (*ElemsAtts)[j].get<3>(), (*ElemsAtts)[j].get<4>());
				break;
			}
		}
		if (found==false) throw new Fatal("SetGeom: Could not find Tag==%d for Element %d in the ElemsAtts list",M->ElemTag(i),i);
	}

	// Set faces boundaries
	if (is3d && FacesBrys!=NULL)
	{
		if (FacesBrys->Size()>0)
		{
			for (size_t b=0; b<M->NElemsBry(); ++b) // loop over all elements on boundary
			{
				int i = M->ElemBry(b);
				for (size_t j=0; j<M->ElemNFTags(i); ++j) // j is the local face id
				{
					int tag = M->ElemFTag(i, j);
					if (tag<0) // this element has a face tag
					{
						bool found = false;
						for (size_t k=0; k<FacesBrys->Size(); ++k)
						{
							if (tag==(*FacesBrys)[k].get<0>())
							{
								found = true;
								G->Ele(i)->FaceBry ((*FacesBrys)[k].get<1>(), (*FacesBrys)[k].get<2>(), j);
								break;
							}
						}
						if (found==false) throw new Fatal("SetGeom: Could not find Tag==%d for Face %d of Element %d in the FacesBrys list",tag,j,i);
					}
				}
			}
		}
	}

	// Set edges boundaries
	if (EdgesBrys!=NULL)
	{
		if (EdgesBrys->Size()>0)
		{
			for (size_t b=0; b<M->NElemsBry(); ++b) // loop over all elements on boundary
			{
				int i = M->ElemBry(b);
				for (size_t j=0; j<M->ElemNETags(i); ++j) // j is the local edge id
				{
					int tag = M->ElemETag(i, j);
					if (tag<0) // this element has an edge tag
					{
						bool found = false;
						for (size_t k=0; k<EdgesBrys->Size(); ++k)
						{
							if (tag==(*EdgesBrys)[k].get<0>())
							{
								found = true;
								G->Ele(i)->EdgeBry ((*EdgesBrys)[k].get<1>(), (*EdgesBrys)[k].get<2>(), j);
								break;
							}
						}
						if (found==false) throw new Fatal("SetGeom: Could not find Tag==%d for Face %d of Element %d in the EdgesBrys list",tag,j,i);
					}
				}
			}
		}
	}

	// Set nodes boundaries
	if (NodesBrys!=NULL)
	{
		if (NodesBrys->Size()>0)
		{
			for (size_t b=0; b<M->NVertsBry(); ++b) // loop over all vertices on boundary
			{
				int i = M->VertBry(b);
				for (size_t j=0; j<NodesBrys->Size(); ++j)
				{
					double x =         (*NodesBrys)[j].get<0>();
					double y =         (*NodesBrys)[j].get<1>();
					double z = (is3d ? (*NodesBrys)[j].get<2>() : 0.0);
					double d = sqrt(pow(x - M->VertX(i),2.0) + pow(y - M->VertY(i),2.0) + (is3d ? pow(z - M->VertZ(i),2.0) : 0.0));
					if (d<Tol) G->Nod(i)->Bry ((*NodesBrys)[j].get<3>(), (*NodesBrys)[j].get<4>());
				}
			}
		}
	}
}

#endif // USE_BOOST || USE_BOOST_PYTHON

}; // namespace FEM

#ifdef USE_BOOST_PYTHON
// {

namespace BPy = boost::python;

void PySetGeom (Mesh::Generic const & M,          ///< In: The mesh
                BPy::list     const & NodesBrys,  ///< In: Give [] when there are no nodes boundary conditions
                BPy::list     const & EdgesBrys,  ///< In: Give [] for 3D mesh without edge boundary conditions
                BPy::list     const & FacesBrys,  ///< In: Give [] for 2D meshes
                BPy::list     const & ElemsAtts,  ///< In: Elements attributes
                FEM::Geom           & G,          ///< Out: The FE geometry
                double                Tol=1.0e-5) ///< In: Tolerance to be used when comparing Nodes
{
	/* Example:
	 *           # Nodes brys
	 *           nbrys = [[L/2., 0.0, 0.0, 'ux', 0.0]] # x,y,z, key, val
	 *
	 *           # Edges brys
	 *           ebrys = [[-10, 'uy', 0.0], # [tag], [key], [val]
	 *                    [-20, 'fy',  -q]] # [tag], [key], [val]
	 *           
	 *           # Faces brys
	 *           fbrys = [[-100, 'uy', 0.0], # [tag], [key], [val]
	 *                    [-200, 'fy',  -q]] # [tag], [key], [val]
	 *           
	 *           # Elements attributes
	 *           eatts = [[-1, 'Quad4PStrain', 'LinElastic', 'E=%f nu=%f'%(E,nu), 'Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0']] # [tag], [type], [model], [prms], [inis]
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

	// Extract list with elements attributes
	FEM::EAtts_T eatts;
	int eatts_size = len(ElemsAtts);
	if (eatts_size>0) eatts.Resize(eatts_size);
	for (int i=0; i<eatts_size; ++i)
	{
		if (len(ElemsAtts[i])==5)
		{
			BPy::list lst = BPy::extract<BPy::list>(ElemsAtts[i])();
			eatts[i] = boost::make_tuple(BPy::extract<int        >(lst[0])(),
			                             BPy::extract<char const*>(lst[1])(),
			                             BPy::extract<char const*>(lst[2])(),
			                             BPy::extract<char const*>(lst[3])(),
			                             BPy::extract<char const*>(lst[4])());
		}
		else throw new Fatal("PySetGeom: Each sublist in ElemsAtts must have 5 items: tag, type, model, prms, inis\n\tExample: ElemsAtts = [[-1, 'Quad4PStrain', 'LinElastic', 'E=207.0 nu=0.3', 'Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0']]");
	}

	// Set geometry
	FEM::SetGeom (&M, nbrys, ebrys, fbrys, &eatts, &G, Tol);

	// Clean up
	if (nbrys!=NULL) delete nbrys;
	if (ebrys!=NULL) delete ebrys;
	if (fbrys!=NULL) delete fbrys;
}

// }
#endif // USE_BOOST_PYTHON

#endif // MECHSYS_FEM_FUNCTIONS_H
