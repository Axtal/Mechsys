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

#ifndef MECHSYS_FEM_DATA_H
#define MECHSYS_FEM_DATA_H

// STL
#include <iostream>
#include <cstring>
#include <cfloat> // for DBL_EPSILON
#include <map>

// Boost
#include <boost/tuple/tuple_io.hpp>

// Boost::Python
#ifdef USE_BOOST_PYTHON
  #include <boost/python.hpp> // this includes everything
  namespace BPy = boost::python;
#endif

// MechSys
#include "fem/node.h"
#include "fem/element.h"
#include "mesh/mesh.h"
#include "util/array.h"
#include "models/model.h"

namespace FEM
{

typedef std::map<String,double> Prop_t; ///< Properties type
typedef char const            * Str_t;

// Tuples used when initializing nodes and elements by groups
typedef Array< boost::tuple<double,double,double, Str_t,double> >           NBrys_T; // Node: x,y,z, key, val
typedef Array< boost::tuple<                 int, Str_t,double> >           EBrys_T; // Edge:   tag, key, val
typedef Array< boost::tuple<                 int, Str_t,double> >           FBrys_T; // Face:   tag, key, val
typedef Array< boost::tuple<int,Str_t,Str_t,Str_t,Str_t,Str_t,Str_t,bool> > EAtts_T; // Elem: tag, geom_t, prob_t, model, prms, inis, props, active

/* FEM Data. */
class Data
{
public:
	/* Constructor */
	Data (int nDim) : _dim(nDim), _tol(1.0e-5), _frame(false) {}

	/* Destructor */
	~Data ();

	// Set by groups methods
	void SetTol        (double Tol)      { _tol   = Tol;   } ///< Tolerance for defining whether X-Y-Z coordinates are in a same plane or not. Also used for comparing coincident nodes
	void SetOnlyFrame  (bool Frame=true) { _frame = Frame; } ///< Only frame (beam/truss) structure ?
	void SetNodesElems (Mesh::Generic const * M,             ///< The mesh
	                    EAtts_T       const * ElemsAtts);    ///< Elements attributes
	void SetBrys       (Mesh::Generic const * M,             ///< The mesh
	                    NBrys_T       const * NodesBrys,     ///< Give NULL when there are no nodes boundary conditions
	                    EBrys_T       const * EdgesBrys,     ///< Give NULL for 3D meshes without edges boundary conditions
	                    FBrys_T       const * FacesBrys);    ///< Give NULL for 2D meshes

	// Set methods
	void      SetNNodes       (size_t NNodes);                                         ///< Set the number of nodes
	void      SetNElems       (size_t NElems);                                         ///< Set the number of elements
	Node    * SetNode         (size_t i, double X, double Y, double Z=0.0, int Tag=0); ///< Set a node
	Element * SetElemAndModel (size_t               i,                                 ///< The ID(index) of the element
	                           Array<long>  const & Conn,                              ///< Connectivity
	                           int                  Tag,                               ///< The tag of new element
	                           Str_t                GeomT,                             ///< Geometry Element type. Ex: Tri3
	                           Str_t                ProbT,                             ///< Problem  Element type. Ex: PStrain
	                           Str_t                MdlName,                           ///< Constitutive model. Ex: "LinElastic"
	                           Str_t                Prms,                              ///< Parameters. Ex: "E=100 nu=0.3"
	                           Str_t                Inis,                              ///< Initial state. Ex: Sx=0.0
	                           Str_t                Props,                             ///< Element properties. Ex: gam=20
	                           bool                 IsAct);                            ///< Active/Inactive?

	// Push new Node and Element
	Node    * PushNode (double X, double Y, double Z, int Tag) { _nodes.Push(new Node);  return SetNode(_nodes.Size()-1,X,Y,Z,Tag); }
	Element * PushElem (Array<long> const & Conn,      ///< Array wih IDs of connected nodes
	                    int                 Tag,       ///< A tag for the new element
	                    Str_t               GeomT,     ///< Geometry type. Ex.: Hex8, Beam, Spring
	                    Str_t               ProbT,     ///< Problem type. Ex.: Equilib, Diffusion
	                    Str_t               MdlName,   ///< Constitutive model name
	                    Str_t               Prms,      ///< Parameters. Ex.: E=100 nu=0.3
	                    Str_t               Inis,      ///< Initial values. Ex.: ZERO, Sx=0
	                    Str_t               Props,     ///< Properties. Ex.: gam=20
	                    bool                IsAct);    ///< Is active?

	// Specific methods
	void AddVolForces () { for (size_t i=0; i<_elems.Size(); ++i) _elems[i]->AddVolForces(); } ///< Apply body forces (equilibrium/coupled problems)
	void ClearDisp    ();                                                                      ///< Clear displacements (equilibrium/coupled problems)
	void Activate     (int ElemTag);                                                           ///< Activate all elements with Tag
	void Deactivate   (int ElemTag);                                                           ///< Activate all elements with Tag

	// Beams
	void      SetNBeams (size_t NBeams) { _beams.Resize(NBeams); _beams.SetValues(NULL); _btags.Resize(NBeams); } ///< Set the number of beams
	void      SetBeam   (size_t iBeam, Element * Beam, int Tag) { _beams[iBeam]=Beam; _btags[iBeam]=Tag; }        ///< Set a beam element
	size_t    NBeams    () const       { return _beams.Size(); }                                                  ///< Return the number of beams
	Element * Beam      (size_t iBeam) { return _beams[iBeam]; }                                                  ///< Access a beam element
	int       BTag      (size_t iBeam) { return _btags[iBeam]; }                                                  ///< Return the tag of a beam element

	// Access methods
	bool                    Check        ();                                                  ///< Check if Nodes and Elements were allocated properly. Should be called before accessing Nodes and Elements, since these may not had been allocated yet (and then causing Segfaults).
	size_t                  NNodes       ()         const { return _nodes.Size(); }           ///< Return the number of nodes
	size_t                  NElems       ()         const { return _elems.Size(); }           ///< Return the number of elements
	size_t                  NDim         ()         const { return _dim;          }           ///< Return the dimension
	Node                  * Nod          (size_t i)       { return _nodes[i];     }           ///< Access (read/write) a node
	Element               * Ele          (size_t i)       { return _elems[i];     }           ///< Access (read/write) an element
	Node            const * Nod          (size_t i) const { return _nodes[i];     }           ///< Access (read-only) a node
	Element         const * Ele          (size_t i) const { return _elems[i];     }           ///< Access (read-only) an element
	size_t                  GetNode      (double X, double Y, double Z=0.0);                  ///< Returns the node ID that matchs the specified coordinates
	Array<Node*>          & Nodes        ()               { return _nodes;        }           ///< Access all nodes (read/write)
	Array<Element*>       & Elems        ()               { return _elems;        }           ///< Access all elements (read/write)
	Array<Node*>    const & Nodes        ()         const { return _nodes;        }           ///< Access all nodes (read-only)
	Array<Element*> const & Elems        ()         const { return _elems;        }           ///< Access all elements (read-only)
	Array<Element*>       & ElemsWithTag (int Tag);                                           ///< Return the elements with for a given tag
	void                    Bounds       (double & MinX, double & MinY,
	                                      double & MaxX, double & MaxY) const;                ///< Return the limits (bounding box) of the geometry
	void                    Bounds       (double & MinX, double & MinY, double & MinZ,
	                                      double & MaxX, double & MaxY, double & MaxZ) const; ///< Return the limits (bounding box) of the geometry

#ifdef USE_BOOST_PYTHON
// {
	// Set by groups methods
	void PySetNodesElems (Mesh::Generic const & M,          ///< The mesh
	                      BPy::list     const & ElemsAtts); ///< Elements attributes
	void PySetBrys       (Mesh::Generic const & M,          ///< The mesh
	                      BPy::list     const & NodesBrys,  ///< Give [] when there are no nodes boundary conditions
	                      BPy::list     const & EdgesBrys,  ///< Give [] for 3D mesh without edge boundary conditions
	                      BPy::list     const & FacesBrys); ///< Give [] for 2D meshes
	void PyAddLinElems   (BPy::dict     const & Edges,      ///< {(n1,n2):tag1, (n3,n4):tag2, ... num edges} n# => node ID
	                      BPy::list     const & EAtts);     ///< Elements attributes

	// Access methods
	Node    const & PyNod          (size_t i)                     { return (*Nod(i));       }
	Element const & PyEle          (size_t i)                     { return (*Ele(i));       }
	size_t          PyGetNode1     (double X, double Y)           { return GetNode  (X,Y);  }
	size_t          PyGetNode2     (double X, double Y, double Z) { return GetNode  (X,Y,Z);}
	void            PyElemsWithTag (int Tag, BPy::list & Elems);
	void            PyBounds2D     (BPy::list & MinXY,  BPy::list & MaxXY ) const;
	void            PyBounds3D     (BPy::list & MinXYZ, BPy::list & MaxXYZ) const;
// }
#endif // USE_BOOST_PYTHON

private:
	// Data
	int                     _dim;    ///< Space dimension
	double                  _tol;    ///< Tolerance for defining whether X-Y-Z coordinates are in a same plane or not. Also used for comparing coincident nodes
	bool                    _frame;  ///< Only frame (beam/truss) structure ?
	Array<Node*>            _nodes;  ///< FE nodes
	Array<Element*>         _elems;  ///< FE elements
	Array<Element*>         _beams;  ///< Beams
	Array<int>              _btags;  ///< Beam tags
	std::map<int,size_t>    _etidx;  ///< Map: Tag=>Idx, where Idx is the index inside _ewtags, _models, and _props
	Array<Array<Element*> > _ewtags; ///< Elements with tags. Size==_etidx.size()
	Array<Model*>           _models; ///< Models.             Size==_etidx.size()
	Array<Prop_t>           _props;  ///< Properties.         Size==_etidx.size()

}; // class Data

// Forward declaration of a method in embedded.h
void AddReinf (double x1, double y1, double z1, double x2, double y2, double z2, Str_t Prms, bool IsAct, int Tag, Data * G);


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public. */

inline Data::~Data()
{
	for (size_t i=0; i<_nodes .Size(); ++i) if (_nodes [i]!=NULL) delete _nodes [i];
	for (size_t i=0; i<_elems .Size(); ++i) if (_elems [i]!=NULL) delete _elems [i];
	for (size_t i=0; i<_models.Size(); ++i) if (_models[i]!=NULL) delete _models[i];
}

inline void Data::SetNodesElems(Mesh::Generic const * M, EAtts_T const * ElemsAtts)
{
	/* Example:
	
		// Elements attributes
		FEM::EAtts_T eatts;
		eatts.Push (make_tuple(-1, "Quad4", "PStrain", "LinElastic", "E=207 nu=0.3", "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0", "gam=20")); // tag, type, model, prms, inis, props
	*/

	// 3D mesh?
	bool is3d = M->Is3D();

	// Set nodes
	size_t nn = M->NVerts();
	SetNNodes (nn);
	double diff_x = 0.0; // used for verification (try to find which plane the problem is defined in)
	double diff_y = 0.0;
	double diff_z = 0.0;
	for (size_t i=0; i<nn; ++i) // loop over all vertices
	{
		// New node
		SetNode (i, M->VertX(i), M->VertY(i), (is3d ? M->VertZ(i) : 0.0));
		          diff_x += fabs(M->VertX(0)-M->VertX(i));
		          diff_y += fabs(M->VertY(0)-M->VertY(i));
		if (is3d) diff_z += fabs(M->VertZ(0)-M->VertZ(i));
	}

	// Check working plane
	if (is3d)
	{
		if (diff_x<_tol || diff_y<_tol || diff_z<_tol)
			throw new Fatal("FEM::SetNodesElems: For 3D problems, vertices cannot be all in a same plane (diff_x=%f, diff_y=%f, diff_z=%f)",diff_x,diff_y,diff_z);
	}
	else
	{
		if (diff_z>_tol)
			throw new Fatal("FEM::SetNodesElems: For 2D problems, only the X and Y coordinates must be used (diff_z=%f)",diff_z);
	}

	// Number of beams (with duplicates)
	Array< boost::tuple<size_t,size_t,size_t,int,int,int> > beams; // elem_id, eatt_id, local_edge_id, beam_tag, v0, v1
	if (_frame==false)
	{
		for (size_t k=0; k<ElemsAtts->Size(); ++k)
		{
			if (strcmp((*ElemsAtts)[k].get<2>(),"Beam")==0)
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
		SetNBeams (beams.Size());
	}

	// Set elements
	SetNElems (M->NElems() + beams.Size());
	for (size_t i=0; i<M->NElems(); ++i)
	{
		// Set element
		bool found = false;
		for (size_t j=0; j<ElemsAtts->Size(); ++j)
		{
			int tag = (*ElemsAtts)[j].get<0>();
			if (M->ElemTag(i)==tag)
			{
				// Connectivity
				Array<long> conn(M->ElemNVerts(i));
				for (size_t k=0; k<M->ElemNVerts(i); ++k) conn[k] = Nod(M->ElemCon(i,k))->GetID();

				// New finite element
				found = true;
				SetElemAndModel (i, conn, tag,
				                 (*ElemsAtts)[j].get<1>(),  // GeomT
				                 (*ElemsAtts)[j].get<2>(),  // ProbT
				                 (*ElemsAtts)[j].get<3>(),  // MdlName
				                 (*ElemsAtts)[j].get<4>(),  // Prms
				                 (*ElemsAtts)[j].get<5>(),  // Inis
				                 (*ElemsAtts)[j].get<6>(),  // Props
				                 (*ElemsAtts)[j].get<7>()); // IsAct

				// Next element
				break;
			}
		}
		if (found==false) throw new Fatal("SetData: Could not find Tag==%d for Element %d in the ElemsAtts list",M->ElemTag(i),i);
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

		// Connectivity
		Array<long> conn(2);
		conn[0] = Nod(M->EdgeToLef(elem_id, local_edge_id))->GetID();
		conn[1] = Nod(M->EdgeToRig(elem_id, local_edge_id))->GetID();

		// New finite element
		FEM::Element * fe = SetElemAndModel (ie, conn, beam_tag,
		                                     (*ElemsAtts)[eatt_id].get<1>(),  // GeomT
		                                     (*ElemsAtts)[eatt_id].get<2>(),  // ProbT
		                                     (*ElemsAtts)[eatt_id].get<3>(),  // MdlName
		                                     (*ElemsAtts)[eatt_id].get<4>(),  // Prms
		                                     (*ElemsAtts)[eatt_id].get<5>(),  // Inis
		                                     (*ElemsAtts)[eatt_id].get<6>(),  // Props
		                                     (*ElemsAtts)[eatt_id].get<7>()); // IsAct

		// Next element
		ie++;

		// Set beam
		SetBeam (i, fe, beam_tag);
	}
}

inline void Data::SetBrys(Mesh::Generic const * M, NBrys_T const * NodesBrys, EBrys_T const * EdgesBrys, FBrys_T const * FacesBrys)
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
							Ele(i)->FaceBry ((*FacesBrys)[k].get<1>(), (*FacesBrys)[k].get<2>(), j);
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
							Ele(i)->EdgeBry ((*EdgesBrys)[k].get<1>(), (*EdgesBrys)[k].get<2>(), j);
						}
					}
				}
			}
			for (size_t b=0; b<NBeams(); ++b)
			{
				if (BTag(b)==(*EdgesBrys)[k].get<0>())
					Beam(b)->EdgeBry ((*EdgesBrys)[k].get<1>(), (*EdgesBrys)[k].get<2>(), 0);
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
				if (d<_tol) Nod(i)->Bry ((*NodesBrys)[j].get<3>(), (*NodesBrys)[j].get<4>());
			}
		}
	}
}

inline void Data::SetNNodes(size_t NNodes)
{
	for (size_t i=0; i<_nodes.Size(); ++i) if (_nodes[i]!=NULL) delete _nodes[i];
	_nodes.Resize(NNodes);
	_nodes.SetValues(NULL);
}

inline void Data::SetNElems(size_t NElems)
{
	for (size_t i=0; i<_elems.Size(); ++i) if (_elems[i]!=NULL) delete _elems[i];
	_elems.Resize(NElems);
	_elems.SetValues(NULL);
}

inline Node * Data::SetNode(size_t i, double X, double Y, double Z, int Tag)
{
   if (_nodes[i]==NULL) _nodes[i] = new Node;
	_nodes[i]->Initialize (i, X,Y,Z, Tag);
	return _nodes[i];
}

inline Element * Data::SetElemAndModel(size_t i, Array<long> const & Conn, int Tag, Str_t GeomT, Str_t ProbT, Str_t MdlName, Str_t Prms, Str_t Inis, Str_t Props, bool IsAct)
{
	// New element
	if (_elems[i]==NULL) _elems[i] = new Element;

	// Initialize constants of element and get geometry index (gi)
	int gi = _elems[i]->InitCtes (_dim, GeomT, ProbT);

	// Set array with elements with tags, allocate models, and set _props
	if (_etidx.count(Tag)==0) // tag not set yet
	{
		// Set elements with tags
		_etidx[Tag] = _ewtags.Size();
		Array<Element*> tmp;
		_ewtags.Push (tmp);
		
		// Allocate model (MUST be after _elems[i]->InitCtes)
		_models.Push (AllocModel(MdlName));
		_models[_models.Size()-1]->SetGeomIdx (gi);        // Set Model geometry index
		_models[_models.Size()-1]->Initialize (Tag, Prms); // Parse parameters

		// Parse properties
		Prop_t prp;
		_props.Push (prp);
		LineParser lp(Props);
		lp.ReadVariables (_elems[i]->NProps(), _elems[i]->Props(), _props[_props.Size()-1], "properties", "Element tag", Tag);
	}
	_ewtags[_etidx[Tag]].Push (_elems[i]);

	// Initialize
	Array<Node*> conn(Conn.Size());
	for (size_t j=0; j<Conn.Size(); ++j) conn[j] = Nod(Conn[j]);
	_elems[i]->Initialize (/*ID*/i, Tag, conn, _models[_etidx[Tag]], Inis, &_props[_etidx[Tag]], IsAct);

	return _elems[i];
}

inline Element * Data::PushElem(Array<long> const & Conn, int Tag, Str_t GeomT, Str_t ProbT, Str_t MdlName, Str_t Prms, Str_t Inis, Str_t Props, bool IsAct)
{
	// Add new element
	_elems.Push (new Element);

	// Set element, model, and properties
	return SetElemAndModel (_elems.Size()-1, Conn, Tag, GeomT, ProbT, MdlName, Prms, Inis, Props, IsAct);
}

inline void Data::ClearDisp()
{
	for (size_t i=0; i<_elems.Size(); ++i) _elems[i]->ClearDisp();
}

inline void Data::Activate(int ElemTag)
{
	Array<Element*> & elems = ElemsWithTag (ElemTag);
	for (size_t i=0; i<elems.Size(); ++i) elems[i]->SetActive (true);
}

inline void Data::Deactivate(int ElemTag)
{
	Array<Element*> & elems = ElemsWithTag (ElemTag);
	for (size_t i=0; i<elems.Size(); ++i) elems[i]->SetActive (false);
}

inline bool Data::Check()
{
	// Check arrays
	if (NNodes()==0 || NElems()==0) return false;

	// Check nodes
	for (size_t i=0; i<NNodes(); ++i) if (_nodes[i]==NULL) return false;

	// Check elements
	for (size_t i=0; i<NElems(); ++i) if (_elems[i]==NULL) return false;

	return true; // OK
}

inline size_t Data::GetNode(double X, double Y, double Z)
{
	for (size_t i=0; i<_nodes.Size(); i++)
	{
		if( _nodes[i]->X()==X && _nodes[i]->Y()==Y && _nodes[i]->Z()==Z) return i;
	//	std::cout << _nodes[i]->X() << " " << _nodes[i]->Y() << " " << _nodes[i]->Z();
	}
	throw new Fatal("Data::GetNode: Node not found (%g, %g, %g)", X, Y, Z);
}

inline Array<Element*> & Data::ElemsWithTag(int Tag)
{
	if (_etidx.count(Tag)==0) throw new Fatal("Data::ElemsWithTag: This Tag==%d was not set for any Element",Tag);
	return _ewtags[_etidx[Tag]];
}

inline void Data::Bounds(double & MinX, double & MinY, double & MaxX, double & MaxY) const
{
	MinX = (NNodes()>0 ? Nod(0)->X() : 0.0);
	MinY = (NNodes()>0 ? Nod(0)->Y() : 0.0);
	MaxX = (NNodes()>0 ? Nod(0)->X() : 0.0);
	MaxY = (NNodes()>0 ? Nod(0)->Y() : 0.0);
	for (size_t i=0; i<NNodes(); ++i)
	{
		if (Nod(i)->X() < MinX) MinX = Nod(i)->X();
		if (Nod(i)->Y() < MinY) MinY = Nod(i)->Y();
		if (Nod(i)->X() > MaxX) MaxX = Nod(i)->X();
		if (Nod(i)->Y() > MaxY) MaxY = Nod(i)->Y();
	}
}

inline void Data::Bounds(double & MinX, double & MinY, double & MinZ, double & MaxX, double & MaxY, double & MaxZ) const
{
	MinX = (NNodes()>0 ? Nod(0)->X() : 0.0);
	MinY = (NNodes()>0 ? Nod(0)->Y() : 0.0);
	MinZ = (NNodes()>0 ? Nod(0)->Z() : 0.0);
	MaxX = (NNodes()>0 ? Nod(0)->X() : 0.0);
	MaxY = (NNodes()>0 ? Nod(0)->Y() : 0.0);
	MaxZ = (NNodes()>0 ? Nod(0)->Z() : 0.0);
	for (size_t i=0; i<NNodes(); ++i)
	{
		if (Nod(i)->X() < MinX) MinX = Nod(i)->X();
		if (Nod(i)->Y() < MinY) MinY = Nod(i)->Y();
		if (Nod(i)->Z() < MinZ) MinZ = Nod(i)->Z();
		if (Nod(i)->X() > MaxX) MaxX = Nod(i)->X();
		if (Nod(i)->Y() > MaxY) MaxY = Nod(i)->Y();
		if (Nod(i)->Z() > MaxZ) MaxZ = Nod(i)->Z();
	}
}


/** Outputs a data structure. */
std::ostream & operator<< (std::ostream & os, FEM::Data const & G)
{
	for (size_t i=0; i<G.NElems(); ++i)
		if (G.Ele(i)!=NULL) os << (*G.Ele(i));
	return os;
}


#ifdef USE_BOOST_PYTHON
// {

inline void Data::PySetNodesElems(Mesh::Generic const & M, BPy::list const & ElemsAtts)
{
	/* Example:
	 *           
	 *           # Elements attributes
	 *           eatts = [[-1, 'Quad4PStrain', 'LinElastic', 'E=%f nu=%f'%(E,nu), 'Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0', 'gam=20', True]] # tag, type, model, prms, inis, props, active?
	 */

	// Extract list with elements attributes
	FEM::EAtts_T eatts;
	int eatts_size = len(ElemsAtts);
	if (eatts_size>0) eatts.Resize(eatts_size);
	for (int i=0; i<eatts_size; ++i)
	{
		if (len(ElemsAtts[i])==7)
		{
			BPy::list lst = BPy::extract<BPy::list>(ElemsAtts[i])();
			eatts[i] = boost::make_tuple(BPy::extract<int        >(lst[0])(),
			                             BPy::extract<Str_t>(lst[1])(),
			                             BPy::extract<Str_t>(lst[2])(),
			                             BPy::extract<Str_t>(lst[3])(),
			                             BPy::extract<Str_t>(lst[4])(),
			                             BPy::extract<Str_t>(lst[5])(), 
			                             BPy::extract<bool>       (lst[6])());
		}
		else throw new Fatal("PySetNodesElems: Each sublist in ElemsAtts must have 7 items: tag, type, model, prms, inis, props, active?\n\tExample: ElemsAtts = [[-1, 'Quad4PStrain', 'LinElastic', 'E=207.0 nu=0.3', 'Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0', 'gam=20', True]]");
	}

	// Set geometry
	SetNodesElems (&M, &eatts);
}

inline void Data::PySetBrys(Mesh::Generic const & M, BPy::list const & NodesBrys, BPy::list const & EdgesBrys, BPy::list const & FacesBrys)
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
			                                BPy::extract<Str_t>(lst[3])(),
			                                BPy::extract<double     >(lst[4])());
		}
		else throw new Fatal("PySetData: Each sublist in NodesBrys must have 5 items: x,y,z, key, val\n\tExample: NodesBrys = [[1.0, 0.0, 0.0, 'ux', 0.0]]");
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
			                                BPy::extract<Str_t>(lst[1])(),
			                                BPy::extract<double     >(lst[2])());
		}
		else throw new Fatal("PySetData: Each sublist in EdgesBrys must have 3 items: tag, key, val\n\tExample: EdgesBrys = [[-10, 'uy', 0.0], [-20, 'fy', -1]]");
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
			                                BPy::extract<Str_t>(lst[1])(),
			                                BPy::extract<double     >(lst[2])());
		}
		else throw new Fatal("PySetData: Each sublist in FacesBrys must have 3 items: tag, key, val\n\tExample: FacesBrys = [[-10, 'uy', 0.0], [-20, 'fy', -1]]");
	}

	// Set geometry
	SetBrys (&M, nbrys, ebrys, fbrys);

	// Clean up
	if (nbrys!=NULL) delete nbrys;
	if (ebrys!=NULL) delete ebrys;
	if (fbrys!=NULL) delete fbrys;
}

inline void Data::PyAddLinElems(BPy::dict const & Edges, BPy::list const & EAtts)
{
	/* Example:
	 *           
	 *           # Elements attributes
	 *           eatts = [[-1, '', 'Spring', 'LinSpring', 'ks=%g', 'ZERO', 'gam=20', True]] # tag, type, model, prms, inis, props, active?
	 */

	// Map element tag to index in EAtts list
	int neatts = BPy::len(EAtts);
	if (neatts<1) throw new Fatal("Data::PyAddLinElems: EAtts (element attributes) must contain at least one element");
	std::map<int,int> tag2idx; 
	for (int i=0; i<neatts; ++i)
	{
		BPy::list const & lst = BPy::extract<BPy::list>(EAtts[i])();
		tag2idx[BPy::extract<int>(lst[0])()] = i;
		if (BPy::len(EAtts[i])!=8) throw new Fatal("Data::PyAddLinElems: Each sublist in EAtts must have 8 items: tag, '', type, model, prms, inis, props, active?\n\tExample: eatts = [[-1, 'Spring', '', 'ks=1e+12', 'ZERO', 'gam=20', True]]\n\tlen(EAtts[i])==%d is invalid.",BPy::len(EAtts[i]));
	}

	// Read edges
	BPy::object const & e_keys = BPy::extract<BPy::dict>(Edges)().iterkeys();
	BPy::object const & e_vals = BPy::extract<BPy::dict>(Edges)().itervalues();
	for (int i=0; i<BPy::len(Edges); ++i)
	{
		// Extract linear element data
		Array<long> conn(2); // connectivity
		BPy::tuple const & edge    = BPy::extract<BPy::tuple> (e_keys.attr("next")())();
		int                tag     = BPy::extract<int>        (e_vals.attr("next")())();
		                   conn[0] = BPy::extract<long>       (edge[0])();
		                   conn[1] = BPy::extract<long>       (edge[1])();
		
		// Find element attributes
		std::map<int,int>::const_iterator iter = tag2idx.find(tag);
		if (iter==tag2idx.end()) throw new Fatal("Data::PyAddLinElems: Could not find tag < %d > in the list of Element Attributes", tag);
		int idx_eatt = iter->second;

		// Add linear element to FE geometry
		PushElem (conn, tag,
		          BPy::extract<Str_t>(EAtts[idx_eatt][1])(),  // GeomT
		          BPy::extract<Str_t>(EAtts[idx_eatt][2])(),  // ProbT
		          BPy::extract<Str_t>(EAtts[idx_eatt][3])(),  // MdlName
		          BPy::extract<Str_t>(EAtts[idx_eatt][4])(),  // Prms
		          BPy::extract<Str_t>(EAtts[idx_eatt][5])(),  // Inis
		          BPy::extract<Str_t>(EAtts[idx_eatt][6])(),  // Props
		          BPy::extract<bool> (EAtts[idx_eatt][7])()); // IsAct
	}
}

inline void Data::PyElemsWithTag(int Tag, BPy::list & Elems)
{
	Array<FEM::Element*> & elems = ElemsWithTag (Tag);
	for (size_t i=0; i<elems.Size(); ++i)
		Elems.append ((*elems[i]));
}

inline void Data::PyBounds2D(BPy::list & MinXY, BPy::list & MaxXY) const
{
	double  minx,miny, maxx,maxy;
	Bounds (minx,miny, maxx,maxy);
	MinXY.append(minx);  MaxXY.append(maxx);
	MinXY.append(miny);  MaxXY.append(maxy);
}

inline void Data::PyBounds3D(BPy::list & MinXYZ, BPy::list & MaxXYZ) const
{
	double  minx,miny,minz, maxx,maxy,maxz;
	Bounds (minx,miny,minz, maxx,maxy,maxz);
	MinXYZ.append(minx);  MaxXYZ.append(maxx);
	MinXYZ.append(miny);  MaxXYZ.append(maxy);
	MinXYZ.append(minz);  MaxXYZ.append(maxz);
}

// }
#endif // USE_BOOST_PYTHON

}; // namespace FEM

#endif // MECHSYS_FEM_DATA_H
