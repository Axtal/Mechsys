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

namespace FEM
{

// Tuples used when initializing nodes and elements by groups
typedef Array< boost::tuple<double,double,double, char const *,double> > NBrys_T; // Node: x,y,z, key, val
typedef Array< boost::tuple<                 int, char const *,double> > EBrys_T; // Edge:   tag, key, val
typedef Array< boost::tuple<                 int, char const *,double> > FBrys_T; // Face:   tag, key, val
typedef Array< boost::tuple<int, char const*, char const*, char const*,
                            char const*, char const*, bool> > EAtts_T; // Elem: tag, type, model, prms, inis, props, active

/* Dataetry */
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
	void      SetNNodes (size_t NNodes);                                         ///< Set the number of nodes
	void      SetNElems (size_t NElems);                                         ///< Set the number of elements
	Node    * SetNode   (size_t i, double X, double Y, double Z=0.0, int Tag=0); ///< Set a node
	Element * SetElem   (size_t i, char const * Type, bool IsActive, int Tag);   ///< Set an element

	// Add methods
	size_t PushNode (Node * N)    { _nodes.Push(N); return _nodes.Size()-1; } ///< Push back a new node
	size_t PushElem (Element * E) { _elems.Push(E); return _elems.Size()-1; } ///< Push back a new element
	size_t PushNode (double X, double Y, double Z=0.0, int Tag=0);            ///< Push back a new node
	size_t PushElem (int Tag,                                                 ///< The tag of new element
	                 char const * Type,                                       ///< Element type. Ex: Tri3PStrain
	                 char const * Model,                                      ///< Constitutive model. Ex: LinElastic
	                 char const * Prms,                                       ///< Model parameters. Ex.: E=200 nu=0.2
	                 char const * Inis,                                       ///< Initial state. Ex: Sx=0.0
	                 char const * Props,                                      ///< Element properties. Ex: gam=20
	                 bool IsActive,                                           ///< Active/Inactive?
	                 Array<int> const & Connectivity);                        ///< Push back a new element

	// Specific methods
	void ApplyBodyForces    () { for (size_t i=0; i<_elems.Size(); ++i) _elems[i]->ApplyBodyForces(); } ///< Apply body forces (equilibrium/coupled problems)
	void ClearDisplacements ();                                                                         ///< Clear displacements (equilibrium/coupled problems)
	void Activate           (int ElemTag);                                                              ///< Activate all elements with Tag
	void Deactivate         (int ElemTag);                                                              ///< Activate all elements with Tag

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
	void PyAddReinfs     (BPy::dict const & Edges,          ///< {(x0,y0, x1,y1):tag1, ... num edges} or {(x0,y0,z0, x1,y1,z1):tag1, ... num edges}
	                      BPy::list const & EAtts);         ///< Elements attributes
	void PyAddLinElems   (BPy::dict const & Edges,          ///< {(n1,n2):tag1, (n3,n4):tag2, ... num edges} n# => node ID
	                      BPy::list const & EAtts);         ///< Elements attributes

	// Set methods
	Node   & PySetNode2D (size_t i, double X, double Y)                       { return (*SetNode(i,X,Y));   }
	Node   & PySetNode3D (size_t i, double X, double Y, double Z)             { return (*SetNode(i,X,Y,Z)); }
	PyElem   PySetElem   (size_t i, BPy::str const & Type, bool Act, int Tag) { return PyElem(SetElem(i,BPy::extract<char const *>(Type)(),Act,Tag)); }

	// Add methods
	size_t PyPushNode1 (double X, double Y)                    { return PushNode (X,Y);       }
	size_t PyPushNode2 (double X, double Y, double Z)          { return PushNode (X,Y,Z);     }
	size_t PyPushNode3 (double X, double Y, double Z, int Tag) { return PushNode (X,Y,Z,Tag); }
	size_t PyPushElem  (int Tag, BPy::str const & Type, BPy::str const & Model, BPy::str const & Prms, BPy::str const & Inis, BPy::str const & Props, bool IsActive, BPy::list const & Connectivity);

	// Access methods
	Node const & PyNod          (size_t i)                     { return (*Nod(i));       }
	PyElem       PyEle          (size_t i)                     { return PyElem(Ele(i));  }
	size_t       PyGetNode1     (double X, double Y)           { return GetNode  (X,Y);  }
	size_t       PyGetNode2     (double X, double Y, double Z) { return GetNode  (X,Y,Z);}
	void         PyElemsWithTag (int Tag, BPy::list & Elems);
	void         PyBounds2D     (BPy::list & MinXY,  BPy::list & MaxXY ) const;
	void         PyBounds3D     (BPy::list & MinXYZ, BPy::list & MaxXYZ) const;
// }
#endif // USE_BOOST_PYTHON

private:
	// Data
	int                     _dim;             ///< Space dimension
	double                  _tol;             ///< Tolerance for defining whether X-Y-Z coordinates are in a same plane or not. Also used for comparing coincident nodes
	bool                    _frame;           ///< Only frame (beam/truss) structure ?
	Array<Node*>            _nodes;           ///< FE nodes
	Array<Element*>         _elems;           ///< FE elements
	Array<Element*>         _beams;           ///< Beams
	Array<int>              _btags;           ///< Beam tags
	std::map<int,size_t>    _elem_tag_idx;    ///< Map Tag => Idx, where Idx is the index inside _elems_with_tags
	Array<Array<Element*> > _elems_with_tags; ///< Element with tags

}; // class Data

// Forward declaration of a method in embedded.h
void AddReinf (double x1, double y1, double z1, double x2, double y2, double z2, char const * Prms, bool IsActive, int Tag, Data * G);


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public. */

inline Data::~Data()
{
	for (size_t i=0; i<_nodes.Size(); ++i) if (_nodes[i]!=NULL) delete _nodes[i];
	for (size_t i=0; i<_elems.Size(); ++i) if (_elems[i]!=NULL) delete _elems[i];
}


// Set by groups methods

inline void Data::SetNodesElems(Mesh::Generic const * M, EAtts_T const * ElemsAtts)
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
	this->SetNNodes (nn);
	double diff_x = 0.0; // used for verification (try to find which plane the problem is defined in)
	double diff_y = 0.0;
	double diff_z = 0.0;
	for (size_t i=0; i<nn; ++i) // loop over all vertices
	{
		// New node
		this->SetNode (i, M->VertX(i), M->VertY(i), (is3d ? M->VertZ(i) : 0.0));
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
		this->SetNBeams (beams.Size());
	}

	// Set elements
	this->SetNElems (M->NElems() + beams.Size());
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
				FEM::Element * fe = this->SetElem (i, (*ElemsAtts)[j].get<1>(), (*ElemsAtts)[j].get<6>(), M->ElemTag(i));

				// Set connectivity
				if (M->ElemNVerts(i)!=fe->NNodes()) throw new Fatal("functions.h::SetNodesElems:: The number of vertices (%d) of mesh object must be compatible with the number of nodes (%d) of the element (%s)",M->ElemNVerts(i),fe->NNodes(),fe->Name());
				for (size_t k=0; k<M->ElemNVerts(i); ++k)
					fe->Connect (k, this->Nod(M->ElemCon(i,k)));

				// Set parameters and initial values
				fe->SetModel ((*ElemsAtts)[j].get<2>(), (*ElemsAtts)[j].get<3>(), (*ElemsAtts)[j].get<4>());

				// Set properties
				fe->SetProps ((*ElemsAtts)[j].get<5>());
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

		// New finite element
		FEM::Element * fe = this->SetElem (ie, (*ElemsAtts)[eatt_id].get<1>(), (*ElemsAtts)[eatt_id].get<6>(), beam_tag);

		// Set connectivity
		fe->Connect (0, this->Nod(M->EdgeToLef(elem_id, local_edge_id)));
		fe->Connect (1, this->Nod(M->EdgeToRig(elem_id, local_edge_id)));

		// Set parameters and initial values
		fe->SetModel ((*ElemsAtts)[eatt_id].get<2>(), (*ElemsAtts)[eatt_id].get<3>(), (*ElemsAtts)[eatt_id].get<4>());

		// Set properties
		fe->SetProps ((*ElemsAtts)[eatt_id].get<5>());
		ie++;

		// Set beam
		this->SetBeam (i, fe, beam_tag);
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
							this->Ele(i)->FaceBry ((*FacesBrys)[k].get<1>(), (*FacesBrys)[k].get<2>(), j);
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
							this->Ele(i)->EdgeBry ((*EdgesBrys)[k].get<1>(), (*EdgesBrys)[k].get<2>(), j);
						}
					}
				}
			}
			for (size_t b=0; b<this->NBeams(); ++b)
			{
				if (this->BTag(b)==(*EdgesBrys)[k].get<0>())
					this->Beam(b)->EdgeBry ((*EdgesBrys)[k].get<1>(), (*EdgesBrys)[k].get<2>(), 0);
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
				if (d<_tol) this->Nod(i)->Bry ((*NodesBrys)[j].get<3>(), (*NodesBrys)[j].get<4>());
			}
		}
	}
}


// Set methods

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
	_nodes[i]->Initialize (i,X,Y,Z, Tag);
	return _nodes[i];
}

inline Element * Data::SetElem(size_t i, char const * Type, bool IsActive, int Tag)
{
	if (_elems[i]==NULL) _elems[i] = AllocElement(Type);
	_elems[i]->Initialize (/*ID*/i, IsActive, _dim, Tag);
	if (Tag!=0)
	{
		if (_elem_tag_idx.count(Tag)==0) // tag not set
		{
			_elem_tag_idx[Tag] = _elems_with_tags.Size();
			Array<Element*> tmp;
			_elems_with_tags.Push (tmp);
		}
		_elems_with_tags[_elem_tag_idx[Tag]].Push (_elems[i]);
	}
	return _elems[i];
}


// Add methods

inline size_t Data::PushNode(double X, double Y, double Z, int Tag)  
{
	// Allocate new node
	Node * new_node;
	new_node = new Node;
	_nodes.Push(new_node);
	size_t ID = _nodes.Size()-1;
	new_node->Initialize(ID, X, Y, Z, Tag);
	return ID;
}

inline size_t Data::PushElem(int Tag, char const * Type, char const * Model, char const * Prms, char const * Inis, char const * Props, bool IsActive, Array<int> const & Connectivity )
{
	// Alocate new element
	Element * new_elem;
	new_elem = AllocElement(Type);
	_elems.Push(new_elem);
	size_t ID = _elems.Size()-1;
	new_elem->Initialize(_elems.Size()-1, IsActive, _dim, Tag);

	// Connector (EmbSpring) connectivities
	if (new_elem->NNodes()!=Connectivity.Size()) throw new Fatal("Data::PushElem: The number of nodes in Connectivity does not match. (Element ID %d)", ID);
	for (size_t i=0; i<new_elem->NNodes(); i++)
		new_elem->Connect(i, _nodes[Connectivity[i]]);

	// Set the model 
	new_elem->SetModel(Model, Prms, Inis);
	new_elem->SetProps(Props);
	return ID;
}


// Specific methods

inline void Data::ClearDisplacements()
{
	for (size_t i=0; i<_elems.Size(); ++i) _elems[i]->ClearDispAndStrains();
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


// Access methods

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
		std::cout << _nodes[i]->X() << " " << _nodes[i]->Y() << " " << _nodes[i]->Z();
	}
	throw new Fatal("Data::GetNode: Node not found (%g, %g, %g)", X, Y, Z);
}

inline Array<Element*> & Data::ElemsWithTag(int Tag)
{
	if (_elem_tag_idx.count(Tag)==0) throw new Fatal("Data::ElemsWithTag: This Tag==%d was not set for any Element",Tag);
	return _elems_with_tags[_elem_tag_idx[Tag]];
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
			                             BPy::extract<char const*>(lst[1])(),
			                             BPy::extract<char const*>(lst[2])(),
			                             BPy::extract<char const*>(lst[3])(),
			                             BPy::extract<char const*>(lst[4])(),
			                             BPy::extract<char const*>(lst[5])(), 
			                             BPy::extract<bool>       (lst[6])());
		}
		else throw new Fatal("PySetNodesElems: Each sublist in ElemsAtts must have 7 items: tag, type, model, prms, inis, props, active?\n\tExample: ElemsAtts = [[-1, 'Quad4PStrain', 'LinElastic', 'E=207.0 nu=0.3', 'Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0', 'gam=20', True]]");
	}

	// Set geometry
	this->SetNodesElems (&M, &eatts);
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
			                                BPy::extract<char const*>(lst[3])(),
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
			                                BPy::extract<char const*>(lst[1])(),
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
			                                BPy::extract<char const*>(lst[1])(),
			                                BPy::extract<double     >(lst[2])());
		}
		else throw new Fatal("PySetData: Each sublist in FacesBrys must have 3 items: tag, key, val\n\tExample: FacesBrys = [[-10, 'uy', 0.0], [-20, 'fy', -1]]");
	}

	// Set geometry
	this->SetBrys (&M, nbrys, ebrys, fbrys);

	// Clean up
	if (nbrys!=NULL) delete nbrys;
	if (ebrys!=NULL) delete ebrys;
	if (fbrys!=NULL) delete fbrys;
}

inline void Data::PyAddReinfs(BPy::dict const & Edges, BPy::list const & EAtts)
{
	/* Example:
	 *           
	 *           # Elements attributes
	 *           eatts = [[-1, 'Reinforcement', '', 'E=%g Ar=%g At=%g ks=%g c=%g phi=%g', 'ZERO', 'gam=20', True]] # tag, type, model, prms, inis, props, active?
	 */

	// Map element tag to index in EAtts list
	int neatts = BPy::len(EAtts);
	if (neatts<1) throw new Fatal("Data::PyAddReinfs: EAtts (element attributes) must contain at least one element");
	std::map<int,int> tag2idx; 
	for (int i=0; i<neatts; ++i)
	{
		BPy::list const & lst = BPy::extract<BPy::list>(EAtts[i])();
		tag2idx[BPy::extract<int>(lst[0])()] = i;
		if (BPy::len(EAtts[i])!=7) throw new Fatal("Data::PyAddReinfs: Each sublist in EAtts must have 7 items: tag, type, model, prms, inis, props, active?\n\tExample: eatts = [[-1, 'Reinforcement', '', 'E=200 Ar=1 At=1 ks=1e+12 c=0 phi=20', 'ZERO', 'gam=20', True]]\n\tlen(EAtts[i])==%d is invalid.",BPy::len(EAtts[i]));
	}

	// Read edges
	double x0=0.0; double y0=0.0; double z0=0.0;
	double x1=0.0; double y1=0.0; double z1=0.0;
	BPy::object const & e_keys = BPy::extract<BPy::dict>(Edges)().iterkeys();
	BPy::object const & e_vals = BPy::extract<BPy::dict>(Edges)().itervalues();
	for (int i=0; i<BPy::len(Edges); ++i)
	{
		// Extract reinforcement data
		BPy::tuple const & xyz = BPy::extract<BPy::tuple> (e_keys.attr("next")())();
		int                tag = BPy::extract<int>        (e_vals.attr("next")())();

		// Check number of coordinates
		size_t nxyz = BPy::len(xyz);
		bool   is3d = (nxyz==4 ? false : true);
		if ((nxyz!=4) || (nxyz!=6)) throw new Fatal("Data::PyAddReinfs: Each edge representing a reinforcement must have either 4(2D) or 6(3D) coordinates corresponding to the start and end points (nxyz==%d is invalid).\n\tEx. {(x0,y0,z0, x1,y1,z1):-100}",nxyz);

		// Extract coordinates
		if (is3d)
		{
			x0 = BPy::extract<double>(xyz[0])();
			y0 = BPy::extract<double>(xyz[1])();
			z0 = BPy::extract<double>(xyz[2])();
			x1 = BPy::extract<double>(xyz[3])();
			y1 = BPy::extract<double>(xyz[4])();
			z1 = BPy::extract<double>(xyz[5])();
		}
		else
		{
			x0 = BPy::extract<double>(xyz[0])();
			y0 = BPy::extract<double>(xyz[1])();
			x1 = BPy::extract<double>(xyz[2])();
			y1 = BPy::extract<double>(xyz[3])();
		}
		std::cout << "Edge: tag="<<tag << ", x0="<<x0 << ", y0="<<y0 << ", z0="<<z0 << ", x1="<<x1 << ", y1="<<y1 << ", z1="<<z1 << std::endl;
		
		// Find element attributes
		std::map<int,int>::const_iterator iter = tag2idx.find(tag);
		if (iter==tag2idx.end()) throw new Fatal("Data::PyAddReinfs: Could not find tag < %d > in the list of Element Attributes", tag);
		int idx_eatt = iter->second;

		// Add reinforcement to FE geometry
		char const * prms = BPy::extract<char const *>(EAtts[idx_eatt][3])();
		bool         act  = BPy::extract<bool>        (EAtts[idx_eatt][6])();
		FEM::AddReinf (x0,y0,z0, x1,y1,z1, prms, act, tag, this);
	}
}

inline void Data::PyAddLinElems(BPy::dict const & Edges, BPy::list const & EAtts)
{
	/* Example:
	 *           
	 *           # Elements attributes
	 *           eatts = [[-1, 'Spring', '', 'ks=%g', 'ZERO', 'gam=20', True]] # tag, type, model, prms, inis, props, active?
	 */

	// Map element tag to index in EAtts list
	int neatts = BPy::len(EAtts);
	if (neatts<1) throw new Fatal("Data::PyAddLinElems: EAtts (element attributes) must contain at least one element");
	std::map<int,int> tag2idx; 
	for (int i=0; i<neatts; ++i)
	{
		BPy::list const & lst = BPy::extract<BPy::list>(EAtts[i])();
		tag2idx[BPy::extract<int>(lst[0])()] = i;
		if (BPy::len(EAtts[i])!=7) throw new Fatal("Data::PyAddLinElems: Each sublist in EAtts must have 7 items: tag, type, model, prms, inis, props, active?\n\tExample: eatts = [[-1, 'Spring', '', 'ks=1e+12', 'ZERO', 'gam=20', True]]\n\tlen(EAtts[i])==%d is invalid.",BPy::len(EAtts[i]));
	}

	// Read edges
	BPy::object const & e_keys = BPy::extract<BPy::dict>(Edges)().iterkeys();
	BPy::object const & e_vals = BPy::extract<BPy::dict>(Edges)().itervalues();
	for (int i=0; i<BPy::len(Edges); ++i)
	{
		// Extract linear element data
		Array<int> conn(2); // connectivity
		BPy::tuple const & edge    = BPy::extract<BPy::tuple> (e_keys.attr("next")())();
		int                tag     = BPy::extract<int>        (e_vals.attr("next")())();
		                   conn[0] = BPy::extract<int>        (edge[0])();
		                   conn[1] = BPy::extract<int>        (edge[1])();
		
		// Find element attributes
		std::map<int,int>::const_iterator iter = tag2idx.find(tag);
		if (iter==tag2idx.end()) throw new Fatal("Data::PyAddLinElems: Could not find tag < %d > in the list of Element Attributes", tag);
		int idx_eatt = iter->second;

		// Add linear element to FE geometry
		char const * type  = BPy::extract<char const *>(EAtts[idx_eatt][1])();
		char const * mdl   = BPy::extract<char const *>(EAtts[idx_eatt][2])();
		char const * prms  = BPy::extract<char const *>(EAtts[idx_eatt][3])();
		char const * inis  = BPy::extract<char const *>(EAtts[idx_eatt][4])();
		char const * props = BPy::extract<char const *>(EAtts[idx_eatt][5])();
		bool         act   = BPy::extract<bool>        (EAtts[idx_eatt][6])();
		PushElem (tag, type, mdl, prms, inis, props, act, conn);
	}
}

inline size_t Data::PyPushElem(int Tag, BPy::str const & Type, BPy::str const & Model, BPy::str const & Prms, BPy::str const & Inis, BPy::str const & Props, bool IsActive, BPy::list const & Connectivity)
{
	size_t nnodes = BPy::len(Connectivity);
	if (nnodes<2) throw new Fatal("Data::PyPushElem: Number of nodes in the Connectivity list must be greater than 1");
	Array<int> conn(nnodes);
	for (size_t i=0; i<nnodes; ++i) conn[i] = BPy::extract<int>(Connectivity[i])();
	return PushElem(Tag,
	                BPy::extract<char const *>(Type)(),
	                BPy::extract<char const *>(Model)(),
	                BPy::extract<char const *>(Prms)(),
	                BPy::extract<char const *>(Inis)(),
	                BPy::extract<char const *>(Props)(),
	                IsActive, conn);
}

inline void Data::PyElemsWithTag(int Tag, BPy::list & Elems)
{
	Array<FEM::Element*> & elems = ElemsWithTag (Tag);
	for (size_t i=0; i<elems.Size(); ++i)
		Elems.append (PyElem(elems[i]));
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
