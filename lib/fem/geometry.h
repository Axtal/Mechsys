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

#ifndef MECHSYS_FEM_GEOMETRY_H
#define MECHSYS_FEM_GEOMETRY_H

// STL
#include <iostream>
#include <cstring>
#include <map>

// Boost::Python
#ifdef USE_BOOST_PYTHON
  //#include <boost/python.hpp> // this includes everything
  namespace BPy = boost::python;
#endif

// MechSys
#include "fem/node.h"
#include "fem/element.h"
#include "util/array.h"

namespace FEM
{

/* Geometry */
class Geom
{
public:
	/* Constructor */
	Geom (int nDim) : _dim(nDim) {}

	/* Destructor */
	~Geom ();

	// Set methods
	void      SetNNodes (size_t NNodes);                                              ///< Set the number of nodes
	void      SetNElems (size_t NElems);                                              ///< Set the number of elements
	Node    * SetNode   (size_t i, double X, double Y, double Z=0.0, int Tag=0);      ///< Set a node
	Element * SetElem   (size_t i, char const * Type, bool IsActive, int Tag);        ///< Set an element

	// Specific methods
	void ApplyBodyForces    () { for (size_t i=0; i<_elems.Size(); ++i) _elems[i]->ApplyBodyForces(); } ///< Apply body forces (equilibrium/coupled problems)
	void ClearDisplacements ();                                                                         ///< Clear displacements (equilibrium/coupled problems)
	void Activate           (int ElemTag);                                                              ///< Activate all elements with Tag
	void Deactivate         (int ElemTag);                                                              ///< Activate all elements with Tag

	// Beam
	void      SetNBeams (size_t NBeams) { _beams.Resize(NBeams); _beams.SetValues(NULL); _btags.Resize(NBeams); }
	void      SetBeam   (size_t iBeam, Element * Beam, int Tag) { _beams[iBeam]=Beam; _btags[iBeam]=Tag; }
	size_t    NBeams    () const       { return _beams.Size(); }
	Element * Beam      (size_t iBeam) { return _beams[iBeam]; }
	int       BTag      (size_t iBeam) { return _btags[iBeam]; }

	// Access methods
	bool                    Check        ();                                        ///< Check if Nodes and Elements were allocated properly. Should be called before accessing Nodes and Elements, since these may not had been allocated yet (and then causing Segfaults).
	size_t                  NNodes       ()         const { return _nodes.Size(); } ///< Return the number of nodes
	size_t                  NElems       ()         const { return _elems.Size(); } ///< Return the number of elements
	size_t                  NDim         ()         const { return _dim;          } ///< Return the dimension
	Node                  * Nod          (size_t i)       { return _nodes[i];     } ///< Access (read/write) a node
	Element               * Ele          (size_t i)       { return _elems[i];     } ///< Access (read/write) an element
	Node            const * Nod          (size_t i) const { return _nodes[i];     } ///< Access (read-only) a node
	Element         const * Ele          (size_t i) const { return _elems[i];     } ///< Access (read-only) an element
	size_t                  PushNode     (Node * N)       { _nodes.Push(N); return _nodes.Size()-1; } ///< Push back a new node
	size_t                  PushElem     (Element * E)    { _elems.Push(E); return _elems.Size()-1; } ///< Push back a new element
	size_t                  PushNode     (double X, double Y, double Z=0.0, int Tag=0);               ///< Push back a new node
	size_t                  PushElem     (int Tag, char const * Type, char const * Model, char const * Prms, char const * Inis, char const * Props, bool IsActive, Array<int> const & Connectivity); ///< Push back a new element
	size_t                  GetNode      (double X, double Y, double Z=0.0);        ///< Returns the node ID that matchs the specified coordinates
	Array<Node*>          & Nodes        ()               { return _nodes;        } ///< Access all nodes (read/write)
	Array<Element*>       & Elems        ()               { return _elems;        } ///< Access all elements (read/write)
	Array<Node*>    const & Nodes        ()         const { return _nodes;        } ///< Access all nodes (read-only)
	Array<Element*> const & Elems        ()         const { return _elems;        } ///< Access all elements (read-only)
	Array<Element*> &       ElemsWithTag (int Tag);                                 ///< Return the elements with for a given tag
	void                    Bounds       (double & MinX, double & MinY, double & MaxX, double & MaxY) const;                               ///< Return the limits of the geometry
	void                    Bounds       (double & MinX, double & MinY, double & MinZ, double & MaxX, double & MaxY, double & MaxZ) const; ///< Return the limits of the geometry

#ifdef USE_BOOST_PYTHON
// {
	Node          & PySetNode2D    (size_t i, double X, double Y)                       { return (*SetNode(i,X,Y));   }
	Node          & PySetNode3D    (size_t i, double X, double Y, double Z)             { return (*SetNode(i,X,Y,Z)); }
	PyElem          PySetElem      (size_t i, BPy::str const & Type, bool Act, int Tag) { return PyElem(SetElem(i,BPy::extract<char const *>(Type)(),Act,Tag)); }
	Node    const & PyNod          (size_t i)                                           { return (*Nod(i)); }
	PyElem          PyEle          (size_t i)                                           { return PyElem(Ele(i)); }
	void            PyBounds2D     (BPy::list & MinXY,  BPy::list & MaxXY ) const;
	void            PyBounds3D     (BPy::list & MinXYZ, BPy::list & MaxXYZ) const;
	void            PyElemsWithTag (int Tag, BPy::list & Elems);
	size_t          PyPushNode1    (double X, double Y)                    { return PushNode (X,Y);       }
	size_t          PyPushNode2    (double X, double Y, double Z)          { return PushNode (X,Y,Z);     }
	size_t          PyPushNode3    (double X, double Y, double Z, int Tag) { return PushNode (X,Y,Z,Tag); }
	size_t          PyGetNode1     (double X, double Y)                    { return GetNode  (X,Y);       }
	size_t          PyGetNode2     (double X, double Y, double Z)          { return GetNode  (X,Y,Z);     }
	size_t          PyPushElem     (int Tag, BPy::str const & Type, BPy::str const & Model, BPy::str const & Prms, BPy::str const & Inis, BPy::str const & Props, bool IsActive, BPy::list const & Connectivity);
	void            PyAddLinElems  (BPy::dict const & Edges,  ///< {(n1,n2):tag1, (n3,n4):tag2, ... num edges} n# => node ID
	                                BPy::list const & EAtts); ///< Elements attributes
// }
#endif // USE_BOOST_PYTHON

private:
	// Data
	int                     _dim;             ///< Space dimension
	Array<Node*>            _nodes;           ///< FE nodes
	Array<Element*>         _elems;           ///< FE elements
	Array<Element*>         _beams;           ///< Beams
	Array<int>              _btags;           ///< Beam tags
	std::map<int,size_t>    _elem_tag_idx;    ///< Map Tag => Idx, where Idx is the index inside _elems_with_tags
	Array<Array<Element*> > _elems_with_tags; ///< Element with tags

}; // class Geom


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Geom::~Geom()
{
	for (size_t i=0; i<_nodes.Size(); ++i) if (_nodes[i]!=NULL) delete _nodes[i];
	for (size_t i=0; i<_elems.Size(); ++i) if (_elems[i]!=NULL) delete _elems[i];
}

inline void Geom::SetNNodes(size_t NNodes)
{
	for (size_t i=0; i<_nodes.Size(); ++i) if (_nodes[i]!=NULL) delete _nodes[i];
	_nodes.Resize(NNodes);
	_nodes.SetValues(NULL);
}

inline void Geom::SetNElems(size_t NElems)
{
	for (size_t i=0; i<_elems.Size(); ++i) if (_elems[i]!=NULL) delete _elems[i];
	_elems.Resize(NElems);
	_elems.SetValues(NULL);
}

inline Node * Geom::SetNode(size_t i, double X, double Y, double Z, int Tag)
{
   if (_nodes[i]==NULL) _nodes[i] = new Node;
	_nodes[i]->Initialize (i,X,Y,Z, Tag);
	return _nodes[i];
}

inline Element * Geom::SetElem(size_t i, char const * Type, bool IsActive, int Tag)
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

inline void Geom::ClearDisplacements()
{
	for (size_t i=0; i<_elems.Size(); ++i) _elems[i]->ClearDispAndStrains();
}

inline void Geom::Activate(int ElemTag)
{
	Array<Element*> & elems = ElemsWithTag (ElemTag);
	for (size_t i=0; i<elems.Size(); ++i) elems[i]->SetActive (true);
}

inline void Geom::Deactivate(int ElemTag)
{
	Array<Element*> & elems = ElemsWithTag (ElemTag);
	for (size_t i=0; i<elems.Size(); ++i) elems[i]->SetActive (false);
}

inline bool Geom::Check()
{
	// Check arrays
	if (NNodes()==0 || NElems()==0) return false;

	// Check nodes
	for (size_t i=0; i<NNodes(); ++i) if (_nodes[i]==NULL) return false;

	// Check elements
	for (size_t i=0; i<NElems(); ++i) if (_elems[i]==NULL) return false;

	return true; // OK
}

inline size_t Geom::PushNode(double X, double Y, double Z, int Tag)  
{
	// Allocate new node
	Node * new_node;
	new_node = new Node;
	_nodes.Push(new_node);
	size_t ID = _nodes.Size()-1;
	new_node->Initialize(ID, X, Y, Z, Tag);
	return ID;
}

inline size_t Geom::PushElem(int Tag, char const * Type, char const * Model, char const * Prms, char const * Inis, char const * Props, bool IsActive, Array<int> const & Connectivity )
{
	// Alocate new element
	Element * new_elem;
	new_elem = AllocElement(Type);
	_elems.Push(new_elem);
	size_t ID = _elems.Size()-1;
	new_elem->Initialize(_elems.Size()-1, IsActive, _dim, Tag);

	// Connector (EmbSpring) connectivities
	if (new_elem->NNodes()!=Connectivity.Size()) throw new Fatal("Geom::PushElem: The number of nodes in Connectivity does not match. (Element ID %d)", ID);
	for (size_t i=0; i<new_elem->NNodes(); i++)
		new_elem->Connect(i, _nodes[Connectivity[i]]);

	// Set the model 
	new_elem->SetModel(Model, Prms, Inis);
	new_elem->SetProps(Props);
	return ID;
}

inline size_t Geom::GetNode(double X, double Y, double Z)
{
	for (size_t i=0; i<_nodes.Size(); i++)
	{
		if( _nodes[i]->X()==X && _nodes[i]->Y()==Y && _nodes[i]->Z()==Z) return i;
		std::cout << _nodes[i]->X() << " " << _nodes[i]->Y() << " " << _nodes[i]->Z();
	}
	throw new Fatal("Geom::GetNode: Node not found (%g, %g, %g)", X, Y, Z);
}

inline Array<Element*> & Geom::ElemsWithTag(int Tag)
{
	if (_elem_tag_idx.count(Tag)==0) throw new Fatal("Geom::ElemsWithTag: This Tag==%d was not set for any Element",Tag);
	return _elems_with_tags[_elem_tag_idx[Tag]];
}

inline void Geom::Bounds(double & MinX, double & MinY, double & MaxX, double & MaxY) const
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

inline void Geom::Bounds(double & MinX, double & MinY, double & MinZ, double & MaxX, double & MaxY, double & MaxZ) const
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


/** Outputs a geometry. */
std::ostream & operator<< (std::ostream & os, FEM::Geom const & G)
{
	for (size_t i=0; i<G.NElems(); ++i)
		if (G.Ele(i)!=NULL) os << (*G.Ele(i));
	return os;
}

#ifdef USE_BOOST_PYTHON
// {

inline void Geom::PyBounds2D(BPy::list & MinXY, BPy::list & MaxXY) const
{
	double  minx,miny, maxx,maxy;
	Bounds (minx,miny, maxx,maxy);
	MinXY.append(minx);  MaxXY.append(maxx);
	MinXY.append(miny);  MaxXY.append(maxy);
}

inline void Geom::PyBounds3D(BPy::list & MinXYZ, BPy::list & MaxXYZ) const
{
	double  minx,miny,minz, maxx,maxy,maxz;
	Bounds (minx,miny,minz, maxx,maxy,maxz);
	MinXYZ.append(minx);  MaxXYZ.append(maxx);
	MinXYZ.append(miny);  MaxXYZ.append(maxy);
	MinXYZ.append(minz);  MaxXYZ.append(maxz);
}

inline void Geom::PyElemsWithTag(int Tag, BPy::list & Elems)
{
	Array<FEM::Element*> & elems = ElemsWithTag (Tag);
	for (size_t i=0; i<elems.Size(); ++i)
		Elems.append (PyElem(elems[i]));
}

inline size_t Geom::PyPushElem(int Tag, BPy::str const & Type, BPy::str const & Model, BPy::str const & Prms, BPy::str const & Inis, BPy::str const & Props, bool IsActive, BPy::list const & Connectivity)
{
	size_t nnodes = BPy::len(Connectivity);
	if (nnodes<2) throw new Fatal("Geom::PyPushElem: Number of nodes in the Connectivity list must be greater than 1");
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

inline void Geom::PyAddLinElems(BPy::dict const & Edges, BPy::list const & EAtts)
{
	/* Example:
	 *           
	 *           # Elements attributes
	 *           eatts = [[-1, 'Spring', '', 'ks=%g', 'ZERO', 'gam=20', True]] # tag, type, model, prms, inis, props, active?
	 */

	// Map element tag to index in EAtts list
	int neatts = BPy::len(EAtts);
	if (neatts<1) throw new Fatal("functions.h::PyAddLinElems: EAtts (element attributes) must contain at least one element");
	std::map<int,int> tag2idx; 
	for (int i=0; i<neatts; ++i)
	{
		BPy::list const & lst = BPy::extract<BPy::list>(EAtts[i])();
		tag2idx[BPy::extract<int>(lst[0])()] = i;
		if (BPy::len(EAtts[i])!=7) throw new Fatal("functions.h::PyAddLinElems: Each sublist in EAtts must have 7 items: tag, type, model, prms, inis, props, active?\n\tExample: eatts = [[-1, 'Spring', '', 'ks=1e+12', 'ZERO', 'gam=20', True]]\n\tlen(EAtts[i])==%d is invalid.",BPy::len(EAtts[i]));
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
		if (iter==tag2idx.end()) throw new Fatal("functions.h::PyAddLinElems: Could not find tag < %d > in the list of Element Attributes", tag);
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

// }
#endif // USE_BOOST_PYTHON

}; // namespace FEM

#endif // MECHSYS_FEM_GEOMETRY_H
