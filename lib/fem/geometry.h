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

	// Set methods
	void      SetNNodes (size_t nNodes);                                   ///< Set the number of nodes
	void      SetNElems (size_t nElems);                                   ///< Set the number of elements
	Node    * SetNode   (size_t i, double X, double Y, double Z=0.0);      ///< Set a node
	Element * SetElem   (size_t i, char const * Type, bool IsActive=true); ///< Set an element

	// Access methods
	size_t                  NNodes    ()         const { return _nodes.Size(); } ///< Return the number of nodes
	size_t                  NElems    ()         const { return _elems.Size(); } ///< Return the number of elements
	Node                  * Nod       (size_t i)       { return _nodes[i];     } ///< Access (read/write) a node
	Element               * Ele       (size_t i)       { return _elems[i];     } ///< Access (read/write) an element
	Node            const * Nod       (size_t i) const { return _nodes[i];     } ///< Access (read-only) a node
	Element         const * Ele       (size_t i) const { return _elems[i];     } ///< Access (read-only) an element
	Array<Node*>          & Nodes     ()               { return _nodes;        } ///< Access all nodes (read/write)
	Array<Element*>       & Elems     ()               { return _elems;        } ///< Access all elements (read/write)
	Array<Node*>    const & Nodes     ()         const { return _nodes;        } ///< Access all nodes (read-only)
	Array<Element*> const & Elems     ()         const { return _elems;        } ///< Access all elements (read-only)

private:
	// Data
	int             _dim;   ///< Space dimension
	Array<Node*>    _nodes; ///< FE nodes
	Array<Element*> _elems; ///< FE elements

}; // class Geom


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void Geom::SetNNodes(size_t nNodes)
{
	for (size_t i=0; i<_nodes.Size(); ++i) if (_nodes[i]!=NULL) delete _nodes[i];
	_nodes.Resize(nNodes);
	_nodes = NULL;
}

inline void Geom::SetNElems(size_t nElems)
{
	for (size_t i=0; i<_elems.Size(); ++i) if (_elems[i]!=NULL) delete _elems[i];
	_elems.Resize(nElems);
	_elems = NULL;
}


inline Node * Geom::SetNode(size_t i, double X, double Y, double Z)
{
   if (_nodes[i]==NULL) _nodes[i] = new Node;
	_nodes[i]->Initialize (i,X,Y,Z);
	return _nodes[i];
}

inline Element * Geom::SetElem(size_t i, char const * Type, bool IsActive)
{
	if (_elems[i]==NULL) _elems[i] = AllocElement(Type);
	_elems[i]->SetID  (i);
	_elems[i]->SetDim (_dim);
	if (IsActive) _elems[i]->Activate  ();
	else          _elems[i]->Deactivate();
	return _elems[i];
}

}; // namespace FEM


#ifdef USE_BOOST_PYTHON
// {

namespace boopy = boost::python;

class PyGeom
{
public:
	// Constructor
	PyGeom(int nDim) : _geom(nDim) {}

	// Set methods
	void   SetNNodes (size_t nNodes)                                    { _geom.SetNNodes(nNodes); }
	void   SetNElems (size_t nElems)                                    { _geom.SetNElems(nElems); }
	PyNode SetNode   (size_t i, double X, double Y)                     { return PyNode(_geom.SetNode(i,X,Y)  ); }
	PyNode SetNode   (size_t i, double X, double Y, double Z)           { return PyNode(_geom.SetNode(i,X,Y,Z)); }
	PyElem SetElem   (size_t i, boopy::str const & Type)                { return PyElem(_geom.SetElem(i,boopy::extract<char const *>(Type)())         ); }
	PyElem SetElem   (size_t i, boopy::str const & Type, bool IsActive) { return PyElem(_geom.SetElem(i,boopy::extract<char const *>(Type)(),IsActive)); }

	// Access methods
	size_t            NNodes    ()         const { return _geom.NNodes(); }
	size_t            NElems    ()         const { return _geom.NElems(); }
	PyNode            Nod       (size_t i)       { return PyNode(_geom.Nod(i)); }
	PyElem            Ele       (size_t i)       { return PyElem(_geom.Ele(i)); }
	FEM::Geom       * GetGeom   ()               { return &_geom; }
	FEM::Geom const * GetGeom   ()         const { return &_geom; }

private:
	// Data
	FEM::Geom _geom;

}; // class PyGeom 

PyNode (PyGeom::*PGSetNode1)(size_t i, double X, double Y)           = &PyGeom::SetNode;
PyNode (PyGeom::*PGSetNode2)(size_t i, double X, double Y, double Z) = &PyGeom::SetNode;

PyElem (PyGeom::*PGSetElem1)(size_t i, boopy::str const & Type)                = &PyGeom::SetElem;
PyElem (PyGeom::*PGSetElem2)(size_t i, boopy::str const & Type, bool IsActive) = &PyGeom::SetElem;

// }
#endif // USE_BOOST_PYTHON


#endif // MECHSYS_FEM_GEOMETRY_H
