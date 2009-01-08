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


/* __ Finite element NODE __

   The Node structure retains all information related to Essential and Natural values.
   It holds the Increments of Essential/Natural values that lead to the new stage.
   In addition, it holds the current values of essential and natural quantities.

   To setup boundary conditions, first the Element must add the nodes, because, during
   the setup of nodes (connectivities), each element will configure the Essential/
   Natural kinds properly.

 */


#ifndef MECHSYS_FEM_NODE_H
#define MECHSYS_FEM_NODE_H

// STL
#include <sstream>
#include <algorithm>

// Blitz++
#include <blitz/tinyvec-et.h>

// Boost::Python
#ifdef USE_BOOST_PYTHON
  #include <boost/python.hpp> // this includes everything
  namespace BPy = boost::python;
#endif

// MechSys
#include "util/array.h"
#include "util/string.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "fem/node.h"

namespace FEM
{

/** %Node information. */
class Node
{
public:
	// Typedefs
	typedef blitz::TinyVector<double,3> TinyVec;

	/** Degrees of freedom. */
	struct DOF
	{
		String EssentialBryName;  ///< Essential boundary name
		String NaturalBryName;    ///< Natural boundary name
		double EssentialBry;      ///< (U) Applied essential boundary condition for a stage
		double NaturalBry;        ///< (F) Applied natural boundary condition for a stage
		bool   IsEssenPresc;      ///< Is EssentialBry Prescribed?
		long   EqID;              ///< Position inside the system of equations
		double EssentialVal;      ///< Calculated essential value for the current state of this node
		double NaturalVal;        ///< Calculated natural value for the current state of this node
	};

	// Methods
	void   Initialize     (int ID, double X, double Y, double Z, int Tag=0);                   ///< Set the ID of this node and its coordinates
	void   AddDOF         (char const * StrEssentialBry, char const * StrNaturalBry);          ///< Add a degree of freedom
	void   RemoveDOF      (char const * Name);                                                 ///< Remove a degree of freedom
	bool   IsEssential    (char const * Name) const;                                           ///< Check if a degree of freedom is essential
	void   SetSharedBy    (int ElementID);                                                     ///< Set a new element which shares this node
	void   RemoveSharedBy (int ElementID);                                                     ///< Remove an element which shares this node
	bool   HasVar         (char const * Name) const { return (_find_var(Name)<0?false:true); } ///< Check if this node has a variable, such as ux, fx, etc., named Name
	void   ClearBryValues ();                                                                  ///< Clear only boundary information, but does NOT change (calculated) EssentialVal and NaturalVal (U and F values)

	// DOFs access methods
	DOF       & DOFVar (char const * Name);                       ///< Access a DOF structure by name (read/write)
	DOF       & DOFVar (int Index)       { return _dofs[Index]; } ///< Access a DOF structure by index (read/write)
	DOF const & DOFVar (int Index) const { return _dofs[Index]; } ///< Access a DOF structure by index (read-only)

	// Access methods
	long   ID        ()          const { return _my_id; }            ///< Return the ID of this node
	int    Tag       ()          const { return _tag;   }            ///< Return the tag of this node
	double Coord     (size_t i)  const { return _coords(i); }        ///< Return coordinate x, y, or z
	double X         ()          const { return _coords(0); }        ///< Return x coordinate
	double Y         ()          const { return _coords(1); }        ///< Return y coordinate
	double Z         ()          const { return _coords(2); }        ///< Return z coordinate
	size_t nDOF      ()          const { return _dofs.Size();      } ///< Return the number of degrees of freedom
	size_t nSharedBy ()          const { return _shared_by.Size(); } ///< Return the array with the elements that share this node
	long   SharedBy  (int Index) const { return _shared_by[Index]; } ///< Return the array with the elements that share this node
	double Val       (char const * Name) const;                      ///< Return the essential or natural value (computed/current). Ex.: Name="ux", "fx", etc.

	// Set methods
	Node * Bry(const char * DOFName, double Value); ///< Set boundary value. If it is Natural, it will be accumulated. If it is Essential, it will just be set.

#ifdef USE_BOOST_PYTHON
// {
	Node & PyBry (BPy::str const & Name, double Value) { return (*Bry(BPy::extract<char const *>(Name)(), Value)); }
	double PyVal (BPy::str const & Name)               { return Val(BPy::extract<char const *>(Name)()); }
// }
#endif // USE_BOOST_PYTHON

private:
	// Data
	long        _my_id;     ///< The ID of this node
	int         _tag;       ///< The tag of this node
	Array<long> _shared_by; ///< IDs of the elements that share this node
	TinyVec     _coords;    ///< Coordinates
	Array<DOF>  _dofs;      ///< Array with degrees of freedom

	// Private methods
	long _find_var(char const * Name) const;

}; // class Node


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Methods

inline void Node::Initialize(int ID, double X, double Y, double Z, int Tag)
{
	_my_id  = ID;
	_coords = X, Y, Z;
	_tag    = Tag;
}

inline void Node::AddDOF(char const * EssentialBryName, char const * NaturalBryName)
{
	if (HasVar(EssentialBryName)==false) // not added yet
	{
		DOF tmp = { EssentialBryName,NaturalBryName,0.0,0.0,false,-1,0.0,0.0 };
		_dofs.Push(tmp);
	}
}

inline void Node::RemoveDOF(char const * Name)
{
	long idx = _find_var(Name); // Index in _dofs array corroesponding to the variable Name
	if (idx==-1)
		throw new Fatal(_("Node::RemoveDOF: Could not find DOF variable name < %s > inside Node (%d)"), Name, _my_id);
	else 
		_dofs.Remove(idx);
}

inline bool Node::IsEssential(char const * Name) const
{
	long idx = _find_var(Name);
	if (idx<0) // not added
		throw new Fatal(_("Node::IsEssential: Could not find DOF variable name < %s > inside Node (%d)"), Name, _my_id);
	return _dofs[idx].EssentialBryName==Name;
}

inline Node::DOF & Node::DOFVar(char const * Name)
{
	long idx = _find_var(Name);
	if (idx<0) // not added
		throw new Fatal(_("Node::DOFVar: Could not find DOF variable name < %s > inside Node (%d)"), Name, _my_id);
	return _dofs[idx];
}

inline double Node::Val(char const * Name) const
{
	long idx = _find_var(Name);
	if (idx<0) // not added
		throw new Fatal(_("Node::DOFVar: Could not find DOF variable name < %s > inside Node (%d)"), Name, _my_id);
	if (_dofs[idx].EssentialBryName==Name) return _dofs[idx].EssentialVal;
	else                                   return _dofs[idx].NaturalVal;
}

inline void Node::SetSharedBy(int ElementID)
{
	long idx = _shared_by.Find(ElementID);
	if (idx<0) // not added yet
		_shared_by.Push(ElementID);
}

inline void Node::RemoveSharedBy(int ElementID)
{
	long idx = _shared_by.Find(ElementID);
	if (idx<0) // not added
		throw new Fatal(_("Node::RemoveSharedBy: Node (%d) does not share Element (%d)"),_my_id,ElementID);
	_shared_by.Remove(idx);
}

inline void Node::ClearBryValues()
{
	for (size_t i=0; i<_dofs.Size(); ++i)
	{
		// Default => (ux=?,fx=0,False,-1)
		_dofs[i].EssentialBry = 0.0;
		_dofs[i].NaturalBry   = 0.0;
		_dofs[i].IsEssenPresc = false;
		_dofs[i].EqID         = -1;
	}
}

// Set methods
inline Node * Node::Bry(const char * DOFName, double Value)
{
	if (_shared_by.Size()==0) return this;
	long idx = _find_var(DOFName);
	if (idx<0) // not added
		throw new Fatal(_("Node::Bry: Could not find DOF variable name < %s > inside Node (ID=%d, X=%f, Y=%f, Z=%f)"), DOFName, ID(), X(), Y(), Z());
	if (_dofs[idx].EssentialBryName==DOFName) // is essential
	{
		_dofs[idx].EssentialBry = Value;
		_dofs[idx].IsEssenPresc = true;
	}
	else _dofs[idx].NaturalBry += Value; // Bry function must not set IsEssenPresc in this case
	return this;
}

// Private
inline long Node::_find_var(char const * Name) const
{
	long found = -1;
	for (size_t i=0; i<_dofs.Size(); i++)
		if (Name==_dofs[i].EssentialBryName || Name==_dofs[i].NaturalBryName) { found=i; break; }
	return found;
}

// operator <<

/** Outputs a DOF. */
std::ostream & operator<< (std::ostream & os, FEM::Node::DOF const & D)
{
	os << "{" << D.EssentialBryName << ","
	          << D.NaturalBryName   << ","
	          << D.EssentialBry     << ","
	          << D.NaturalBry       << ","
	          << D.IsEssenPresc     << ","
	          << D.EqID             << ","
	          << D.EssentialVal     << ","
	          << D.NaturalVal       << "}";
	return os;
}

/** Outputs a node. */
std::ostream & operator<< (std::ostream & os, FEM::Node const & N)
{
	os << "[" << N.ID() << "] (" << Util::_8_4<<N.Coord(0) << "," << Util::_8_4<<N.Coord(1) << "," << Util::_8_4<<N.Coord(2) << ")";
	for (size_t i=0; i<N.nDOF(); ++i)
		os << " : " << N.DOFVar(i);
	return os;
}

}; //namespace FEM

#endif // MECHSYS_FEM_NODE_H
