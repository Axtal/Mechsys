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
#include <map>
#include <sstream>
#include <algorithm>

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
	/** Degrees of freedom. */
	struct DOF
	{
		String EssentialBryName;  ///< Essential boundary name
		String NaturalBryName;    ///< Natural boundary name
		double EssentialBry;      ///< (U) Applied essential boundary condition for a stage
		double NaturalBry;        ///< (F) Applied natural boundary condition for a stage
		bool   IsEssenPresc;      ///< Is EssentialBry Prescribed?
		int    EqID;              ///< Position inside the system of equations
		double EssentialVal;      ///< Calculated essential value for the current state of this node
		double NaturalVal;        ///< Calculated natural value for the current state of this node
	};

	// Methods
	void   Initialize     (int ID, double X, double Y, double Z);                              ///< Set the ID of this node and its coordinates
	void   AddDOF         (char const * StrEssentialBry, char const * StrNaturalBry);          ///< TODO
	bool   IsEssential    (char const * Name) const;                                           ///< TODO
	void   SetSharedBy    (int ElementID);                                                     ///< Set a new element which shares this node
	void   RemoveSharedBy (int ElementID);                                                     ///< Remove an element which shares this node
	bool   HasVar         (char const * Name) const { return (_find_var(Name)<0?false:true); } ///< TODO
	void   ClearBryValues ();                                                                  ///< Clear only boundary information, but does NOT change (calculated) EssentialVal and NaturalVal (U and F values)

	// DOFs access methods
	DOF & DOFVar(char const * Name);
	DOF & DOFVar(int Index) { return _dofs[Index]; }

	// Access methods
	int    GetID     ()          const { return _my_id; }            ///< Return the ID of this node
	double X         ()          const { return _x;     }            ///< X coordinate
	double Y         ()          const { return _y;     }            ///< Y coordinate
	double Z         ()          const { return _z;     }            ///< Z coordinate
	size_t nDOF      ()          const { return _dofs.Size();      } ///< TODO
	size_t nSharedBy ()          const { return _shared_by.Size(); } ///< Return the array with the elements that share this node
	int    SharedBy  (int Index) const { return _shared_by[Index]; } ///< Return the array with the elements that share this node

	// Set methods
	Node * Bry(const char * DOFName, double Value);

private:
	// Data
	int        _my_id;     ///< The ID of this node
	Array<int> _shared_by; ///< IDs of the elements that share this node
	double     _x;         ///< X coordinate
	double     _y;         ///< Y coordinate
	double     _z;         ///< Z coordinate
	Array<DOF> _dofs;      ///< TODO

	// Private methods
	long _find_var(char const * Name) const;

}; // class Node

Array<Node*> Nodes; ///< Array with all nodes


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Methods

inline void Node::Initialize(int ID, double X, double Y, double Z)
{
	_my_id = ID;
	_x     = X;
	_y     = Y;
	_z     = Z;
}

inline void Node::AddDOF(char const * EssentialBryName, char const * NaturalBryName)
{
	if (HasVar(EssentialBryName)==false) // not added yet
	{
		DOF tmp = { EssentialBryName,NaturalBryName,0.0,0.0,false,-1,0.0,0.0 };
		_dofs.Push(tmp);
	}
}

inline bool Node::IsEssential(char const * Name) const
{
	long idx = _find_var(Name);
	if (idx<0) // not added
		throw new Fatal(_("Node::IsEssential: Could not find DOF variable name < %s > inside Node"), Name);
	return _dofs[idx].EssentialBryName==Name;
}

inline Node::DOF & Node::DOFVar(char const * Name)
{
	long idx = _find_var(Name);
	if (idx<0) // not added
		throw new Fatal(_("Node::DOFVar: Could not find DOF variable name < %s > inside Node"), Name);
	return _dofs[idx];
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
	long idx = _find_var(DOFName);
	if (idx<0) // not added
		throw new Fatal(_("Node::Bry: Could not find DOF variable name < %s > inside Node"), DOFName);
	if (_dofs[idx].EssentialBryName==DOFName) // is essential
	{
		_dofs[idx].EssentialBry = Value;
		_dofs[idx].IsEssenPresc = true;
	}
	else
	{
		_dofs[idx].NaturalBry   = Value;
		_dofs[idx].IsEssenPresc = false;
	}
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

/** Outputs nodes. */
std::ostream & operator<< (std::ostream & os, Array<FEM::Node> const & nodes)
{
	os << "Number of nodes = " << nodes.Size() << std::endl;
	os << Util::_6 << "Name" << Util::_8s << "EssentialBry" << Util::_8s << "NaturalBry" << Util::_a << "Epresc?" << Util::_6 << "EqID" << Util::_8s << "EssentialVal" << Util::_8s << "NaturalVal" << std::endl;
	for (size_t i=0; i<nodes.Size(); ++i)
	{
		os << "Node #" << i << " X=" << nodes[i].X() << " Y=" << nodes[i].Y() << " Z=" << nodes[i].Z() << std::endl;
	}
	return os;
}

}; //namespace FEM

#endif // MECHSYS_FEM_NODE_H
