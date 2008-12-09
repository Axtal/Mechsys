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

#ifndef MECHSYS_FEM_TRI6BIOT_H
#define MECHSYS_FEM_TRI6BIOT_H

// MechSys
#include "fem/biotelem.h"
#include "fem/elems/tri6.h"
#include "util/exception.h"

namespace FEM
{

class Tri6Biot : public Tri6, public BiotElem
{
public:
	// Derived methods
	char const * Name() const { return "Tri6Biot"; }

private:
	// Private methods
	int  _geom     () const { return 2; } ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void _initialize()
	{
		if (_ndim<2) throw new Fatal("Tri6Biot::_initialize: For this element, _ndim must be greater than or equal to 2 (%d is invalid)",_ndim);
		_d  = _ndim-1;
		_nd = BiotElem::ND[_d];
		_nl = BiotElem::NL[_geom()-1];
	}

}; // class Tri6Biot


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Tri6Biot element
Element * Tri6BiotMaker()
{
	return new Tri6Biot();
}

// Register a Tri6Biot element into ElementFactory array map
int Tri6BiotRegister()
{
	ElementFactory["Tri6Biot"] = Tri6BiotMaker;
	return 0;
}

// Execute the autoregistration
int __Tri6Biot_dummy_int  = Tri6BiotRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_TRI6BIOT_H
