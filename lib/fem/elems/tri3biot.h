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

#ifndef MECHSYS_FEM_QUAD4CONSOLID_H
#define MECHSYS_FEM_QUAD4CONSOLID_H

// MechSys
#include "fem/biotelem.h"
#include "fem/elems/tri3.h"
#include "util/exception.h"

namespace FEM
{

class Tri3Biot : public Tri3, public BiotElem
{
public:
	// Derived methods
	char const * Name() const { return "Tri3Biot"; }

private:
	// Private methods
	int  _geom     () const { return 2; } ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void _set_ndim (int nDim)             ///< Set space dimension
	{
		if (nDim<2) throw new Fatal("Tri3Biot::_set_ndim: For this element, nDim must be greater than or equal to 2 (%d is invalid)",nDim);
		_ndim = nDim;
		_d    = _ndim-1;
		_nd   = BiotElem::ND[_d];
	}

}; // class Tri3Biot


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Tri3Biot element
Element * Tri3BiotMaker()
{
	return new Tri3Biot();
}

// Register a Tri3Biot element into ElementFactory array map
int Tri3BiotRegister()
{
	ElementFactory["Tri3Biot"] = Tri3BiotMaker;
	return 0;
}

// Execute the autoregistration
int __Tri3Biot_dummy_int  = Tri3BiotRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_QUAD4CONSOLID_H
