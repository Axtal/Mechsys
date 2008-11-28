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

#ifndef MECHSYS_FEM_QUAD4BIOT_H
#define MECHSYS_FEM_QUAD4BIOT_H

// MechSys
#include "fem/biotelem.h"
#include "fem/elems/quad4.h"
#include "util/exception.h"

namespace FEM
{

class Quad4Biot : public Quad4, public BiotElem
{
public:
	// Derived methods
	char const * Name() const { return "Quad4Biot"; }

private:
	// Private methods
	int  _geom     () const { return 2; } ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void _set_ndim (int nDim)             ///< Set space dimension
	{
		if (nDim<2) throw new Fatal("Quad4Biot::_set_ndim: For this element, nDim must be greater than or equal to 2 (%d is invalid)",nDim);
		_ndim = nDim;
		_d    = _ndim-1;
		_nd   = BiotElem::ND[_d];
	}

}; // class Quad4Biot


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Quad4Biot element
Element * Quad4BiotMaker()
{
	return new Quad4Biot();
}

// Register a Quad4Biot element into ElementFactory array map
int Quad4BiotRegister()
{
	ElementFactory["Quad4Biot"] = Quad4BiotMaker;
	return 0;
}

// Execute the autoregistration
int __Quad4Biot_dummy_int  = Quad4BiotRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_QUAD4BIOT_H
