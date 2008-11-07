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

#ifndef MECHSYS_FEM_QUAD4PSTRESS_H
#define MECHSYS_FEM_QUAD4PSTRESS_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/quad4.h"
#include "util/exception.h"

namespace FEM
{

class Quad4PStress : public Quad4, public EquilibElem
{
public:
	// Constants
	static char const * NAME;

	// Derived methods
	char const * Name() const { return NAME; };

private:
	// Private methods
	int  _geom     () const { return 5; } ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void _set_ndim (int nDim)             ///< Set space dimension
	{
		if (nDim<2) throw new Fatal("Quad4PStress::_set_ndim: For this element, nDim must be greater than or equal to 2 (%d is invalid)",nDim);
		_ndim = nDim;
		_d    = _ndim-1;
		_nd   = EquilibElem::ND[_d];
	}

}; // class Quad4PStress

// Quad4PStress constants
char const * Quad4PStress::NAME = "Quad4PStress";


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Quad4PStress element
Element * Quad4PStressMaker()
{
	return new Quad4PStress();
}

// Register a Quad4PStress element into ElementFactory array map
int Quad4PStressRegister()
{
	ElementFactory[Quad4PStress::NAME] = Quad4PStressMaker;
	return 0;
}

// Execute the autoregistration
int __Quad4PStress_dummy_int  = Quad4PStressRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_QUAD4PSTRESS_H
