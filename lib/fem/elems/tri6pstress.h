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

#ifndef MECHSYS_FEM_TRI6PSTRESS_H
#define MECHSYS_FEM_TRI6PSTRESS_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/tri6.h"
#include "util/exception.h"

namespace FEM
{

class Tri6PStress : public Tri6, public EquilibElem
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
		if (nDim<2) throw new Fatal("Tri6PStress::_set_ndim: For this element, nDim must be greater than or equal to 2 (%d is invalid)",nDim);
		_ndim = nDim;
		_d    = _ndim-1;
		_nd   = EquilibElem::ND[_d];
	}

}; // class Tri6PStress

// Tri6PStress constants
char const * Tri6PStress::NAME = "Tri6PStress";


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Tri6PStress element
Element * Tri6PStressMaker()
{
	return new Tri6PStress();
}

// Register a Tri6PStress element into ElementFactory array map
int Tri6PStressRegister()
{
	ElementFactory[Tri6PStress::NAME] = Tri6PStressMaker;
	return 0;
}

// Execute the autoregistration
int __Tri6PStress_dummy_int  = Tri6PStressRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_TRI6PSTRESS_H
