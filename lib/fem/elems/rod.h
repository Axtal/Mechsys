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

#ifndef MECHSYS_FEM_ROD_H
#define MECHSYS_FEM_ROD_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/lin2.h"
#include "util/exception.h"

namespace FEM
{

class Rod : public Lin2, public EquilibElem
{
public:
	// Constants
	static char const * NAME;

	// Derived methods
	char const * Name() const { return NAME; };

private:
	// Private methods
	int  _geom     () const { return 1; } ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void _set_ndim (int nDim)             ///< Set space dimension
	{
		if (nDim<1) throw new Fatal("Rod::_set_ndim: For this element, nDim must be greater than or equal to 1 (%d is invalid)",nDim);
		_ndim = nDim;
		_d    = _ndim-1;
		_nd   = EquilibElem::ND[_d];
	}

}; // class Rod

// Rod constants
char const * Rod::NAME = "Rod";


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Rod element
Element * RodMaker()
{
	return new Rod();
}

// Register a Rod element into ElementFactory array map
int RodRegister()
{
	ElementFactory[Rod::NAME] = RodMaker;
	return 0;
}

// Execute the autoregistration
int __Rod_dummy_int  = RodRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_LIN2EQUILIB_H
