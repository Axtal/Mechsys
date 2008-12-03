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

#ifndef MECHSYS_FEM_QUAD4PSTRAIN_H
#define MECHSYS_FEM_QUAD4PSTRAIN_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/quad4.h"
#include "util/exception.h"

namespace FEM
{

class Quad4PStrain : public Quad4, public EquilibElem
{
public:
	// Constants
	static char const * NAME;

	// Derived methods
	char const * Name() const { return NAME; };

private:
	// Private methods
	int  _geom     () const { return 2; } ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void _initialize()
	{
		if (_ndim<2) throw new Fatal("Quad4PStrain::_initialize: For this element, _ndim must be greater than or equal to 2 (%d is invalid)",_ndim);
		_d  = _ndim-1;
		_nd = EquilibElem::ND[_d];
		_nl = EquilibElem::NL[_geom()-1];
	}

}; // class Quad4PStrain

// Quad4PStrain constants
char const * Quad4PStrain::NAME = "Quad4PStrain";


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Quad4PStrain element
Element * Quad4PStrainMaker()
{
	return new Quad4PStrain();
}

// Register a Quad4PStrain element into ElementFactory array map
int Quad4PStrainRegister()
{
	ElementFactory[Quad4PStrain::NAME] = Quad4PStrainMaker;
	return 0;
}

// Execute the autoregistration
int __Quad4PStrain_dummy_int  = Quad4PStrainRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_QUAD4PSTRAIN_H
