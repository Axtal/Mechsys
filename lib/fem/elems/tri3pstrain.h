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

#ifndef MECHSYS_FEM_TRI3PSTRAIN_H
#define MECHSYS_FEM_TRI3PSTRAIN_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/tri3.h"
#include "util/exception.h"

namespace FEM
{

class Tri3PStrain : public Tri3, public EquilibElem
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
		if (_ndim<2) throw new Fatal("Tri3PStrain::_initialize: For this element, _ndim must be greater than or equal to 2 (%d is invalid)",_ndim);
		_d  = _ndim-1;
		_nd = EquilibElem::ND[_d];
		_nl = EquilibElem::NL[_geom()-1];
	}

}; // class Tri3PStrain

// Tri3PStrain constants
char const * Tri3PStrain::NAME = "Tri3PStrain";


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Tri3PStrain element
Element * Tri3PStrainMaker()
{
	return new Tri3PStrain();
}

// Register a Tri3PStrain element into ElementFactory array map
int Tri3PStrainRegister()
{
	ElementFactory[Tri3PStrain::NAME] = Tri3PStrainMaker;
	return 0;
}

// Execute the autoregistration
int __Tri3PStrain_dummy_int  = Tri3PStrainRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_TRI3PSTRAIN_H
