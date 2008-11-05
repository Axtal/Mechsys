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

#ifndef MECHSYS_FEM_TRI3DIFFUSION_H
#define MECHSYS_FEM_TRI3DIFFUSION_H

// MechSys
#include "fem/diffusionelem.h"
#include "fem/elems/tri3.h"
#include "util/exception.h"

namespace FEM
{

class Tri3Diffusion : public Tri3, public DiffusionElem
{
public:
	// Constants
	static char const * NAME;

	// Derived methods
	char const * Name() const { return NAME; };

private:
	// Private methods
	int  _geom     () const { return 2; } ///< Geometry of the element: 1:1D, 2:2D, 3:3D
	void _set_ndim (int nDim)             ///< Set space dimension
	{
		if (nDim<2) throw new Fatal("Tri3Diffusion::_set_ndim: For this element, nDim must be greater than or equal to 2 (%d is invalid)",nDim);
		_ndim = nDim;
	}

}; // class Tri3Diffusion

// Tri3Diffusion constants
char const * Tri3Diffusion::NAME = "Tri3Diffusion";


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Tri3Diffusion element
Element * Tri3DiffusionMaker()
{
	return new Tri3Diffusion();
}

// Register a Tri3Diffusion element into ElementFactory array map
int Tri3DiffusionRegister()
{
	ElementFactory[Tri3Diffusion::NAME] = Tri3DiffusionMaker;
	return 0;
}

// Execute the autoregistration
int __Tri3Diffusion_dummy_int  = Tri3DiffusionRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_TRI3DIFFUSION_H
