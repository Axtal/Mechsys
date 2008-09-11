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

#ifndef MECHSYS_FEM_TRI3HEAT_H
#define MECHSYS_FEM_TRI3HEAT_H

// MechSys
#include "fem/heatelem.h"
#include "fem/elems/tri3.h"

namespace FEM
{

class Tri3Heat : public Tri3, public HeatElem
{
public:
	// Constants
	static char const * NAME;

	// Derived methods
	char const * Name() const { return NAME; };

private:
	// Private methods
	int _geom() const { return 2; } ///< Geometry of the element: 1:1D, 2:2D, 3:3D

}; // class Tri3Heat

// Tri3Heat constants
char const * Tri3Heat::NAME = "Tri3Heat";

///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Tri3Heat element
Element * Tri3HeatMaker()
{
	return new Tri3Heat();
}

// Register a Tri3Heat element into ElementFactory array map
int Tri3HeatRegister()
{
	ElementFactory[Tri3Heat::NAME] = Tri3HeatMaker;
	return 0;
}

// Execute the autoregistration
int __Tri3Heat_dummy_int  = Tri3HeatRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_TRI3HEAT_H
