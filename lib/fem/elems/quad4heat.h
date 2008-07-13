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

#ifndef MECHSYS_FEM_QUAD4HEAT_H
#define MECHSYS_FEM_QUAD4HEAT_H

// MechSys
#include "fem/heatelem.h"
#include "fem/elems/quad4.h"

namespace FEM
{

class Quad4Heat : public Quad4, public HeatElem
{
public:
	// Constants
	static char const * NAME;

	// Derived methods
	char const * Name() const { return NAME; };

private:
	// Private methods
	int _geom() const { return 2;} ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)

}; // class Quad4Heat

// Quad4Heat constants
char const * Quad4Heat::NAME = "Quad4Heat";

///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Quad4Heat element
Element * Quad4HeatMaker()
{
	return new Quad4Heat();
}

// Register a Quad4Heat element into ElementFactory array map
int Quad4HeatRegister()
{
	ElementFactory[Quad4Heat::NAME] = Quad4HeatMaker;
	return 0;
}

// Execute the autoregistration
int __Quad4Heat_dummy_int  = Quad4HeatRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_QUAD4HEAT_H
