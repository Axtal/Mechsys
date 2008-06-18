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

#ifndef MECHSYS_FEM_LIN2EQUILIB_H
#define MECHSYS_FEM_LIN2EQUILIB_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/lin2.h"

namespace FEM
{

class Lin2Equilib : public Lin2, public EquilibElem
{
public:
	// Constants
	static String NAME;

	// Derived methods
	String Name() const { return NAME; };

private:
}; // class Lin2Equilib

// Lin2Equilib constants
String Lin2Equilib::NAME = "Lin2Equilib";

///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Lin2Equilib element
Element * Lin2EquilibMaker()
{
	return new Lin2Equilib();
}

// Register a Lin2Equilib element into ElementFactory array map
int Lin2EquilibRegister()
{
	ElementFactory[Lin2Equilib::NAME] = Lin2EquilibMaker;
	return 0;
}

// Execute the autoregistration
int __Lin2Equilib_dummy_int  = Lin2EquilibRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_LIN2EQUILIB_H
