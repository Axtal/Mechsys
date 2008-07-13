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

#ifndef MECHSYS_MODEL_H
#define MECHSYS_MODEL_H

// MechSys
#include "util/string.h"
#include "util/util.h"

class Model
{
public:
	// Destructor
	virtual ~Model () {}

	// Methods that MUST be derived
	virtual void SetPrms (char const * Prms) =0;
	virtual void SetInis (char const * Inis) =0;

	// Access methods that MUST be derived
	virtual double       Val  (char const * Name) const =0; ///< Value: Sx, Sy, Ex, Ey, Wp, z0, etc.
	virtual char const * Name ()                  const =0; ///< Return the name of the constitutive model

}; // class Model


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Define a pointer to a function that makes (allocate) a new Model
typedef Model * (*ModelMakerPtr)();

// Typdef of the array map that contains all the pointers to the functions that makes Models
typedef std::map<String, ModelMakerPtr, std::less<String> > ModelFactory_t;

// Instantiate the array map that contains all the pointers to the functions that makes Models
ModelFactory_t ModelFactory;

// Allocate a new Models according to a string giving the name of the Model
Model * AllocModel(char const * Name)
{
	// Check if there is Name model implemented
	ModelMakerPtr ptr=NULL;
	ptr = ModelFactory[Name];
	if (ptr==NULL)
		throw new Fatal(_("FEM::AllocModel: There is no < %s > implemented in this library"), Name);

	return (*ptr)();
}


#endif // MECHSYS_MODEL_H
