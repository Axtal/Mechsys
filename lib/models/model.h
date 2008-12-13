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
#include "util/array.h"
#include "util/util.h"

typedef char const * Str_t;

class Model
{
public:
	// Typedefs
	typedef const char    PrmName_t[8]; ///< Parameters names. Ex: "E", "nu", "lam", "kap", ...
	typedef Array<double> IntVals;      ///< Internal values (specific volume, yield surface size, etc.)

	// Constructor
	Model () : _gi(-1), _prms(NULL) { STOL().dTini().mMin().mMax().maxSS(); }

	// Destructor
	virtual ~Model () {}

	// Integration constants
	Model & STOL  (double Val=1.0e-5) { _STOL =Val; return (*this); }
	Model & dTini (double Val=1.0   ) { _dTini=Val; return (*this); }
	Model & mMin  (double Val=0.1   ) { _mMin =Val; return (*this); }
	Model & mMax  (double Val=10.0  ) { _mMax =Val; return (*this); }
	Model & maxSS (size_t Val=2000  ) { _maxSS=Val; return (*this); }

	// Methods
	void BackupState  () { _ivs_bkp=_ivs; _backup_state (); } ///< Backup internal state
	void RestoreState () { _ivs=_ivs_bkp; _restore_state(); } ///< Restore internal state

	// Methods that MUST be derived
	virtual Str_t       Name        () const                                              =0; ///< Model name
	virtual size_t      NPrms       () const                                              =0; ///< Return the number of parameters
	virtual PrmName_t * GetPrmNames () const                                              =0; ///< Get array with parameter names
	virtual void        Initialize  (int GeomIdx, Array<double> const * Prms, Str_t Inis) =0; ///< Initialize model
	virtual double      Val         (Str_t Name) const                                    =0; ///< Value: Sx, Sy, Ex, Ey, Wp, z0, etc.

protected:
	// Data
	int                   _gi;      ///< Geometry index: Equilib: 3D=0, PStrain=1, PStress=2, Axis=3. Others: 3D=0, 2D=1
	Array<double> const * _prms;    ///< Parameters
	IntVals               _ivs;     ///< Internal values
	IntVals               _ivs_bkp; ///< Backup internal values

	// Constants for the stress update algorithm
	double _STOL;
	double _dTini;
	double _mMin;
	double _mMax;
	size_t _maxSS;

	// Private methods that MUST be derived
	virtual void _backup_state  () =0;
	virtual void _restore_state () =0;

}; // class Model


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Define a pointer to a function that makes (allocate) a new Model
typedef Model * (*ModelMakerPtr)();

// Typdef of the array map that contains all the pointers to the functions that makes Models
typedef std::map<String, ModelMakerPtr, std::less<String> > ModelFactory_t;

// Instantiate the array map that contains all the pointers to the functions that makes Models
ModelFactory_t ModelFactory;

// Allocate a new Models according to a string giving the name of the Model
Model * AllocModel(Str_t Name)
{
	// Check if there is Name model implemented
	ModelMakerPtr ptr=NULL;
	ptr = ModelFactory[Name];
	if (ptr==NULL)
		throw new Fatal(_("FEM::AllocModel: There is no < %s > implemented in this library"), Name);

	return (*ptr)();
}


#endif // MECHSYS_MODEL_H
