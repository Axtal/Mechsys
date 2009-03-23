/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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

#ifndef MECHSYS_LINDIFFUSION_H
#define MECHSYS_LINDIFFUSION_H

// MechSys
#include "models/diffusionmodel.h"
#include "util/string.h"
#include "util/util.h"
#include "util/lineparser.h"

class LinDiffusion : public DiffusionModel
{
public:
	// Constants
	static const char LINDIFFUSION_PN[1][8];

	// Destructor
	virtual ~LinDiffusion () {}

	// Derived methods
	int         NPrms () const { return 1;               }
	PrmName_t * Prms  () const { return LINDIFFUSION_PN; }
	Str_t       Name  () const { return "LinDiffusion";  }

private:
	// Data
	Mat3_t _K; ///< Constant conductivity

	// Private methods
	void _initialize ();
	void _cond       (Vec3_t const & DGra, Vec3_t const & Vel, Vec3_t const & Gra, IntVals const & Ivs,  Mat3_t & D, Array<Vec3_t> & B) const;

}; // class LinDiffusion

const char LinDiffusion::LINDIFFUSION_PN[1][8] = {"k"};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void LinDiffusion::_initialize()
{
	// Parameters
	double k = Prm("k");

	// Check
	if (k<=0.0) throw new Fatal("LinDiffusion::_initialize: Tag=%d: Permeability (k) must be positive. k==%f is invalid",_tag,k);

	// Set conductivity
	_K =   k, 0.0, 0.0,
	     0.0,   k, 0.0,
	     0.0, 0.0,   k;
}

inline void LinDiffusion::_cond(Vec3_t const & DGra, Vec3_t const & Vel, Vec3_t const & Gra, IntVals const & Ivs,  Mat3_t & D, Array<Vec3_t> & B) const
{
	D = _K;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new model
Model * LinDiffusionMaker() { return new LinDiffusion(); }

// Register model
int LinDiffusionRegister() { ModelFactory["LinDiffusion"]=LinDiffusionMaker;  return 0; }

// Call register
int __LinDiffusion_dummy_int = LinDiffusionRegister();


#endif // MECHSYS_LINDIFFUSION_H
