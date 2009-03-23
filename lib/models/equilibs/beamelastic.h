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

#ifndef MECHSYS_BEAMELASTIC_H
#define MECHSYS_BEAMELASTIC_H

// MechSys
#include "models/equilibmodel.h"
#include "tensors/tensors.h"
#include "util/string.h"
#include "util/util.h"
#include "util/lineparser.h"

class BeamElastic : public EquilibModel
{
public:
	// Constants
	static const char BEAMELASTIC_PN[3][8];

	// Destructor
	virtual ~BeamElastic () {}

	// Derived methods
	int         NPrms () const { return 3;              }
	PrmName_t * Prms  () const { return BEAMELASTIC_PN; }
	Str_t       Name  () const { return "BeamElastic";  }

private:
	// Private methods
	void _initialize ();
	void _stiff      (Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B, bool First) const {}

}; // class BeamElastic

const char BeamElastic::BEAMELASTIC_PN[3][8] = {"E", "A", "Izz"};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void BeamElastic::_initialize()
{
	// 3D:gi=0, 2D:gi=1
	if (_gi<0 || _gi>1) throw new Fatal("BeamElastic::_initialize: Tag=%d. This model is only available for 3D(gi==0) and 2D(plane-strain, gi==1) problems. (gi==%d is invalid)",_tag,_gi);

	// Parameters
	double E   = Prm("E");
	double A   = Prm("A");
	double Izz = Prm("Izz");

	// Check
	if (E<=0.0)   throw new Fatal("BeamElastic::_initialize: Tag=%d. Young modulus (E) must be positive). E==%f is invalid",_tag,E);
	if (A<=0.0)   throw new Fatal("BeamElastic::_initialize: Tag=%d. Cross sectional area (A) must be positive). A==%f is invalid",_tag,A);
	if (Izz<=0.0) throw new Fatal("BeamElastic::_initialize: Tag=%d. Inertia (Izz) must be positive). Izz==%f is invalid",_tag,Izz);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new model
Model * BeamElasticMaker() { return new BeamElastic(); }

// Register model
int BeamElasticRegister() { ModelFactory["BeamElastic"]=BeamElasticMaker;  return 0; }

// Call register
int __BeamElastic_dummy_int = BeamElasticRegister();


#endif // MECHSYS_BEAMELASTIC_H
