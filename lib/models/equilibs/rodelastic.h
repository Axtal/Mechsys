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

#ifndef MECHSYS_RODELASTIC_H
#define MECHSYS_RODELASTIC_H

// MechSys
#include "models/equilibmodel.h"
#include "tensors/tensors.h"
#include "util/string.h"
#include "util/util.h"
#include "util/lineparser.h"

class RodElastic : public EquilibModel
{
public:
	// Constants
	static const char RODELASTIC_PN[2][8];

	// Destructor
	virtual ~RodElastic () {}

	// Derived methods
	int         NPrms () const { return 2;             }
	PrmName_t * Prms  () const { return RODELASTIC_PN; }
	Str_t       Name  () const { return "RodElastic";  }

private:
	// Private methods
	void _initialize ();
	void _stiff      (Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B) const {}

}; // class RodElastic

const char RodElastic::RODELASTIC_PN[2][8] = {"E", "A"};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void RodElastic::_initialize()
{
	// 3D:gi=0, 2D:gi=1
	if (_gi<0 || _gi>1) throw new Fatal("RodElastic::_initialize: Tag=%d. This model is only available for 3D(gi==0) and 2D(plane-strain, gi==1) problems. (gi==%d is invalid)",_tag,_gi);

	// Parameters
	double E = Prm("E");
	double A = Prm("A");

	// Check
	if (E<=0.0) throw new Fatal("RodElastic::_initialize: Tag=%d. Young modulus (E) must be positive. E==%f is invalid",_tag,E);
	if (A<=0.0) throw new Fatal("RodElastic::_initialize: Tag=%d. Cross sectional area (A) must be positive. A==%f is invalid",_tag,A);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new model
Model * RodElasticMaker() { return new RodElastic(); }

// Register model
int RodElasticRegister() { ModelFactory["RodElastic"]=RodElasticMaker;  return 0; }

// Call register
int __RodElastic_dummy_int = RodElasticRegister();


#endif // MECHSYS_RODELASTIC_H
