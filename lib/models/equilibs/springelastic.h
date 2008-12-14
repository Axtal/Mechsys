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

#ifndef MECHSYS_SPRINGELASTIC_H
#define MECHSYS_SPRINGELASTIC_H

// MechSys
#include "models/equilibmodel.h"
#include "tensors/tensors.h"
#include "util/string.h"
#include "util/util.h"
#include "util/lineparser.h"

class SpringElastic : public EquilibModel
{
public:
	// Constants
	static const char SPRINGELASTIC_PN[1][8];

	// Destructor
	virtual ~SpringElastic () {}

	// Derived methods
	int         NPrms () const { return 1;                }
	PrmName_t * Prms  () const { return SPRINGELASTIC_PN; }
	Str_t       Name  () const { return "SpringElastic";  }

private:
	// Private methods
	void _initialize ();
	void _stiff      (Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B) const {}

}; // class SpringElastic

const char SpringElastic::SPRINGELASTIC_PN[1][8] = {"ks"};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void SpringElastic::_initialize()
{
	// 3D:gi=0, 2D:gi=1
	if (_gi<0 || _gi>1) throw new Fatal("SpringElastic::_initialize: Tag=%d. This model is only available for 3D(gi==0) and 2D(plane-strain, gi==1) problems. (gi==%d is invalid)",_tag,_gi);

	// Parameters
	double ks = Prm("ks");

	// Check
	if (ks<=0.0) throw new Fatal("SpringElastic::_initialize: Tag=%d. Spring constant ks must be positive. ks==%f is invalid",_tag,ks);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new model
Model * SpringElasticMaker() { return new SpringElastic(); }

// Register model
int SpringElasticRegister() { ModelFactory["SpringElastic"]=SpringElasticMaker;  return 0; }

// Call register
int __SpringElastic_dummy_int = SpringElasticRegister();


#endif // MECHSYS_SPRINGELASTIC_H
