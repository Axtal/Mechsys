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

#ifndef MECHSYS_LINELASTIC_H
#define MECHSYS_LINELASTIC_H

// MechSys
#include "models/equilibmodel.h"
#include "tensors/tensors.h"
#include "util/string.h"
#include "util/util.h"
#include "util/lineparser.h"

class LinElastic : public EquilibModel
{
public:
	// Constants
	static const char LINELASTIC_PN[2][8];
	static double     LINELASTIC_DP[2];

	// Destructor
	virtual ~LinElastic () {}

	// Derived methods
	int         NPrms   () const { return 2;             }
	PrmName_t * Prms    () const { return LINELASTIC_PN; }
	double    * DefPrms () const { return LINELASTIC_DP; }
	Str_t       Name    () const { return "LinElastic";  }

private:
	// Data
	Tensor4 _De; ///< Constant tangent stiffness

	// Private methods
	void _initialize ();
	void _stiff      (Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs, Tensor2 const & DEps, Tensor4 & D, Array<Tensor2> & B) const { D = _De; }

}; // class LinElastic

const char LinElastic::LINELASTIC_PN[2][8] = {"E", "nu"};
double     LinElastic::LINELASTIC_DP[2]    = {1.0, 0.2};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void LinElastic::_initialize()
{
	// Parameters
	double E  = Prm("E");
	double nu = Prm("nu");

	// Check
	if (E<=0.0)             throw new Fatal("LinElastic::_initialize: Tag=%d: Young modulus (E) must be positive. E==%f is invalid",_tag,E);
	if (nu<0.0 || nu>0.499) throw new Fatal("LinElastic::_initialize: Tag=%d: Poisson ratio (nu) must be in the range: 0 < nu < 0.5. nu==%f is invalid",_tag,nu);

	// Set stiffness
	double c  = (_gi==2 ? E/(1.0-nu*nu)  : E/((1.0+nu)*(1.0-2.0*nu)) ); // (2)plane-stress != (plane-strain=3D)
	double c1 = (_gi==2 ? c*1.0          : c*(1.0-nu)                ); // (2)plane-stress != (plane-strain=3D)
	double c2 = (_gi==2 ? c*0.5*(1.0-nu) : c*(1.0-2.0*nu)/2.0        ); // (2)plane-stress != (plane-strain=3D)
	double c3 = c*nu;
	_De = c1     , c3     , c3     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
	      c3     , c1     , c3     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
	      c3     , c3     , c1     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
	      0.0*SQ2, 0.0*SQ2, 0.0*SQ2, c2 *2.0, 0.0*2.0, 0.0*2.0,
	      0.0*SQ2, 0.0*SQ2, 0.0*SQ2, 0.0*2.0, c2 *2.0, 0.0*2.0,
	      0.0*SQ2, 0.0*SQ2, 0.0*SQ2, 0.0*2.0, 0.0*2.0, c2 *2.0; // In Mandel's basis
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new model
Model * LinElasticMaker() { return new LinElastic(); }

// Register model
int LinElasticRegister() { ModelFactory["LinElastic"]=LinElasticMaker;  return 0; }

// Call register
int __LinElastic_dummy_int = LinElasticRegister();


#endif // MECHSYS_LINELASTIC_H
