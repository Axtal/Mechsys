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

#ifndef MECHSYS_TOTO_H
#define MECHSYS_TOTO_H

// MechSys
#include "models/equilibmodel.h"
#include "tensors/tensors.h"
#include "tensors/operators.h"
#include "tensors/functions.h"
#include "util/string.h"
#include "util/util.h"
#include "util/lineparser.h"

using Tensors::AddScaled;
using Tensors::Psd;
using Tensors::IdyI;
using Tensors::IIsym;
using Tensors::Stress_p_q;
using Tensors::Strain_Ev_Ed;

class Toto : public EquilibModel
{
public:
	// Constants
	static const char TOTO_PN[16][8];
	
	// Destructor
	virtual ~Toto () {}

	// Derived methods
	int         NPrms () const { return 16;      }
	PrmName_t * Prms  () const { return TOTO_PN; }
	Str_t       Name  () const { return "Toto";  }

private:
	// Data
	Tensor4 _De; ///< Constant tangent stiffness

	// Private methods
	void _initialize ();
	void _stiff      (Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs, Tensor2 const & DEps, Tensor4 & D, Array<Tensor2> & B) const;

}; // class Toto

const char Toto::TOTO_PN[16][8] = { "k1","l1","b1","psi1",  "k2","l2","b2","psi2",  "k3","l3","b3","ev3",  "k4","l4","b4","ev4" };


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void Toto::_initialize()
{
}

inline void Toto::_stiff(Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs, Tensor2 const & DEps, Tensor4 & D, Array<Tensor2> & B) const
{
	if (false)
	{
		double E  = 1.0e+15;
		double nu = 0.0;
		AddScaled (E/(1.0+nu), IIsym, nu*E/((1.0+nu)*(1.0-2.0*nu)), IdyI, D); // Elastic tangent tensor
		return;
	}
	double p, q, Ev, Ed, psi;
	Stress_p_q   (Sig, p, q);
	Strain_Ev_Ed (Eps, Ev, Ed);
	psi = q/p;
	double F1 = Prm("k1")+(Prm("l1")-Prm("k1")) * exp(-Prm("b1")*( Prm("psi1")+Prm("l1")*Ed-psi ));
	double F2 = Prm("k2")+(Prm("l2")-Prm("k2")) * exp(-Prm("b2")*( Prm("psi2")+Prm("l2")*Ev-psi ));
	double F3 = Prm("k3")+(Prm("l3")-Prm("k3")) * exp(-Prm("b3")*( Ev-Prm("ev3")-Prm("l3")*p    ));
	double F4 = Prm("k4")+(Prm("l4")-Prm("k4")) * exp(-Prm("b4")*( Ev-Prm("ev4")+Prm("l4")*Ed   ));
	double f1 = -p*F1-(p*F2+psi/F3)*F4;
	AddScaled (f1,Psd, 1.0/F3,IdyI, D); // Elastic tangent tensor

	double DEv, DEd;
	Strain_Ev_Ed (DEps, DEv, DEd);

	//std::cout << "DEv = " << DEv << ",  DEd = " << DEd << std::endl;
	//std::cout << D << std::endl;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new model
Model * TotoMaker() { return new Toto(); }

// Register model
int TotoRegister() { ModelFactory["Toto"]=TotoMaker;  return 0; }

// Call register
int __Toto_dummy_int = TotoRegister();


#endif // MECHSYS_TOTO_H
