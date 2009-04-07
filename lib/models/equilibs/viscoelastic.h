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

#ifndef MECHSYS_VISCOELASTIC_H
#define MECHSYS_VISCOELASTIC_H

// MechSys
#include "models/equilibmodel.h"
#include "tensors/tensors.h"
#include "util/string.h"
#include "util/util.h"
#include "util/lineparser.h"

class ViscoElastic : public EquilibModel
{
public:
	// Constants
	static const char LINELASTIC_PN[4][8];
	static double     LINELASTIC_DP[4];

	// Destructor
	virtual ~ViscoElastic () {}

	// Derived methods
	int         NPrms   () const { return 4;              }
	PrmName_t * Prms    () const { return LINELASTIC_PN;  }
	double    * DefPrms () const { return LINELASTIC_DP;  }
	Str_t       Name    () const { return "ViscoElastic"; }

	// Delta sigma star (due to viscosity)
	bool HasDSigStar  () const { return true; }
	void CalcDSigStar (double t, double Dt, MechState const & State, Vec_t & DSigStar) const;

private:
	// Data
	Tensor4 _De; ///< Constant tangent stiffness

	// Derived methods
	void _initialize ();
	void _stiff      (Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs, Tensor2 const & DEps, Tensor4 & D, Array<Tensor2> & B) const { D = _De; }
	void _dsig_star  (double Time, double Dt, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs, Tensor2 & DSs) const;

}; // class ViscoElastic

const char ViscoElastic::LINELASTIC_PN[4][8] = {"E", "nu", "alp", "bet"};
double     ViscoElastic::LINELASTIC_DP[4]    = {1.0, 0.2,  1.0,   1.0};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void ViscoElastic::CalcDSigStar(double Time, double Dt, MechState const & State, Vec_t & DSs) const
{
	Tensor2 dss;
	_dsig_star (Time, Dt, State.Sig, State.Eps, State.Ivs, dss);
	Tensor2ToVector (_gi,dss, DSs);
}


/* private */

inline void ViscoElastic::_initialize()
{
	// Parameters
	double E  = Prm("E");
	double nu = Prm("nu");

	// Check
	if (E<=0.0)             throw new Fatal("ViscoElastic::_initialize: Tag=%d: Young modulus (E) must be positive. E==%f is invalid",_tag,E);
	if (nu<0.0 || nu>0.499) throw new Fatal("ViscoElastic::_initialize: Tag=%d: Poisson ratio (nu) must be in the range: 0 < nu < 0.5. nu==%f is invalid",_tag,nu);

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

inline void ViscoElastic::_dsig_star(double Time, double Dt, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs, Tensor2 & DSs) const
{
	double alp      = Prm("alp");
	double bet      = Prm("bet");
	double norm_sig = Tensors::Norm(Sig);

	/*
	DSs = (0.01/norm_sig) * Sig;
	*/

	Tensor2 deps_star;
	deps_star = (alp*exp(-bet*Time)/norm_sig) * Sig;
	Tensors::Dot (_De,deps_star, DSs); // DSs <- De:deps_star
	//std::cout << Time << ", " << deps_star << ", " << DSs << std::endl;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new model
Model * ViscoElasticMaker() { return new ViscoElastic(); }

// Register model
int ViscoElasticRegister() { ModelFactory["ViscoElastic"]=ViscoElasticMaker;  return 0; }

// Call register
int __ViscoElastic_dummy_int = ViscoElasticRegister();


#endif // MECHSYS_VISCOELASTIC_H
