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

#ifndef MECHSYS_LINELASTIC_H
#define MECHSYS_LINELASTIC_H

// MechSys
#include "models/equilibmodel.h"
#include "tensors/tensors.h"
#include "util/string.h"
#include "util/util.h"
#include "util/lineparser.h"

typedef char const * Str_t;

const char LINELASTIC_PN[2][8] = {"E", "nu"};

class LinElastic : public EquilibModel
{
public:
	// Destructor
	virtual ~LinElastic () {}

	// Derived methods
	Str_t       Name        () const { return "LinElastic"; }
	size_t      NPrms       () const { return 2; }
	PrmName_t * GetPrmNames () const { return LINELASTIC_PN; }
	void        Initialize  (int GeomIdx, Array<double> const * Prms, Str_t Inis);

private:
	// Data
	Tensor4 _De; ///< Constant tangent stiffness

	// Private methods
	void   _stiff (Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B) const; ///< Tangent or secant stiffness
	double _val   (Str_t Name) const { throw new Fatal("LinElastic::_val The Name==%s is invalid",Name); }                                       ///< Return internal values

}; // class LinElastic


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void LinElastic::Initialize(int GeomIdx, Array<double> const * Prms, Str_t Inis)
{
	// Set geometry index
	_gi = GeomIdx;
	if (_gi<0) throw new Fatal("LinElastic::Initialize: Geometry type:\n\t[0==3D, 1==2D(plane-strain), 2==2D(plane-stress), 3==2D(axis-symmetric)] must be set via SetGeom before calling this method");

	// Parameters
	_prms = Prms;
	double E  = (*_prms)[0];
	double nu = (*_prms)[1];

	// Check
	if (E<=0.0)             throw new Fatal("LinElastic::Initialize: Young modulus (E) must be provided (and positive). E==%f is invalid",E);
	if (nu<0.0 || nu>0.499) throw new Fatal("LinElastic::Initialize: Poisson ratio (nu) must be provided (and in the range: 0 < nu < 0.5). nu==%f is invalid",nu);

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

	// Initial values
	LineParser lp(Inis);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions (names,values);

	// Parse input
	_sig = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	_eps = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	for (size_t i=0; i<names.Size(); i++)
	{
		     if (names[i]=="ZERO")                   break;
		else if (names[i]=="Sx")                     _sig(0) = values[i];
		else if (names[i]=="Sy")                     _sig(1) = values[i];
		else if (names[i]=="Sz")                     _sig(2) = values[i];
		else if (names[i]=="Sxy" || names[i]=="Syx") _sig(3) = values[i]*SQ2;
		else if (names[i]=="Syz" || names[i]=="Szy") _sig(4) = values[i]*SQ2;
		else if (names[i]=="Szx" || names[i]=="Sxz") _sig(5) = values[i]*SQ2;
		else throw new Fatal("LinElastic::Initialize: '%s' component of stress is invalid",names[i].CStr());
	}
}


/* private */

inline void LinElastic::_stiff(Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B) const
{
	D = _De;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new LinElastic model
Model * LinElasticMaker()
{
	return new LinElastic();
}

// Register LinElastic model into ModelFactory array map
int LinearElasticRegister()
{
	ModelFactory["LinElastic"] = LinElasticMaker;
	return 0;
}

// Execute the autoregistration
int __LinearElastic_dummy_int = LinearElasticRegister();


#endif // MECHSYS_LINELASTIC_H
