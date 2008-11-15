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

class LinElastic : public EquilibModel
{
public:
	// Destructor
	virtual ~LinElastic () {}

	// Derived Methods
	void         SetPrms (char const * Prms);
	void         SetInis (char const * Inis);
	char const * Name    () const { return "LinElastic"; }

private:
	// Data
	Tensor4 _De; ///< Constant tangent stiffness

	// Private methods
	void   _stiff (Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B) const; ///< Tangent or secant stiffness
	double _val   (char const * Name) const { throw new Fatal("LinElastic::_val The Name==%s is invalid",Name); }                                ///< Return internal values

}; // class LinElastic


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void LinElastic::SetPrms(char const * Prms)
{
	if (_geom<0) throw new Fatal("LinElastic::SetPrms: Geometry type:\n\t[1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)] must be set via SetGeom before calling this method");

	/* "E=20000.0 nu=0.2" */
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Set
	double E  = -1.0;
	double nu = -1.0;
	for (size_t i=0; i<names.Size(); ++i)
	{
			 if (names[i]=="E" ) E  = values[i];
		else if (names[i]=="nu") nu = values[i];
		else if (names[i]=="A" ) _A = values[i];
	}
	if (_geom==1)
	{
		_De      = 0.0;
		_De(0,0) = E*_A;
	}
	else
	{
		double c  = (_geom==5 ? E/(1.0-nu*nu)  : E/((1.0+nu)*(1.0-2.0*nu)) ); // plane-stress != (plane-strain=3D)
		double c1 = (_geom==5 ? c*1.0          : c*(1.0-nu)                ); // plane-stress != (plane-strain=3D)
		double c2 = (_geom==5 ? c*0.5*(1.0-nu) : c*(1.0-2.0*nu)/2.0        ); // plane-stress != (plane-strain=3D)
		double c3 = c*nu;
		_De = c1     , c3     , c3     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
		      c3     , c1     , c3     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
		      c3     , c3     , c1     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
		      0.0*SQ2, 0.0*SQ2, 0.0*SQ2, c2 *2.0, 0.0*2.0, 0.0*2.0,
		      0.0*SQ2, 0.0*SQ2, 0.0*SQ2, 0.0*2.0, c2 *2.0, 0.0*2.0,
		      0.0*SQ2, 0.0*SQ2, 0.0*SQ2, 0.0*2.0, 0.0*2.0, c2 *2.0; // In Mandel's basis
	}
	if (E <=0.0) throw new Fatal("LinElastic::SetPrms: Young modulus (E) must be provided (and positive). E==%f is invalid",E);
	if (nu<=0.0 || nu>0.499999999) throw new Fatal("LinElastic::SetPrms: Poisson ratio (nu) must be provided (and in the range: 0 < nu < 0.5). nu==%f is invalid",nu);
}

inline void LinElastic::SetInis(char const * Inis)
{
	/* "Sx=0.0 Sy=0.0 Sxy=0.0 ..." or "ZERO" */
	LineParser lp(Inis);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

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
		else throw new Fatal("LinElastic::SetInis: '%s' component of stress is invalid",names[i].CStr());
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
