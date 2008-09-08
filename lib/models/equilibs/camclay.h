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

#ifndef MECHSYS_CAMCLAY_H
#define MECHSYS_CAMCLAY_H

// MechSys
#include "models/equilibmodel.h"
#include "tensors/tensors.h"
#include "util/string.h"
#include "util/util.h"
#include "util/lineparser.h"

class CamClay : public EquilibModel
{
public:
	// Destructor
	virtual ~CamClay () {}

	// Derived Methods
	void         SetPrms (char const * Prms);
	void         SetInis (char const * Inis);
	char const * Name    () const { return "CamClay"; }

private:
	// Private data
	double _lam;
	double _kap;
	double _Mcs;
	double _w;
	double _G;

	// Private methods
	void   _stiffness (Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B) const;
	double _calc_M    (double const & t)    const;
	double _calc_z0   (Tensor2 const & Sig) const;
	double _val       (char const * Name)   const;

}; // class CamClay


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void CamClay::SetPrms(char const * Prms)
{
	if (_geom<0) throw new Fatal("CamClay::SetPrms: Geometry type:\n\t[1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)] must be set via SetGeom before calling this method");

	/* "lam=0.1 kap=0.01 phics=30 G=100" */
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Set
	if (names.Size()!=4) throw new Fatal("CamClay::SetPrms: The number of paramaters (%d) is incorrect (must be equal to 4). Ex.: lam=0.1 kap=0.01 phics=30 G=100",names.Size());
	double sinphi;
	for (size_t i=0; i<names.Size(); ++i)
	{
			 if (names[i]=="lam")    _lam   =     values[i];
		else if (names[i]=="kap")    _kap   =     values[i];
		else if (names[i]=="phics")  sinphi = sin(values[i]*Util::PI/180.0);
		else if (names[i]=="G")      _G     =     values[i];
		else throw new Fatal("CamClay::SetPrms: Parameter label==%s is invalid. Ex.: lam=0.1 kap=0.01 phics=30 G=100",names[i].CStr());
	}
	_Mcs = 6.0*sinphi/(3.0-sinphi);
	_w   = pow((3.0-sinphi)/(3.0+sinphi),4.0);
}

inline void CamClay::SetInis(char const * Inis)
{
	/* "Sx=0.0 Sy=0.0 Sxy=0.0 v=2.2" */
	LineParser lp(Inis);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Check
	_sig = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	_eps = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	_ivs.Resize(2); // z0 (yield surface size), and v (specific volume)
	for (size_t i=0; i<names.Size(); i++)
	{
		     if (names[i]=="ZERO")                   break;
		else if (names[i]=="Sx")                     _sig(0) = values[i];
		else if (names[i]=="Sy")                     _sig(1) = values[i];
		else if (names[i]=="Sz")                     _sig(2) = values[i];
		else if (names[i]=="Sxy" || names[i]=="Syx") _sig(3) = values[i]*SQ2;
		else if (names[i]=="Syz" || names[i]=="Szy") _sig(4) = values[i]*SQ2;
		else if (names[i]=="Szx" || names[i]=="Sxz") _sig(5) = values[i]*SQ2;
		else if (names[i]=="v")                      _ivs[1] = values[i];
		else throw new Fatal("CamClay::SetInis: '%s' component of stress is invalid",names[i].CStr());
	}
	_ivs[0] = _calc_z0 (_sig);
}


/* private */

inline void CamClay::_stiffness(Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B) const
{
	// Intersection (elastic/elastoplastic transition)
	double alp = 0;

	// Stiffness
	if (alp+1.0e+5>=1.0) // pure elastic
	{
		double p = (Sig(0)+Sig(1)+Sig(2))/3.0;
		double K = p*Ivs[1]/_kap; // Bulk modulus
		Tensors::AddScaled (2.0*_G,Tensors::Psd, K,Tensors::IdyI, D);
		B.SetValues (0.0);
	}
	else
	{
		// Constants
		const double MIN_Q = 1.0e-5;

		// Invariants
		double   p,q,t;
		Tensor2  S;
		Tensors::Stress_p_q_S_t (Sig, p,q,S,t);

		// Gradients
		Tensor2 r;  // Plastic strain direction
		double y0;  // Derivative of f w.r.t z0
		double z0, M, MM;
		z0 = Ivs[0];
		M  = _calc_M (t);
		MM = M*M;
		r  = (MM*(2.0*p-z0)/3.0)*Tensors::I;
		if (q>MIN_Q)
		{
			r += 3.0*S;
		}
		y0 = -MM*p;

		// Hardening
		double HH0;           // Hardening moduli
		double hp;            // Plastic coeficient
		double v    = Ivs[1]; // specific volume
		double chi  = (_lam - _kap)/v;
		double trr  = r(0)+r(1)+r(2);
		       HH0 = z0*trr/chi;
		       hp  = -y0*HH0;

		// Elastic stiffness
		Tensor4 De;
		double K = p*v/_kap; // Bulk modulus
		Tensors::AddScaled (2.0*_G,Tensors::Psd, K,Tensors::IdyI, De);

		// Phi and B0
		double Phi = Tensors::Reduce(r,De,r) + hp;
		B.Resize(1); // B0
		Tensors::DotScaled (HH0/Phi,r,De, B[0]); // B[0] = (HH0/Phi)* (r:De)

		// Elastoplastic/secant stiffness
		if (alp<1.0e-5) // elastoplastic only
		{
			Tensors::GerX (-1.0/Phi,De,r,r,De,De, D); // D <- Dep = (-1/Phi) * (De:r)dy(r:De) + De
		}
		else // secant
		{
			Tensor4 Dep;
			Tensors::GerX (-1.0/Phi,De,r,r,De,De, Dep);  // Dep = (-1/Phi) * (De:r)dy(r:De) + De
			Tensors::AddScaled (alp,De, 1.0-alp,Dep, D); // D <- alp*De + (1-alp)*Dep
			B[0] = (1.0-alp)*B[0];
		}
	}
}

inline double CamClay::_calc_M(double const & t) const
{
	return _Mcs*pow( 2.0*_w/(1.0+_w+(_w-1.0)*t) ,0.25);
}

inline double CamClay::_calc_z0(Tensor2 const & Sig) const
{
	double   p,q,t;
	Tensor2  S;
	Tensors::Stress_p_q_S_t (Sig, p,q,S,t);
	double M = _calc_M (t);
	return p+q*q/(M*M*p);
}

inline double CamClay::_val(char const * Name) const
{
	     if (strcmp(Name,"z0")==0) return _ivs[0];
	else if (strcmp(Name,"v" )==0) return _ivs[1];
	else throw new Fatal("CamClay::_val: There is no value named (%s) in this model", Name);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new CamClay model
Model * CamClayMaker()
{
	return new CamClay();
}

// Register CamClay model into ModelFactory array map
int CamClayRegister()
{
	ModelFactory["CamClay"] = CamClayMaker;
	return 0;
}

// Execute the autoregistration
int __CamClay_dummy_int = CamClayRegister();


#endif // MECHSYS_CAMCLAY_H
