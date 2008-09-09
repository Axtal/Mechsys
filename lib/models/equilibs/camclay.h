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
#include "numerical/brentroot.h"

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
	double _calc_M  (double const & t)    const;                                    ///< Variable CSL slope
	double _calc_f  (Tensor2 const & Sig, IntVals const & Ivs) const;               ///< Yield function
	double _calc_z0 (Tensor2 const & Sig) const;                                    ///< Yield surface size for given stress state
	void   _df_dsig (Tensor2 const & Sig, IntVals const & Ivs,  Tensor2 & V) const; ///< Calculate the V=df_dsig derivative for the (Sig,Ivs) state
	double _f_alpha (double Alpha, void const * Sig, void const * DSig, void const * Ivs) const;

	// Derived private methods
	void   _stiff (Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B) const; ///< Tangent or secant stiffness
	double _val   (char const * Name) const;                                                                                                     ///< Return internal values

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

inline double CamClay::_calc_M(double const & t) const
{
	return _Mcs*pow( 2.0*_w/(1.0+_w+(_w-1.0)*t) ,0.25);
}

inline double CamClay::_calc_f(Tensor2 const & Sig, IntVals const & Ivs) const
{
	// Invariants
	double   p,q,t;
	Tensor2  S;
	Tensors::Stress_p_q_S_t (Sig, p,q,S,t);

	// Yield function
	double z0 = Ivs[0];
	double M  = _calc_M (t);
	return q*q-M*M*p*(z0-p);
}

inline double CamClay::_calc_z0(Tensor2 const & Sig) const
{
	// Invariants
	double   p,q,t;
	Tensor2  S;
	Tensors::Stress_p_q_S_t (Sig, p,q,S,t);

	// z0
	double M = _calc_M (t);
	return p+q*q/(M*M*p);
}

inline void CamClay::_df_dsig(Tensor2 const & Sig, IntVals const & Ivs,  Tensor2 & V) const
{
	// Constants
	const double MIN_Q = 1.0e-5;

	// Invariants
	double   p,q,t;
	Tensor2  S,dqds,dtds;
	Tensors::Stress_p_q_S_t (Sig, p,q,S,t,dqds,dtds);

	// Gradients
	double z0 = Ivs[0];
	double M  = _calc_M (t);
	double m  = M*M*(2.0*p-z0)/3.0;
	if (q>MIN_Q)
	{
		double n    = 2.0*M*p*(p-z0);
		double dMdt = 0.25*M*(1.0-_w)/(1.0+_w-(1.0-_w)*t);
		Tensor2 dMds; dMds = dMdt*dtds;
		V = 3.0*S + m*Tensors::I + n*dMds;
	}
	else V = m*Tensors::I;
}

inline double CamClay::_f_alpha(double Alpha, void const * Sig, void const * DSig, void const * Ivs) const
{
	Tensor2 const &  sig = (*static_cast<Tensor2 const *>( Sig));
	Tensor2 const & dsig = (*static_cast<Tensor2 const *>(DSig));
	Tensor2 siga;
	siga = sig + Alpha*dsig;
	return _calc_f (siga, (*static_cast<IntVals const *>(Ivs)));
}

inline void CamClay::_stiff(Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B) const
{
	// B tensors: dzk = Bk:DEps
	B.Resize    (2); // 2 internal values
	B.SetValues (0.0);

	// Elastic tangent stiffness
	double v = Ivs[1];                     // specific volume
	double p = (Sig(0)+Sig(1)+Sig(2))/3.0; // mean stress
	double K = p*v/_kap;                   // bulk modulus
	Tensors::AddScaled (2.0*_G,Tensors::Psd, K,Tensors::IdyI, D);

	// Elastic trial state
	Tensor2 dsig_tr;                // trial stress increment
	Tensor2  sig_tr;                // elastic trial state
	Tensors::Dot (D,DEps, dsig_tr); // dsig_tr = De:DEps
	sig_tr = Sig + dsig_tr;

	// Initial and trial yield values
	double f_ini = _calc_f (Sig,    Ivs); // initial yield function value
	double f_tr  = _calc_f (sig_tr, Ivs); // trial yield function value

	// Elastic/Elastoplastic ?
	bool   intersect = false; // elastic/elastoplastic intersection ?
	double alpha     = 0.0;   // intersection position
	if (f_ini<0.0) // stress state inside yield surface
	{
		if (f_tr>0.0) intersect = true; // going outside => elastic/elastoplastic intersection
		else return;                    // still inside => pure elastic
	}
	else // on the yield surface or slightly outside
	{
		// Loading/unloading?
		Tensor2 V;  _df_dsig (Sig,Ivs, V);          // V = df_dsig derivative
		double num_dL = Tensors::Reduce (V,D,DEps); // num_dL = V:De:DEps
		if (num_dL<=0.0) return;                    // going inside (or tangent) => pure elastic
	}

	// Set current state on yield surface (after intersection)
	Tensor2 sig;  sig  = Sig;
	Tensor2 eps;  eps  = Eps;
	Tensor2 deps; deps = DEps;
	if (intersect)
	{
		// Find intersection
		Numerical::BrentRoot<CamClay> br(this, &CamClay::_f_alpha);
		alpha = br.Solve (0.0, 1.0, &Sig, &dsig_tr, &Ivs);

		// Update stress state up to the intersection
		Tensor2 deps_elastic;  deps_elastic = alpha*DEps;
		Tensor2 dsig_elastic;
		Tensors::Dot (D,deps_elastic, dsig_elastic); // dsig_elastic = De:deps_elastic
		sig += dsig_elastic;
		eps += deps_elastic;
		deps = (1.0-alpha)*DEps; // remaining strain increment (elastoplastic)
	}

	// Invariants
	double   q,t;
	Tensor2  S;
	Tensors::Stress_p_q_S_t (sig, p,q,S,t);

	// Yield surface gradient
	Tensor2 V;             // V = df_dsig derivative
	_df_dsig (sig,Ivs, V); // V <- df_dsig(Sig,Ivs)

	// Hardening
	double M   = _calc_M (t);     // variable CSL slope
	double z0  = Ivs[0];          // z0 (yield surface size)
	double y0  = -M*M*p;          // y0 = df_dz0
	double chi = (_lam - _kap)/v; // Chi
	double trr = V(0)+V(1)+V(2);  // trace of (r=V)
	double HH0 = z0*trr/chi;      // hardening moduli
	double hp  = -y0*HH0;         // plastic coeficient

	// Elastic tangent stiffness
	Tensor4 De;
	K = p*v/_kap;
	Tensors::AddScaled (2.0*_G,Tensors::Psd, K,Tensors::IdyI, De);

	// Elastoplastic tangent stiffness
	Tensor4 Dep;
	double Phi = Tensors::Reduce(V,De,V) + hp;       // Phi = V:De:V + hp
	Tensors::GerX      (-1.0/Phi,De,V,V,De,De, Dep); // Dep = (-1/Phi) * (De:V)dy(V:De) + De
	Tensors::DotScaled (HH0/Phi,V,De, B[0]);         // B[0] = (HH0/Phi)*(V:De)

	// Secant stiffness
	if (alpha>0.0)
	{
		Tensors::AddScaled (alpha,D, 1.0-alpha,Dep, D); // D <- alp*De + (1-alp)*Dep
		B[0] = (1.0-alpha)*B[0];                        // B <-          (1-alp)*B
	}
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
