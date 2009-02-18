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

using Tensors::AddScaled;
using Tensors::IIsym;
using Tensors::IdyI;

class CamClay : public EquilibModel
{
public:
	// Constants
	static const char CAMCLAY_PN[4][8]; ///< Parameters

	// Destructor
	virtual ~CamClay () {}

	// Derived methods
	int         NPrms   () const { return 4;          }
	PrmName_t * Prms    () const { return CAMCLAY_PN; }
	Str_t       Name    () const { return "CamClay";  }
	void        InitIVS (Ini_t const & Ini, Tensor2 const & Sig, Tensor2 const & Eps, IntVals & Ivs) const;

private:
	// Private data
	double _Mcs;
	double _w;

	// Derived private methods
	void _initialize ();
	void _stiff      (Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B, bool First) const;

	// Private methods
	double _calc_M  (double const & t)    const;                                    ///< Variable CSL slope
	double _calc_f  (Tensor2 const & Sig, IntVals const & Ivs) const;               ///< Yield function
	double _calc_z0 (Tensor2 const & Sig) const;                                    ///< Yield surface size for given stress state
	void   _df_dsig (Tensor2 const & Sig, IntVals const & Ivs,  Tensor2 & V) const; ///< Calculate the V=df_dsig derivative for the (Sig,Ivs) state
	double _f_alpha (double Alpha, void const * Sig, void const * DSig, void const * Ivs) const;

}; // class CamClay

const char CamClay::CAMCLAY_PN[4][8] = {"lam", "kap", "nu", "phics"};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void CamClay::InitIVS(Ini_t const & Ini, Tensor2 const & Sig, Tensor2 const & Eps, IntVals & Ivs) const
{
	Ini_t::const_iterator it = Ini.find("v");
	if (it==Ini.end()) throw new Fatal("CamClay::InitIVS: Specific volume (v) must be provided with the array of Initial Values");
	Ivs.Resize(2); // z0 (yield surface size), and v (specific volume)
	Ivs[0] = _calc_z0 (Sig);
	Ivs[1] = it->second;
}


/* private */

inline void CamClay::_initialize()
{
	// Parameters
	double lam   = Prm("lam");
	double kap   = Prm("kap");
	double nu    = Prm("nu");
	double phics = Prm("phics");
	
	// Check
	if (lam<=0.0)           throw new Fatal("LinElastic::_initialize: Tag=%d: Compressibility index (lam) must be positive. lam==%f is invalid",_tag,lam);
	if (kap<=0.0)           throw new Fatal("LinElastic::_initialize: Tag=%d: Rebounding index (kap) must be positive. kap==%f is invalid",_tag,kap);
	if (nu<0.0 || nu>0.499) throw new Fatal("LinElastic::_initialize: Tag=%d: Poisson ratio (nu) must be provided (and in the range: 0 < nu < 0.5). nu==%f is invalid",_tag,nu);
	if (phics<=0.0)         throw new Fatal("LinElastic::_initialize: Tag=%d: Friction angle at critical state (phics) must be positive. phics==%f is invalid",_tag,phics);
	
	// Derived constants
	double sinphi = sin(phics*Util::PI/180.0);
	_Mcs = 6.0*sinphi/(3.0-sinphi);
	_w   = pow((3.0-sinphi)/(3.0+sinphi),4.0);
}

inline void CamClay::_stiff(Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B, bool First) const
{
	// Parameters
	double lam = Prm("lam");
	double kap = Prm("kap");
	double nu  = Prm("nu");

	// B tensors: dzk = Bk:DEps
	B.Resize    (2); // 2 internal values
	B.SetValues (0.0);

	// Elastic tangent stiffness
	double v = Ivs[1];                                                    // specific volume
	double E = -(Sig(0)+Sig(1)+Sig(2))*(1.0-2.0*nu)*v/kap;                // Young modulus
	AddScaled (E/(1.0+nu), IIsym, nu*E/((1.0+nu)*(1.0-2.0*nu)), IdyI, D); // Elastic tangent tensor

	// Return if first stiffness
	if (First) return;

	// Elastic trial state
	Tensor2 dsig_tr;                // trial stress increment
	Tensor2  sig_tr;                // elastic trial state
	Tensors::Dot (D,DEps, dsig_tr); // dsig_tr = De:DEps
	sig_tr = Sig + dsig_tr;

	// Initial and trial yield values
	double f_ini = _calc_f (Sig,    Ivs); // initial yield function value
	double f_tr  = _calc_f (sig_tr, Ivs); // trial yield function value

	//std::cout << f_ini << std::endl;

	// Elastic/Elastoplastic ?
	bool   intersect = false; // elastic/elastoplastic intersection ?
	double alpha     = 0.0;   // intersection position
	if (f_ini<0.0) // stress state inside yield surface
	{
		if (f_tr>0.0) intersect = true; // going outside => elastic/elastoplastic intersection
		else
		{
			//std::cout << "pure elastic" << std::endl;
			return;                    // still inside => pure elastic
		}
	}
	else // on the yield surface or slightly outside
	{
		// Loading/unloading?
		Tensor2 V;  _df_dsig (Sig,Ivs, V);          // V = df_dsig derivative
		//std::cout << "df_dsig = " << V << std::endl;
		//std::cout << "D = " << D << std::endl;
		//std::cout << "DEps = " << DEps << std::endl;
		double num_dL = Tensors::Reduce (V,D,DEps); // num_dL = V:De:DEps
		//std::cout << "num_dL = " << num_dL << std::endl;
		if (num_dL<=0.0)
		{
			//std::cout << "going inside pure elastic" << std::endl;
			return;                    // going inside (or tangent) => pure elastic
		}
	}

	//std::cout << "intersect = " << intersect << std::endl;

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
	double   p,q,t;
	Tensor2  S;
	Tensors::Stress_p_q_S_t (sig, p,q,S,t);

	// Yield surface gradient
	Tensor2 V;             // V = df_dsig derivative
	_df_dsig (sig,Ivs, V); // V <- df_dsig(Sig,Ivs)

	// Hardening
	double M   = _calc_M (t);     // variable CSL slope
	double z0  = Ivs[0];          // z0 (yield surface size)
	double y0  = -M*M*p;          // y0 = df_dz0
	double chi = -(lam - kap)/v;  // Chi
	double trr = V(0)+V(1)+V(2);  // trace of (r=V)
	double HH0 = z0*trr/chi;      // hardening moduli
	double hp  = -y0*HH0;         // plastic coeficient

	// Elastic tangent stiffness
	Tensor4 De;
	E = -3.0*p*(1.0-2.0*nu)*v/kap;
	AddScaled (E/(1.0+nu), IIsym, nu*E/((1.0+nu)*(1.0-2.0*nu)), IdyI, De);

	// Elastoplastic tangent stiffness
	Tensor4 Dep;
	double Phi = Tensors::Reduce(V,De,V) + hp;       // Phi = V:De:V + hp
	Tensors::GerX      (-1.0/Phi,De,V,V,De,De, Dep); // Dep = (-1/Phi) * (De:V)dy(V:De) + De
	Tensors::DotScaled (HH0/Phi,V,De, B[0]);         // B[0] = (HH0/Phi)*(V:De)

	//std::cout << "alpha = " << alpha << ", B[0]=" << B[0] << std::endl;
	//std::cout << "DEps = " << _8s<<DEps(0) << _8s<<DEps(1) << _8s<<DEps(2) << _8s<<DEps(3) << _8s<<DEps(4) << _8s<<DEps(5) << std::endl;
	//std::cout << "B[0] = " << _8s<<B[0](0) << _8s<<B[0](1) << _8s<<B[0](2) << _8s<<B[0](3) << _8s<<B[0](4) << _8s<<B[0](5) << std::endl;
	//std::cout << "B[0]:DEps = " << Tensors::Dot(B[0],DEps) << std::endl;

	// Secant stiffness
	//if (alpha>0.0)
	//{
		//Tensors::AddScaled (alpha,D, 1.0-alpha,Dep, D); // D <- alp*De + (1-alp)*Dep
		//B[0] = (1.0-alpha)*B[0];                        // B <-          (1-alp)*B
	//}
	
	D = Dep;
	//std::cout << D << std::endl;
}

inline double CamClay::_calc_M(double const & t) const
{
	return _Mcs*pow( 2.0*_w/(1.0+_w+(_w-1.0)*(-t)) ,0.25);
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
	std::cout << "t = " << t << ", M = " << M << std::endl;
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
	//std::cout << "t = " << t << ", M = " << M << std::endl;
	if (q>MIN_Q)
	{
		double n    = 2.0*M*p*(p-z0);
		double dMdt = -0.25*M*(1.0-_w)/(1.0+_w-(1.0-_w)*(-t));
		Tensor2 dMds; dMds = dMdt*dtds;
		V = 3.0*S + m*Tensors::I;// + n*dMds;
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


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new model
Model * CamClayMaker() { return new CamClay(); }

// Register model
int CamClayRegister() { ModelFactory["CamClay"]=CamClayMaker;  return 0; }

// Call register
int __CamClay_dummy_int = CamClayRegister();


#endif // MECHSYS_CAMCLAY_H
