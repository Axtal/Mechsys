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

#ifndef MECHSYS_FEM_ROD3_H
#define MECHSYS_FEM_ROD3_H

// MechSys
#include "fem/equilibelem.h"
#include "util/string.h"
#include "util/util.h"
#include "util/exception.h"

namespace FEM
{

class Rod3: public EquilibElem
{
public:
	// Constants
	static const size_t NL_ROD3;         ///< Number of labels
	static const char   LB_ROD3  [3][4]; ///< Name of labels
	static const char   ROD3_PROP[1][8]; ///< Properties

	// Methods related to PROBLEM
	int         InitCtes     (int nDim);
	void        AddVolForces ();
	int         NProps       () const { return 1; } ///< "gam"
	ProName_t * Props        () const { return ROD3_PROP; }
	void        ClearDisp    ();
	void        CalcDeps     () const;
	double      Val          (int iNod, Str_t Key) const;
	double      Val          (          Str_t Key) const;
	void        Update       (double h, Vec_t const & dU, Vec_t & dFint);
	void        OutInfo      (std::ostream & os) const;
	void        CMatrix      (size_t Idx, Mat_t & M) const;

protected:
	// Depedent variables (calculated by CalcDeps)
	mutable Vector<double> _eps_np; ///< Strain at nodal points
	mutable Vector<double> _eps_ip; ///< Strain at integration points

	// Private methods
	void _excavate   () { throw new Fatal("Rod3::_excavate: Method not available"); }
	void _initialize (Str_t Inis) {}                                       ///< Initialize the element
	void _B_mat      (Mat_t const & dN, Mat_t const & J, Mat_t & B) const; ///< Calculate B matrix

}; // class Rod3

const size_t Rod3::NL_ROD3         = 3;
const char   Rod3::LB_ROD3  [3][4] = {"Ea", "Sa", "N"};
const char   Rod3::ROD3_PROP[1][8] = {"gam"};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline int Rod3::InitCtes(int nDim)
{
	if (nDim==2)
	{
		_gi = 1;
		_nd = ND_EQUILIB_2D;
		UD  = UD_EQUILIB_2D;
		FD  = FD_EQUILIB_2D;
		_nl = NL_ROD3; 
		LB  = LB_ROD3;
	}
	else if (nDim==3)
	{
		_gi = 0;
		_nd = ND_EQUILIB_3D;
		UD  = UD_EQUILIB_3D;
		FD  = FD_EQUILIB_3D;
		_nl = NL_ROD3; 
		LB  = LB_ROD3;
	}
	else throw new Fatal("Rod3::InitCtes: nDim==%d is invalid",nDim);

	// Return geometry index
	return _gi;
}

inline void Rod3::AddVolForces() 
{
	// Verify if element is active
	if (IsActive==false) return;

	// Weight
	double dx = _ge->Conn[1]->X()-_ge->Conn[0]->X();
	double dy = _ge->Conn[1]->Y()-_ge->Conn[0]->Y();
	double L  = sqrt(dx*dx+dy*dy);
	double W  = _mdl->Prm("A")*L*Prop("gam");

	// Set boundary conditions
	if (_ge->NDim==2)
	{
		_ge->Conn[0]->Bry("fy", -W/2.0);
		_ge->Conn[1]->Bry("fy", -W/2.0);
	}
	else if (_ge->NDim==3)
	{
		_ge->Conn[0]->Bry("fz", -W/2.0);
		_ge->Conn[1]->Bry("fz", -W/2.0);
	}
}

inline void Rod3::ClearDisp()
{
	if (IsActive==false) return;

	// Clear displacements
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		_ge->Conn[i]->DOFVar(UD[j]).EssentialVal = 0.0;
}

inline void Rod3::CalcDeps() const
{
	if (IsActive==false) throw new Fatal("Rod3::CalcDeps: This element is inactive");

	// Assemble (local/element) displacements vector
	Vec_t u(_nd*_ge->NNodes);
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		u(i*_nd+j) = _ge->Conn[i]->Val(UD[j]);

	// Calculate strain
	_eps_ip.Resize    (_ge->NIPs);
	_eps_ip.SetValues (0.0);
	Vec_t eps;
	Mat_t dN,J,B;
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_ge->Derivs   (_ge->IPs[i].r, 0, 0, dN);
		_ge->Jacobian (dN, J);
		_B_mat        (dN, J, B);
		eps = B*u;
		_eps_ip(i) = eps(0);
	}
	_ge->Extrap(_eps_ip, _eps_np);
}

inline double Rod3::Val(int iNod, Str_t Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[j])==0) return _ge->Conn[iNod]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[j])==0) return _ge->Conn[iNod]->DOFVar(Name).NaturalVal;

	if (_eps_np.Size()<1) throw new Fatal("Rod3::Val: Please, call CalcDeps() before calling this method");
	     if (strcmp(Name,"Ea")==0) return _eps_np(iNod);
	else if (strcmp(Name,"Sa")==0) return _eps_np(iNod)*_mdl->Prm("E");
	else if (strcmp(Name,"N" )==0) return _eps_np(iNod)*_mdl->Prm("E")*_mdl->Prm("A");
	else throw new Fatal("Rod3::Val: This element does not have a Val named < %s >",Name);
}

inline double Rod3::Val(char const * Name) const
{
	throw new Fatal("Rod3::Val: Value at CG: Method not available yet");
}

inline void Rod3::Update(double h, Vec_t const & dU, Vec_t & dFint)
{
	// Allocate (local/element) displacements vector
	Vec_t du(_nd*_ge->NNodes); // Delta disp. of this element

	// Assemble (local/element) displacements vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		du(i*_nd+j) = dU(_ge->Conn[i]->DOFVar(UD[j]).EqID);

	// Allocate (local/element) internal force vector
	Vec_t df(_nd*_ge->NNodes); // Delta internal force of this element
	df.SetValues (0.0);

	Mat_t Ke;
	CMatrix (0,Ke);
	df = Ke * du;

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		dFint(_ge->Conn[i]->DOFVar(UD[j]).EqID) += df(i*_nd+j);
}

inline void Rod3::OutInfo(std::ostream & os) const
{
	CalcDeps();
	for (size_t i=0; i<_ge->NIPs; i++)
		os << "IP # " << i << ": Ea=" << Util::_8s << _eps_ip(i) << ", Sa=" << Util::_8s << _eps_ip(i)*_mdl->Prm("E") << "  ";
}

inline void Rod3::CMatrix(size_t Idx, Mat_t & Ke) const
{
	/* Stiffness:
	   ==========
	                 /    T
	        [Ke]  =  | [B]  * [D] * [B]  * dV
	                 /
	*/

	double E = _mdl->Prm("E");
	double A = _mdl->Prm("A");
	Ke.Resize    (_nd*_ge->NNodes, _nd*_ge->NNodes);
	Ke.SetValues (0.0);
	Mat_t dN,J,B;
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_ge->Derivs       (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, dN);
		_ge->Jacobian     (dN, J);
		_B_mat            (dN, J, B);
		Ke += trn(B)*E*A*B*det(J)*_ge->IPs[i].w;
	}
}


/* private */

inline void Rod3::_B_mat(Mat_t const & dN, Mat_t const & J, Mat_t & B) const
{
	int dim = _ge->NDim;
	Mat_t dC(dim, dim*_ge->NNodes); // cartesian derivatives
	double det_J = det(J);

	dC.SetValues(0.0);
	for(size_t i=0; i<_ge->NNodes; i++)
	for(int    j=0; j<dim;         j++)
		dC(j, i*dim+j) = dN(0,i);

	// B matrix for a rod element (2D and 3D)
	B = 1.0/(det_J*det_J)*J*dC;
}

///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
ProbElem * Rod3Maker() { return new Rod3(); }

// Register element
int Rod3Register() { ProbElemFactory["Rod3"]=Rod3Maker;  return 0; }

// Call register
int __Rod3_dummy_int  = Rod3Register();


}; // namespace FEM

#endif // MECHSYS_FEM_ROD3_H
