/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo, Fernando Alonso                   *
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

#ifndef MECHSYS_FEM_ROD_H
#define MECHSYS_FEM_ROD_H

// MechSys
#include "fem/equilibelem.h"
#include "util/string.h"
#include "util/util.h"
#include "util/exception.h"

namespace FEM
{

class Rod: public EquilibElem
{
public:
	// Constants
	static const size_t NL_ROD;         ///< Number of labels
	static const char   LB_ROD  [3][4]; ///< Name of labels
	static const char   ROD_PROP[1][8]; ///< Properties

	// Methods related to PROBLEM
	int         InitCtes     (int nDim);
	void        AddVolForces ();
	int         NProps       () const { return 1; } ///< "gam"
	ProName_t * Props        () const { return ROD_PROP; }
	void        ClearDisp    ();
	void        CalcDeps     () const;
	double      Val          (int iNod, Str_t Key) const;
	double      Val          (          Str_t Key) const;
	void        Update       (double h, Vec_t const & dU, Vec_t & dFint);
	void        CMatrix      (size_t Idx, Mat_t & M) const;

	// Methods
	double N (double l) const; ///< Axial force (0 < l < 1) (Must be used after CalcDeps())

protected:
	// Depedent variables (calculated by CalcDeps)
	mutable double _L;  ///< Rod length
	mutable Vec_t  _uL; ///< Rod-Local displacements/rotations

	// Private methods
	void _excavate   () { throw new Fatal("Rod::_excavate: Method not available"); }
	void _initialize (Str_t Inis) {}    ///< Initialize the element
	void _transf_mat (Mat_t & T) const; ///< Calculate transformation matrix

}; // class Rod

const size_t Rod::NL_ROD         = 3;
const char   Rod::LB_ROD  [3][4] = {"Ea", "Sa", "N"};
const char   Rod::ROD_PROP[1][8] = {"gam"};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline int Rod::InitCtes(int nDim)
{
	if (nDim==2)
	{
		_gi = 1;
		_nd = ND_EQUILIB_2D;
		UD  = UD_EQUILIB_2D;
		FD  = FD_EQUILIB_2D;
		_nl = NL_ROD; 
		LB  = LB_ROD;
	}
	else if (nDim==3)
	{
		_gi = 0;
		_nd = ND_EQUILIB_3D;
		UD  = UD_EQUILIB_3D;
		FD  = FD_EQUILIB_3D;
		_nl = NL_ROD; 
		LB  = LB_ROD;
	}
	else throw new Fatal("Rod::InitCtes: nDim==%d is invalid",nDim);

	// Return geometry index
	return _gi;
}

inline void Rod::AddVolForces() 
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

inline void Rod::ClearDisp()
{
	if (IsActive==false) return;

	// Clear displacements
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		_ge->Conn[i]->DOFVar(UD[j]).EssentialVal = 0.0;
}

inline void Rod::CalcDeps() const
{
	if (IsActive==false) throw new Fatal("Rod::CalcDeps: This element is inactive");

	// Element displacements vector
	_uL.Resize(_nd*_ge->NNodes);
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		_uL(i*_nd+j) = _ge->Conn[i]->DOFVar(UD[j]).EssentialVal;

	// Transform to rod-local coordinates
	Mat_t T;
	_transf_mat(T);
	_uL = T * _uL;
}

inline double Rod::Val(int iNod, Str_t Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[j])==0) return _ge->Conn[iNod]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[j])==0) return _ge->Conn[iNod]->DOFVar(Name).NaturalVal;

	if (_uL.Size()<1) throw new Fatal("Rod::Val: Please, call CalcDeps() before calling this method");
	double l = (iNod==0 ? 0 : 1.0);
	     if (strcmp(Name,"N" )==0) return N(l);
	else if (strcmp(Name,"Ea")==0) return                (_uL(_nd)-_uL(0))/_L;
	else if (strcmp(Name,"Sa")==0) return _mdl->Prm("E")*(_uL(_nd)-_uL(0))/_L;
	else throw new Fatal("Rod::Val: This element does not have a Val named < %s >",Name);
}

inline double Rod::Val(char const * Name) const
{
	if (_uL.Size()<1) throw new Fatal("Beam::Val: Please, call CalcDeps() before calling this method");
	double l = 0.5;
	     if (strcmp(Name,"N" )==0) return N(l);
	else if (strcmp(Name,"Ea")==0) return                (_uL(_nd)-_uL(0))/_L;
	else if (strcmp(Name,"Sa")==0) return _mdl->Prm("E")*(_uL(_nd)-_uL(0))/_L;
	else throw new Fatal("Rod::Val: This element does not have a Val named < %s >",Name);
}

inline void Rod::Update(double h, Vec_t const & dU, Vec_t & dFint)
{
	// Allocate (local/element) displacements vector
	Vec_t du(_nd*_ge->NNodes); // Delta disp. of this element

	// Assemble (local/element) displacements vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		du(i*_nd+j) = dU(_ge->Conn[i]->DOFVar(UD[j]).EqID);

	// Calculate internal force vector
	Vec_t df;
	Mat_t Ke;
	CMatrix (0,Ke);
	df = Ke * du;

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		dFint(_ge->Conn[i]->DOFVar(UD[j]).EqID) += df(i*_nd+j);
}

inline void Rod::CMatrix(size_t Idx, Mat_t & Ke) const
{
	if (_ge->NDim==2)
	{
		double E  = _mdl->Prm("E");
		double A  = _mdl->Prm("A");
		double dx = _ge->Conn[1]->X()-_ge->Conn[0]->X();
		double dy = _ge->Conn[1]->Y()-_ge->Conn[0]->Y();
		double LL = dx*dx+dy*dy;
		      _L  = sqrt(LL);
		double c  = dx/_L;
		double s  = dy/_L;
		double c1 = E*A*c*c/_L;
		double c2 = E*A*s*c/_L;
		double c3 = E*A*s*s/_L;
		Ke.Resize(_nd*_ge->NNodes, _nd*_ge->NNodes);
		Ke =   c1,  c2, -c1, -c2,
		       c2,  c3, -c2, -c3,
		      -c1, -c2,  c1,  c2,
		      -c2, -c3,  c2,  c3;
	}
	else throw new Fatal("Rod::CMatrix: Feature not available for nDim==%d",_ge->NDim);
}

inline double Rod::N(double l) const
{
	return _mdl->Prm("E")*_mdl->Prm("A")*(_uL(_nd)-_uL(0))/_L;
}


/* private */

inline void Rod::_transf_mat(Mat_t & T) const
{
	// Transformation matrix
	if (_ge->NDim==2)
	{
		double dx = _ge->Conn[1]->X()-_ge->Conn[0]->X();
		double dy = _ge->Conn[1]->Y()-_ge->Conn[0]->Y();
		double LL = dx*dx+dy*dy;
		      _L  = sqrt(LL);
		double c  = dx/_L;
		double s  = dy/_L;
		T.Resize(4,4);
		T =    c,   s, 0.0, 0.0,
		      -s,   c, 0.0, 0.0,
		     0.0, 0.0,   c,   s,
		     0.0, 0.0,  -s,   c;
	}
	else throw new Fatal("Rod::_transf_mat: Feature not available for nDim==%d",_ge->NDim);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
ProbElem * RodMaker() { return new Rod(); }

// Register element
int RodRegister() { ProbElemFactory["Rod"]=RodMaker;  return 0; }

// Call register
int __Rod_dummy_int  = RodRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_ROD_H
