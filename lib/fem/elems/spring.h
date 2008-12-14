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

#ifndef MECHSYS_FEM_SPRING_H
#define MECHSYS_FEM_SPRING_H

// MechSys
#include "fem/equilibelem.h"
#include "util/string.h"
#include "util/util.h"
#include "util/exception.h"

namespace FEM
{

class Spring: public EquilibElem
{
public:
	// Constants
	static const size_t NL_SPRING;         ///< Number of labels
	static const char   LB_SPRING  [1][4]; ///< Name of labels
	static const char   SPRING_PROP[1][8]; ///< Properties

	// Methods related to PROBLEM
	int         InitCtes     (int nDim);
	int         NProps       () const { return 0; }
	ProName_t * Props        () const { return SPRING_PROP; }
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
	mutable Vec_t _uL; ///< Spring-Local displacements/rotations

	// Private methods
	void _excavate   () { throw new Fatal("Spring::_excavate: Method not available"); }
	void _initialize (Str_t Inis) {}    ///< Initialize the element
	void _transf_mat (Mat_t & T) const; ///< Calculate transformation matrix

}; // class Spring

const size_t Spring::NL_SPRING         = 1;
const char   Spring::LB_SPRING  [1][4] = {"N"};
const char   Spring::SPRING_PROP[1][8] = {""};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline int Spring::InitCtes(int nDim)
{
	if (nDim==2)
	{
		_gi = 1;
		_nd = ND_EQUILIB_2D;
		UD  = UD_EQUILIB_2D;
		FD  = FD_EQUILIB_2D;
		_nl = NL_SPRING; 
		LB  = LB_SPRING;
	}
	else if (nDim==3)
	{
		_gi = 0;
		_nd = ND_EQUILIB_3D;
		UD  = UD_EQUILIB_3D;
		FD  = FD_EQUILIB_3D;
		_nl = NL_SPRING; 
		LB  = LB_SPRING;
	}
	else throw new Fatal("Spring::InitCtes: nDim==%d is invalid",nDim);

	// Return geometry index
	return _gi;
}

inline void Spring::ClearDisp()
{
	if (IsActive==false) return;

	// Clear displacements
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		_ge->Conn[i]->DOFVar(UD[j]).EssentialVal = 0.0;
}

inline void Spring::CalcDeps() const
{
	if (IsActive==false) throw new Fatal("Spring::CalcDeps: This element is inactive");

	// Element displacements vector
	_uL.Resize(_nd*_ge->NNodes);
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		_uL(i*_nd+j) = _ge->Conn[i]->Val(UD[j]);

	// Transform to local coordinates
	Mat_t T;
	_transf_mat (T);
	_uL = T * _uL;
}

inline double Spring::Val(int iNod, Str_t Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[j])==0) return _ge->Conn[iNod]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[j])==0) return _ge->Conn[iNod]->DOFVar(Name).NaturalVal;

	if (_uL.Size()<1) throw new Fatal("Spring::Val: Please, call CalcDeps() before calling this method");
	double l = (iNod==0 ? 0 : 1.0);
	if (strcmp(Name,"N" )==0) return N(l);
	else throw new Fatal("Spring::Val: This element does not have a Val named < %s >",Name);
}

inline double Spring::Val(char const * Name) const
{
	if (_uL.Size()<1) throw new Fatal("Beam::Val: Please, call CalcDeps() before calling this method");
	double l = 0.5;
	if (strcmp(Name,"N" )==0) return N(l);
	else throw new Fatal("Spring::Val: This element does not have a Val named < %s >",Name);
}

inline void Spring::Update(double h, Vec_t const & dU, Vec_t & dFint)
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

inline void Spring::CMatrix(size_t Idx, Mat_t & Ke) const
{
	/*          T   T
	     K = [T0]*[B]*k0*[B]*[T0]*Area
	*/

	Mat_t B(1,2); B = 1, -1;
	Mat_t T;
	_transf_mat (T);
	Ke = trn(T)*trn(B)*_mdl->Prm("ks")*B*T;
}

inline double Spring::N(double l) const
{
	return _mdl->Prm("ks")*(_uL(_nd)-_uL(0));
}


/* private */

inline void Spring::_transf_mat(Mat_t & T) const
{
	double dx = _ge->Conn[1]->X() - _ge->Conn[0]->X();
	double dy = _ge->Conn[1]->Y() - _ge->Conn[0]->Y();
	double dz = _ge->Conn[1]->Z() - _ge->Conn[0]->Z();
	double LL = dx*dx + dy*dy + dz*dz;
	double L  = sqrt(LL);
	double l  = dx/L;
	double m  = dy/L;
	double n  = dz/L;
	T.Resize(2, _ge->NDim*2);
	if (_ge->NDim==2)
	{
		T = l, m, 0, 0,
		    0, 0, l, m;
	}
	else
	{
		T = l, m, n, 0, 0, 0,
		    0, 0, 0, l, m, n;
	}
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
ProbElem * SpringMaker() { return new Spring(); }

// Register element
int SpringRegister() { ProbElemFactory["Spring"]=SpringMaker;  return 0; }

// Call register
int __Spring_dummy_int  = SpringRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_SPRING_H
