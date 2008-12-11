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

#ifndef MECHSYS_FEM_ROD_H
#define MECHSYS_FEM_ROD_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/probelem.h"
#include "util/exception.h"
#include "fem/elems/vtkCellType.h"

namespace FEM
{

class Rod: public EquilibElem
{
public:
	// Constructor
	Rod () : _gam(0.0), _E(-1), _A(-1) { };

	// Destructor
	~Rod () { }

	// Methods related to PROBLEM
	void    AddVolForces (); 
	void    ClearDisp    ();
	void    CalcDeps     () const;
	Str_t   ModelName    () const { return "__no_model__"; }
	double  Val          (int iNod, Str_t Key) const;
	double  Val          (          Str_t Key) const;
	void    SetProps     (Str_t Properties);
	void    SetModel     (Str_t ModelName, Str_t Prms, Str_t Inis);
	void    Update       (double h, Vec_t const & dU, Vec_t & dFint);
	void    Backup       () {}
	void    Restore      () {}
	void    CMatrix      (size_t Index, Mat_t & Ke) const;

	// Derived methods
	bool CheckModel () const;

protected:
	// Data
	double _gam; ///< Specific weigth
	double _E; ///< Young modulus
	double _A; ///< Cross-sectional area

	// Derived methods
	void _initialize ();

	// Private methods
	void _B_mat           (Mat_t const & dN, Mat_t const & J, Mat_t & B) const; ///< Calculate B matrix
	void _initial_state   ();
	
	// Depedent variables (calculated by CalcDepVars)
	mutable double         _L;  ///< Rod length
	mutable Vector<double> _uL; ///< Rod-Local displacements/rotations

	void    _transf_mat (Mat_t & T) const; ///< Calculate transformation matrix
	double  N           (double l)  const;

}; // class Rod


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void Rod::AddVolForces() 
{
	// Verify if element is active
	if (IsActive==false) return;

	// Weight
	double dx = _ge->Conn[1]->X()-_ge->Conn[0]->X();
	double dy = _ge->Conn[1]->Y()-_ge->Conn[0]->Y();
	double L  = sqrt(dx*dx+dy*dy);
	double W  = _A*L*_gam;

	// Set boundary conditions
	if (_ge->NDim==1) throw new Fatal("Rod::ApplyBodyForces: feature not available for NDim==1");
	else if (_ge->NDim==2)
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
	// Element displacements vector
	_uL.Resize(_nd*_ge->NNodes);
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		_uL(i*_nd+j) = _ge->Conn[i]->DOFVar(UD[j]).EssentialVal;

	// Transform to rod-local coordinates
	Mat_t T;
	_transf_mat(T);
	_uL = T * _uL;
}

inline double Rod::Val(int iNodeLocal, char const * Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[j])==0) return _ge->Conn[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[j])==0) return _ge->Conn[iNodeLocal]->DOFVar(Name).NaturalVal;

	if (_uL.Size()<1) throw new Fatal("Rod::Val: Please, call CalcDepVars() before calling this method");
	double l = (iNodeLocal==0 ? 0 : 1.0);
	     if (strcmp(Name,"N" )==0) return N(l);
	else if (strcmp(Name,"Ea")==0) return    (_uL(_nd)-_uL(0))/_L;
	else if (strcmp(Name,"Sa")==0) return _E*(_uL(_nd)-_uL(0))/_L;
	else throw new Fatal("Rod::Val: This element does not have a Val named %s",Name);
}

inline double Rod::Val(char const * Name) const
{
	throw new Fatal("Rod::Val: Feature not available");
}

inline void Rod::SetProps(char const * Properties)
{
	/* "gam=20 */
	LineParser lp(Properties);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Set
	for (size_t i=0; i<names.Size(); ++i)
	{
		 if (names[i]=="gam") _gam = values[i];
	}
}

inline void Rod::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// Check NDim
	if (_ge->NDim<1) throw new Fatal("Rod::SetModel: The space dimension (SetDim) must be set before calling this method");
	if (_ge->CheckConn()==false) throw new Fatal("Rod::SetModel: Connectivity is not correct. Connectivity MUST be set before calling this method");

	/* "E=1 A=1" */
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Set
	for (size_t i=0; i<names.Size(); ++i)
	{
		     if (names[i]=="E") _E = values[i];
		else if (names[i]=="A") _A = values[i];
		else throw new Fatal("Rod::SetModel: Parameter name (%s) is invalid",names[i].CStr());
	}
}

inline void Rod::Update(double TimeInc, Vec_t const & dUglobal, Vec_t & dFint)
{
	// Allocate (local/element) displacements vector
	Vec_t du(_nd*_ge->NNodes); // Delta disp. of this element

	// Assemble (local/element) displacements vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		du(i*_nd+j) = dUglobal(_ge->Conn[i]->DOFVar(UD[j]).EqID);

	// Allocate (local/element) internal force vector
	Vec_t df(_nd*_ge->NNodes); // Delta internal force of this element
	df.SetValues(0.0);

	Mat_t Ke;
	CMatrix(0,Ke);
	df = Ke * du;

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		dFint(_ge->Conn[i]->DOFVar(UD[j]).EqID) += df(i*_nd+j);
}

inline void Rod::CMatrix(size_t Index, Mat_t & Ke) const
{
	if (_ge->NDim==2)
	{
		double dx = _ge->Conn[1]->X()-_ge->Conn[0]->X();
		double dy = _ge->Conn[1]->Y()-_ge->Conn[0]->Y();
		double LL = dx*dx+dy*dy;
		      _L  = sqrt(LL);
		double c  = dx/_L;
		double s  = dy/_L;
		double c1 = _E*_A*c*c/_L;
		double c2 = _E*_A*s*c/_L;
		double c3 = _E*_A*s*s/_L;
		Ke.Resize(_nd*_ge->NNodes, _nd*_ge->NNodes);
		Ke =   c1,  c2, -c1, -c2,
		       c2,  c3, -c2, -c3,
		      -c1, -c2,  c1,  c2,
		      -c2, -c3,  c2,  c3;
	}
	else throw new Fatal("Rod::CMatrix: Feature not available for nDim==%d",_ge->NDim);
}

inline void Rod::_B_mat(Mat_t const & derivs, Mat_t const & J, Mat_t & B) const
{
	throw new Fatal("Rod::_B_mat: Feature not available");
}

inline double Rod::N(double l) const
{
	return _E*_A*(_uL(_nd)-_uL(0))/_L;
}

/* private */

inline void Rod::_initialize()
{
	if (_ge->NDim==2)
	{
		_gi = 1;
		_nd = ND_EQUILIB_2D;
		UD  = UD_EQUILIB_2D;
		FD  = FD_EQUILIB_2D;
		_nl = NL_ROD; 
		LB  = LB_ROD;
	}
	else if (_ge->NDim==3)
	{
		_gi = 0;
		_nd = ND_EQUILIB_3D;
		UD  = UD_EQUILIB_3D;
		FD  = FD_EQUILIB_3D;
		_nl = NL_ROD; 
		LB  = LB_ROD;
	}
	
}

inline void Rod::_initial_state()
{
	throw new Fatal("Rod::_initial_state: Feature not available");
}

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

inline bool Rod::CheckModel() const
{
	if (_E<0.0 || _A<0.0) return false;
	return true;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
ProbElem * RodMaker()
{
	return new Rod();
}

// Register element
int RodRegister()
{
	ProbElemFactory["Rod"] = RodMaker;
	return 0;
}

// Call register
int __Rod_dummy_int  = RodRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_ROD_H
