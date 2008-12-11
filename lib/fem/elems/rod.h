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
#include "util/exception.h"
#include "fem/elems/vtkCellType.h"

namespace FEM
{

class Rod : public EquilibElem
{
public:
	// Constructor
	Rod () : _gam(0.0), _E(-1), _A(-1) { NNodes=2; Conn.Resize(NNodes); Conn.SetValues(NULL); }

	// Derived methods
	char const * Name() const { return "Rod"; }

	// Derived methods
	bool   CheckModel   () const;
	void   SetModel     (char const * ModelName, char const * Prms, char const * Inis);
	void   SetProps     (char const * Properties);
	void   UpdateState  (double TimeInc, Vec_t const & dUglobal, Vec_t & dFint);
	void   ApplyBodyForces ();
	void   CalcDepVars  () const;
	double Val          (int iNodeLocal, char const * Name) const;
	double Val          (char const * Name) const;
	void   Order1Matrix (size_t Index, Mat_t & Ke) const; ///< Stiffness
	void   B_Matrix     (Mat_t const & derivs, Mat_t const & J, Mat_t & B) const;
	int    VTKCellType  () const { return VTK_LINE; }
	void   VTKConn   (String & Nodes) const { Nodes.Printf("%d %d",Conn[0]->GetID(),Conn[1]->GetID()); }

	// Methods
	double N(double l) const; ///< Axial force (0 < l < 1) (Must be used after CalcDepVars())

private:
	// Data
	double _gam; ///< Specific weigth
	double _E; ///< Young modulus
	double _A; ///< Cross-sectional area

	// Depedent variables (calculated by CalcDepVars)
	mutable double         _L;  ///< Rod length
	mutable Vector<double> _uL; ///< Rod-Local displacements/rotations

	// Private methods
	int  _geom                        () const { return 1; }              ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void _initialize                  ();                                 ///< Initialize the element
	void _calc_initial_internal_state ();                                 ///< Calculate initial internal state
	void _transf_mat                  (Mat_t & T) const; ///< Calculate transformation matrix

}; // class Rod


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline bool Rod::CheckModel() const
{
	if (_E<0.0 || _A<0.0) return false;
	return true;
}

inline void Rod::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// Check _ndim
	if (_ndim<1) throw new Fatal("Rod::SetModel: The space dimension (SetDim) must be set before calling this method");
	if (CheckConnect()==false) throw new Fatal("Rod::SetModel: Connectivity is not correct. Connectivity MUST be set before calling this method");

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

inline void Rod::UpdateState(double TimeInc, Vec_t const & dUglobal, Vec_t & dFint)
{
	// Allocate (local/element) displacements vector
	Vec_t du(_nd*NNodes); // Delta disp. of this element

	// Assemble (local/element) displacements vector
	for (size_t i=0; i<NNodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		du(i*_nd+j) = dUglobal(Conn[i]->DOFVar(UD[_d][j]).EqID);

	// Allocate (local/element) internal force vector
	Vec_t df(_nd*NNodes); // Delta internal force of this element
	df.SetValues(0.0);

	Mat_t Ke;
	Order1Matrix(0,Ke);
	df = Ke * du;

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<NNodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		dFint(Conn[i]->DOFVar(UD[_d][j]).EqID) += df(i*_nd+j);
}

inline void Rod::ApplyBodyForces() 
{
	// Verify if element is active
	if (_is_active==false) return;

	// Weight
	double dx = Conn[1]->X()-Conn[0]->X();
	double dy = Conn[1]->Y()-Conn[0]->Y();
	double L  = sqrt(dx*dx+dy*dy);
	double W  = _A*L*_gam;

	// Set boundary conditions
	if (_ndim==1) throw new Fatal("Rod::ApplyBodyForces: feature not available for NDim==1");
	else if (_ndim==2)
	{
		Conn[0]->Bry("fy", -W/2.0);
		Conn[1]->Bry("fy", -W/2.0);
	}
	else if (_ndim==3)
	{
		Conn[0]->Bry("fz", -W/2.0);
		Conn[1]->Bry("fz", -W/2.0);
	}
}

inline void Rod::CalcDepVars() const
{
	// Element displacements vector
	_uL.Resize(_nd*NNodes);
	for (size_t i=0; i<NNodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		_uL(i*_nd+j) = Conn[i]->DOFVar(UD[_d][j]).EssentialVal;

	// Transform to rod-local coordinates
	Mat_t T;
	_transf_mat(T);
	_uL = T * _uL;
}

inline double Rod::Val(int iNodeLocal, char const * Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[_d][j])==0) return Conn[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[_d][j])==0) return Conn[iNodeLocal]->DOFVar(Name).NaturalVal;

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

inline void Rod::Order1Matrix(size_t Index, Mat_t & Ke) const
{
	if (_ndim==2)
	{
		double dx = Conn[1]->X()-Conn[0]->X();
		double dy = Conn[1]->Y()-Conn[0]->Y();
		double LL = dx*dx+dy*dy;
		      _L  = sqrt(LL);
		double c  = dx/_L;
		double s  = dy/_L;
		double c1 = _E*_A*c*c/_L;
		double c2 = _E*_A*s*c/_L;
		double c3 = _E*_A*s*s/_L;
		Ke.Resize(_nd*NNodes, _nd*NNodes);
		Ke =   c1,  c2, -c1, -c2,
		       c2,  c3, -c2, -c3,
		      -c1, -c2,  c1,  c2,
		      -c2, -c3,  c2,  c3;
	}
	else throw new Fatal("Rod::Order1Matrix: Feature not available for nDim==%d",_ndim);
}

inline void Rod::B_Matrix(Mat_t const & derivs, Mat_t const & J, Mat_t & B) const
{
	throw new Fatal("Rod::B_Matrix: Feature not available");
}

inline double Rod::N(double l) const
{
	return _E*_A*(_uL(_nd)-_uL(0))/_L;
}


/* private */

inline void Rod::_initialize()
{
	if (_ndim<1) throw new Fatal("Rod::_initialize: For this element, _ndim must be greater than or equal to 1 (%d is invalid)",_ndim);
	_d  = _ndim-1;
	_nd = EquilibElem::ND[_d];
	_nl = EquilibElem::NL[_geom()-1];
}

inline void Rod::_calc_initial_internal_state()
{
	throw new Fatal("Rod::_calc_initial_internal_state: Feature not available");
}

inline void Rod::_transf_mat(Mat_t & T) const
{
	// Transformation matrix
	if (_ndim==2)
	{
		double dx = Conn[1]->X()-Conn[0]->X();
		double dy = Conn[1]->Y()-Conn[0]->Y();
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
	else throw new Fatal("Rod::_transf_mat: Feature not available for nDim==%d",_ndim);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Rod element
Element * RodMaker()
{
	return new Rod();
}

// Register a Rod element into ElementFactory array map
int RodRegister()
{
	ElementFactory["Rod"] = RodMaker;
	return 0;
}

// Execute the autoregistration
int __Rod_dummy_int  = RodRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_ROD_H
