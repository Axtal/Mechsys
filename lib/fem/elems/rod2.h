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

#ifndef MECHSYS_FEM_ROD2_H
#define MECHSYS_FEM_ROD2_H

// MechSys
#include "fem/equilibelem.h"
#include "util/exception.h"
#include "fem/elems/vtkCellType.h"

namespace FEM
{

class Rod2 : public EquilibElem
{
public:
	// Constructor
	Rod2 () : _gam(0.0), _E(-1), _A(-1) { _n_nodes=2; _connects.Resize(_n_nodes); _connects.SetValues(NULL); }

	// Derived methods
	char const * Name() const { return "Rod2"; }

	// Derived methods
	bool   CheckModel   () const;
	void   SetModel     (char const * ModelName, char const * Prms, char const * Inis);
	void   SetProps     (char const * Properties);
	void   UpdateState  (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
	void   ApplyBodyForces ();
	void   CalcDepVars  () const;
	double Val          (int iNodeLocal, char const * Name) const;
	double Val          (char const * Name) const;
	void   Order1Matrix (size_t Index, LinAlg::Matrix<double> & Ke) const; ///< Stiffness
	void   B_Matrix     (LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const;
	int    VTKCellType  () const { return VTK_LINE; }
	void   VTKConnect   (String & Nodes) const { Nodes.Printf("%d %d",_connects[0]->GetID(),_connects[1]->GetID()); }

	// Methods
	double N(double l) const; ///< Axial force (0 < l < 1) (Must be used after CalcDepVars())

private:
	// Data
	double _gam; ///< Specific weigth
	double _E; ///< Young modulus
	double _A; ///< Cross-sectional area

	// Depedent variables (calculated by CalcDepVars)
	mutable double         _L;  ///< Rod2 length
	mutable Vector<double> _uL; ///< Rod2-Local displacements/rotations

	// Private methods
	int  _geom                        () const { return 1; }              ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void _initialize                  ();                                 ///< Initialize the element
	void _calc_initial_internal_state ();                                 ///< Calculate initial internal state
	void _transf_mat                  (LinAlg::Matrix<double> & T) const; ///< Calculate transformation matrix

}; // class Rod2


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline bool Rod2::CheckModel() const
{
	if (_E<0.0 || _A<0.0) return false;
	return true;
}

inline void Rod2::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// Check _ndim
	if (_ndim<1) throw new Fatal("Rod2::SetModel: The space dimension (SetDim) must be set before calling this method");
	if (CheckConnect()==false) throw new Fatal("Rod2::SetModel: Connectivity is not correct. Connectivity MUST be set before calling this method");

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
		else throw new Fatal("Rod2::SetModel: Parameter name (%s) is invalid",names[i].CStr());
	}
}

inline void Rod2::SetProps(char const * Properties)
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

inline void Rod2::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
{
	// Allocate (local/element) displacements vector
	LinAlg::Vector<double> du(_nd*_n_nodes); // Delta disp. of this element

	// Assemble (local/element) displacements vector
	for (size_t i=0; i<_n_nodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		du(i*_nd+j) = dUglobal(_connects[i]->DOFVar(UD[_d][j]).EqID);

	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> df(_nd*_n_nodes); // Delta internal force of this element
	df.SetValues(0.0);

	LinAlg::Matrix<double> Ke;
	Order1Matrix(0,Ke);
	df = Ke * du;

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_n_nodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		dFint(_connects[i]->DOFVar(UD[_d][j]).EqID) += df(i*_nd+j);
}

inline void Rod2::ApplyBodyForces() 
{
	// Verify if element is active
	if (_is_active==false) return;

	// Weight
	double dx = _connects[1]->X()-_connects[0]->X();
	double dy = _connects[1]->Y()-_connects[0]->Y();
	double L  = sqrt(dx*dx+dy*dy);
	double W  = _A*L*_gam;

	// Set boundary conditions
	if (_ndim==1) throw new Fatal("Rod2::ApplyBodyForces: feature not available for NDim==1");
	else if (_ndim==2)
	{
		_connects[0]->Bry("fy", -W/2.0);
		_connects[1]->Bry("fy", -W/2.0);
	}
	else if (_ndim==3)
	{
		_connects[0]->Bry("fz", -W/2.0);
		_connects[1]->Bry("fz", -W/2.0);
	}
}

inline void Rod2::CalcDepVars() const
{
	// Element displacements vector
	_uL.Resize(_nd*_n_nodes);
	for (size_t i=0; i<_n_nodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		_uL(i*_nd+j) = _connects[i]->DOFVar(UD[_d][j]).EssentialVal;

	// Transform to rod-local coordinates
	LinAlg::Matrix<double> T;
	_transf_mat(T);
	_uL = T * _uL;
}

inline double Rod2::Val(int iNodeLocal, char const * Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).NaturalVal;

	if (_uL.Size()<1) throw new Fatal("Rod2::Val: Please, call CalcDepVars() before calling this method");
	double l = (iNodeLocal==0 ? 0 : 1.0);
	     if (strcmp(Name,"N" )==0) return N(l);
	else if (strcmp(Name,"Ea")==0) return    (_uL(_nd)-_uL(0))/_L;
	else if (strcmp(Name,"Sa")==0) return _E*(_uL(_nd)-_uL(0))/_L;
	else throw new Fatal("Rod2::Val: This element does not have a Val named %s",Name);
}

inline double Rod2::Val(char const * Name) const
{
	throw new Fatal("Rod2::Val: Feature not available");
}

inline void Rod2::Order1Matrix(size_t Index, LinAlg::Matrix<double> & Ke) const
{
	if (_ndim==2)
	{
		double dx = _connects[1]->X()-_connects[0]->X();
		double dy = _connects[1]->Y()-_connects[0]->Y();
		double LL = dx*dx+dy*dy;
		      _L  = sqrt(LL);
		double c  = dx/_L;
		double s  = dy/_L;
		double c1 = _E*_A*c*c/_L;
		double c2 = _E*_A*s*c/_L;
		double c3 = _E*_A*s*s/_L;
		Ke.Resize(_nd*_n_nodes, _nd*_n_nodes);
		Ke =   c1,  c2, -c1, -c2,
		       c2,  c3, -c2, -c3,
		      -c1, -c2,  c1,  c2,
		      -c2, -c3,  c2,  c3;
	}
	else throw new Fatal("Rod2::Order1Matrix: Feature not available for nDim==%d",_ndim);
}

inline void Rod2::B_Matrix(LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const
{
	throw new Fatal("Rod2::B_Matrix: Feature not available");
}

inline double Rod2::N(double l) const
{
	return _E*_A*(_uL(_nd)-_uL(0))/_L;
}


/* private */

inline void Rod2::_initialize()
{
	if (_ndim<1) throw new Fatal("Rod2::_initialize: For this element, _ndim must be greater than or equal to 1 (%d is invalid)",_ndim);
	_d  = _ndim-1;
	_nd = EquilibElem::ND[_d];
	_nl = EquilibElem::NL[_geom()-1];
}

inline void Rod2::_calc_initial_internal_state()
{
	throw new Fatal("Rod2::_calc_initial_internal_state: Feature not available");
}

inline void Rod2::_transf_mat(LinAlg::Matrix<double> & T) const
{
	// Transformation matrix
	if (_ndim==2)
	{
		double dx = _connects[1]->X()-_connects[0]->X();
		double dy = _connects[1]->Y()-_connects[0]->Y();
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
	else throw new Fatal("Rod2::_transf_mat: Feature not available for nDim==%d",_ndim);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Rod2 element
Element * Rod2Maker()
{
	return new Rod2();
}

// Register a Rod2 element into ElementFactory array map
int Rod2Register()
{
	ElementFactory["Rod2"] = Rod2Maker;
	return 0;
}

// Execute the autoregistration
int __Rod2_dummy_int  = Rod2Register();

}; // namespace FEM

#endif // MECHSYS_FEM_ROD2_H