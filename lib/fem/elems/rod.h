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
	Rod () : _E(-1), _A(-1) { _n_nodes=2; _connects.Resize(_n_nodes); _connects.SetValues(NULL); }

	// Derived methods
	char const * Name() const { return "Rod"; }

	// Derived methods
	bool   CheckModel   () const;
	void   SetModel     (char const * ModelName, char const * Prms, char const * Inis);
	void   UpdateState  (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
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
	double _E; ///< Young modulus
	double _A; ///< Cross-sectional area

	// Depedent variables (calculated by CalcDepVars)
	mutable double         _L;  ///< Rod length
	mutable Vector<double> _uL; ///< Rod-Local displacements/rotations

	// Private methods
	int  _geom                        () const { return 1; }              ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void _set_ndim                    (int nDim);                         ///< Set space dimension
	void _calc_initial_internal_state ();                                 ///< Calculate initial internal state
	void _transf_mat                  (LinAlg::Matrix<double> & T) const; ///< Calculate transformation matrix

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

inline void Rod::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
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

inline void Rod::CalcDepVars() const
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

inline double Rod::Val(int iNodeLocal, char const * Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).NaturalVal;

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

inline void Rod::Order1Matrix(size_t Index, LinAlg::Matrix<double> & Ke) const
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
	else throw new Fatal("Rod::Order1Matrix: Feature not available for nDim==%d",_ndim);
}

inline void Rod::B_Matrix(LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const
{
	throw new Fatal("Rod::B_Matrix: Feature not available");
}

inline double Rod::N(double l) const
{
	return _E*_A*(_uL(_nd)-_uL(0))/_L;
}


/* private */

inline void Rod::_set_ndim(int nDim)
{
	if (nDim<1) throw new Fatal("Rod::_set_ndim: For this element, nDim must be greater than or equal to 1 (%d is invalid)",nDim);
	_ndim = nDim;
	_d    = _ndim-1;
	_nd   = EquilibElem::ND[_d];
}

inline void Rod::_calc_initial_internal_state()
{
	throw new Fatal("Rod::_calc_initial_internal_state: Feature not available");
}

inline void Rod::_transf_mat(LinAlg::Matrix<double> & T) const
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
