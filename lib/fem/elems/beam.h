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

#ifndef MECHSYS_FEM_BEAM_H
#define MECHSYS_FEM_BEAM_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/lin2.h"
#include "util/exception.h"

namespace FEM
{

class Beam : public Lin2, public EquilibElem
{
public:
	// Constants
	static char const * NAME;

	// Constructor
	Beam () : _E(-1), _A(-1), _Izz(-1) {}

	// Derived methods
	char const * Name() const { return NAME; };

	// Methods
	bool   CheckModel   () const;
	void   SetModel     (char const * ModelName, char const * Prms, char const * Inis);
	void   UpdateState  (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
	void   BackupState  ();
	void   RestoreState ();
	void   CalcDepVars  () const;
	double Val          (int iNodeLocal, char const * Name) const;
	double Val          (char const * Name) const;
	void   Order1Matrix (size_t Index, LinAlg::Matrix<double> & Ke) const; ///< Stiffness
	void   B_Matrix     (LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const;

private:
	// Private methods
	int  _geom () const { return 1; }    ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	bool _beam () const { return true; } ///< This is a beam element

	// Private methods
	void _set_ndim                    (int nDim); ///< Set space dimension
	void _calc_initial_internal_state ();         ///< Calculate initial internal state

	// Data
	double  _E;   ///< Young modulus
	double  _A;   ///< Cross-sectional area
	double  _Izz; ///< Cross-sectional inertia

}; // class Beam

// Beam constants
char const * Beam::NAME = "Beam";


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline bool Beam::CheckModel() const
{
	if (_E<0.0 || _A<0.0 || _Izz<0.0) return false;
	return true;
}

inline void Beam::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// Check _ndim
	if (_ndim<1) throw new Fatal("Beam::SetModel: The space dimension (SetDim) must be set before calling this method");
	if (CheckConnect()==false) throw new Fatal("Beam::SetModel: Connectivity is not correct. Connectivity MUST be set before calling this method");

	/* "E=1 A=1 Izz=1" */
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Set
	for (size_t i=0; i<names.Size(); ++i)
	{
		     if (names[i]=="E")    _E   = values[i];
		else if (names[i]=="A")    _A   = values[i];
		else if (names[i]=="Izz" ) _Izz = values[i];
		else throw new Fatal("Beam::SetModel: Parameter name (%s) is invalid",names[i].CStr());
	}
}

inline void Beam::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
{
	// Allocate (local/element) displacements vector
	LinAlg::Vector<double> du(_nd*_n_nodes); // Delta disp. of this element

	// Assemble (local/element) displacements vector
	for (size_t i=0; i<_n_nodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		du(i*_ndim+j) = dUglobal(_connects[i]->DOFVar(UD[_d][j]).EqID);

	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> df(_nd*_n_nodes); // Delta internal force of this element
	df.SetValues(0.0);



	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_n_nodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		dFint(_connects[i]->DOFVar(UD[_d][j]).EqID) += df(i*_ndim+j);
}

inline void Beam::BackupState()
{
}

inline void Beam::RestoreState()
{
}

inline void Beam::CalcDepVars() const
{
}

inline double Beam::Val(int iNodeLocal, char const * Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).NaturalVal;

	throw new Fatal("Beam::Val: This element does not have a Val named %s",Name);
}

inline double Beam::Val(char const * Name) const
{
	throw new Fatal("Beam::Val: Feature not available");
}

inline void Beam::Order1Matrix(size_t Index, LinAlg::Matrix<double> & Ke) const
{
	if (_ndim==2)
	{
		double dx = _connects[1]->X()-_connects[0]->X();
		double dy = _connects[1]->Y()-_connects[0]->Y();
		double LL = dx*dx+dy*dy;
		double L  = sqrt(LL);
		double c  = dx/L;
		double s  = dy/L;
		double c1 = _E*(_A*c*c+12.0*_Izz*s*s/LL)/L;
		double c2 = _E*((_A-12.0*_Izz/LL)*c*s)/L;
		double c3 = _E*(6.0*_Izz*s/L)/L;
		double c4 = _E*(_A*s*s+12.0*_Izz*c*c/LL)/L;
		double c5 = _E*(6.0*_Izz*c/L)/L;
		double c6 = _E*(4.0*_Izz)/L;
		double c7 = _E*(2.0*_Izz)/L;
		Ke.Resize(_nd*_n_nodes, _nd*_n_nodes);
		Ke =  c1,  c2, -c3, -c1, -c2, -c3,
		      c2,  c4,  c5, -c2, -c4,  c5,
		     -c3,  c5,  c6,  c3, -c5,  c7,
		     -c1, -c2,  c3,  c1,  c2,  c3,
		     -c2, -c4, -c5,  c2,  c4, -c5,
		     -c3,  c5,  c7,  c3, -c5,  c6;
	}
	else throw new Fatal("Beam::Order1Matrix: Feature no available for nDim==%d",_ndim);
}

inline void Beam::B_Matrix(LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const
{
	throw new Fatal("Beam::B_Matrix: Feature not available");
}


/* private */

inline void Beam::_set_ndim(int nDim)
{
	if (nDim<1) throw new Fatal("Beam::_set_ndim: For this element, nDim must be greater than or equal to 1 (%d is invalid)",nDim);
	_ndim = nDim;
	_d    = _ndim-1;
	_nd   = EquilibElem::NDB[_d];
}

inline void Beam::_calc_initial_internal_state()
{
	throw new Fatal("Beam::_calc_initial_internal_state: Feature not available");
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Beam element
Element * BeamMaker()
{
	return new Beam();
}

// Register a Beam element into ElementFactory array map
int BeamRegister()
{
	ElementFactory[Beam::NAME] = BeamMaker;
	return 0;
}

// Execute the autoregistration
int __Beam_dummy_int  = BeamRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_BEAM_H
