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
#include "util/exception.h"
#include "fem/elems/vtkCellType.h"

namespace FEM
{

class Spring : public EquilibElem
{
public:
	// Constructor
	Spring() : _ks(-1) { _n_nodes=2; _connects.Resize(_n_nodes); _connects.SetValues(NULL);}
	

	// Derived methods
	char const * Name() const { return "Spring"; }

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
	void   VTKConnect   (String & Nodes) const;
	void   OutInfo(std::ostream & os) const;
	

private:
	// Data
	double _ks;  ///< Spring stiffness

	// Private methods
	int  _geom     () const { return 1; }     ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void _initialize();

	void _calc_initial_internal_state();
	void _mount_T_matrix(Matrix<double> & T) const;

}; // class Spring


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */
inline bool Spring::CheckModel() const
{
	if (_ks<0.0) return false;
	return true;
}

inline void Spring::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// Check _ndim
	if (_ndim<1) throw new Fatal("Spring::SetModel: The space dimension (SetDim) must be set before calling this method");
	if (CheckConnect()==false) throw new Fatal("Spring::SetModel: Connectivity is not correct. Connectivity MUST be set before calling this method");

	/* "ks=1" */
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Set
	for (size_t i=0; i<names.Size(); ++i)
	{
		     if (names[i]=="ks" )  _ks = values[i];
		else throw new Fatal("Spring::SetModel: Parameter name (%s) is invalid",names[i].CStr());
	}
}

inline void Spring::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
{
	// Allocate (local/element) displacements vector
	LinAlg::Vector<double> du(_nd*_n_nodes); // Delta disp. of this element

	// Assemble (local/element) displacements vector
	for (size_t i=0; i<_n_nodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		du(i*_nd+j) = dUglobal(_connects[i]->DOFVar(UD[_d][j]).EqID);

	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> df; // Delta internal force of this element

	LinAlg::Matrix<double> Ke;
	Order1Matrix(0,Ke);
	df = Ke * du;

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_n_nodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		dFint(_connects[i]->DOFVar(UD[_d][j]).EqID) += df(i*_nd+j);
}

inline void Spring::CalcDepVars() const
{
}

inline double Spring::Val(int iNodeLocal, char const * Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).NaturalVal;

	     if (strcmp(Name,"Ea")==0)  return 0.0;
	else if (strcmp(Name,"Sa")==0)  return 0.0;
	else if (strcmp(Name,"N" )==0)
	{
		// Allocate (local/element) displacements vector
		LinAlg::Vector<double> du(_nd*_n_nodes); // Delta disp. of this element
		// Assemble (local/element) displacements vector
		for (size_t i=0; i<_n_nodes; ++i)
		for (int    j=0; j<_nd;      ++j)
		du(i*_nd+j) = _connects[i]->Val(UD[_d][j]);

		Matrix<double> T; _mount_T_matrix(T);
		Vector<double> D; D = T*du;
		return (D(1)-D(0))*_ks;
	}
	else throw new Fatal("Rod3::Val: This element does not have a Val named %s",Name);
}

inline double Spring::Val(char const * Name) const
{
	throw new Fatal("Spring::Val: Feature not available");
}

inline void Spring::Order1Matrix(size_t Index, LinAlg::Matrix<double> & Ke) const
{
	//              T   T                    
	//       K = [T0]*[B]*k0*[B]*[T0]*Area 
	//        

	// Mount B Matrix
	Matrix<double> B(1,2); B = 1, -1;
	// Mount T Matrix
	Matrix<double> T; _mount_T_matrix(T);

	Ke = trn(T)*trn(B)*_ks*B*T;
}

inline void Spring::_mount_T_matrix(Matrix<double> & T) const
{
	double x0 = _connects[0]->X(), y0 = _connects[0]->Y(), z0 = _connects[0]->Z();
	double x1 = _connects[1]->X(), y1 = _connects[1]->Y(), z1 = _connects[1]->Z();
	double L  = sqrt(pow(x0-x1,2)+pow(y0-y1,2)+pow(z0-z1,2));
	double l = (x0-x1)/L;
	double m = (y0-y1)/L;
	double n = (z0-z1)/L;

	// Mount T Matrix
	T.Resize(2, _ndim*2);
	if (_ndim==2)
		T = l, m, 0, 0, 
	        0, 0, l, m;
	else 
		T = l, m, n, 0, 0, 0,
	        0, 0, 0, l, m, n;
}

inline void Spring::B_Matrix(LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const
{
	throw new Fatal("Spring::B_Matrix: Feature not available");
}

/* private */

inline void Spring::_initialize()
{
	if (_ndim<1) throw new Fatal("Spring::_initialize: For this element, _ndim must be greater than or equal to 1 (%d is invalid)",_ndim);
	_d  = _ndim-1; // Not used
	_nd = _ndim;   // Not used
	_nl = 3;
}

inline void Spring::_calc_initial_internal_state()
{
	throw new Fatal("Spring::_calc_initial_internal_state: Feature not available");
}

inline void Spring::OutInfo(std::ostream & os) const
{
}

inline void Spring::VTKConnect(String & Nodes) const 
{ 
	std::ostringstream os;
	for(size_t i=0; i<_connects.Size(); i++)
		os << _connects[i]->GetID() << " ";
	Nodes.Printf(os.str().c_str());
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Spring element
Element * SpringMaker()
{
	return new Spring();
}

// Register a Spring element into ElementFactory array map
int SpringRegister()
{
	ElementFactory["Spring"] = SpringMaker;
	return 0;
}

// Execute the autoregistration
int __Spring_dummy_int  = SpringRegister();


}; // namespace FEM

#endif // MECHSYS_FEM_SPRING_H
