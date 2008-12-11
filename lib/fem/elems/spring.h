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
	Spring() : _ks(-1) { NNodes=2; Conn.Resize(NNodes); Conn.SetValues(NULL);}
	

	// Derived methods
	char const * Name() const { return "Spring"; }

	// Derived methods
	bool   CheckModel   () const;
	void   SetModel     (char const * ModelName, char const * Prms, char const * Inis);
	void   UpdateState  (double TimeInc, Vec_t const & dUglobal, Vec_t & dFint);
	void   CalcDepVars  () const;
	double Val          (int iNodeLocal, char const * Name) const;
	double Val          (char const * Name) const;
	void   Order1Matrix (size_t Index, Mat_t & Ke) const; ///< Stiffness
	void   B_Matrix     (Mat_t const & derivs, Mat_t const & J, Mat_t & B) const;
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

inline void Spring::UpdateState(double TimeInc, Vec_t const & dUglobal, Vec_t & dFint)
{
	// Allocate (local/element) displacements vector
	Vec_t du(_nd*NNodes); // Delta disp. of this element

	// Assemble (local/element) displacements vector
	for (size_t i=0; i<NNodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		du(i*_nd+j) = dUglobal(Conn[i]->DOFVar(UD[_d][j]).EqID);

	// Allocate (local/element) internal force vector
	Vec_t df; // Delta internal force of this element

	Mat_t Ke;
	Order1Matrix(0,Ke);
	df = Ke * du;

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<NNodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		dFint(Conn[i]->DOFVar(UD[_d][j]).EqID) += df(i*_nd+j);
}

inline void Spring::CalcDepVars() const
{
}

inline double Spring::Val(int iNodeLocal, char const * Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[_d][j])==0) return Conn[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[_d][j])==0) return Conn[iNodeLocal]->DOFVar(Name).NaturalVal;

	     if (strcmp(Name,"Ea")==0)  return 0.0;
	else if (strcmp(Name,"Sa")==0)  return 0.0;
	else if (strcmp(Name,"N" )==0)
	{
		// Allocate (local/element) displacements vector
		Vec_t du(_nd*NNodes); // Delta disp. of this element
		// Assemble (local/element) displacements vector
		for (size_t i=0; i<NNodes; ++i)
		for (int    j=0; j<_nd;      ++j)
		du(i*_nd+j) = Conn[i]->Val(UD[_d][j]);

		Matrix<double> T; _mount_T_matrix(T);
		Vector<double> D; D = T*du;
		return (D(1)-D(0))*_ks;
	}
	else throw new Fatal("Spring::Val: This element does not have a Val named %s",Name);
}

inline double Spring::Val(char const * Name) const
{
	if (strcmp(Name,"N" )==0)
	{
		// Allocate (local/element) displacements vector
		Vec_t du(_nd*NNodes); // Delta disp. of this element
		// Assemble (local/element) displacements vector
		for (size_t i=0; i<NNodes; ++i)
		for (int    j=0; j<_nd;      ++j)
		du(i*_nd+j) = Conn[i]->Val(UD[_d][j]);

		Matrix<double> T; _mount_T_matrix(T);
		Vector<double> D; D = T*du;
		return (D(1)-D(0))*_ks;
	}
	else throw new Fatal("Spring::Val: This element does not have a Val named %s",Name);
}

inline void Spring::Order1Matrix(size_t Index, Mat_t & Ke) const
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
	double x0 = Conn[0]->X(); double y0 = Conn[0]->Y(); double z0 = Conn[0]->Z();
	double x1 = Conn[1]->X(); double y1 = Conn[1]->Y(); double z1 = Conn[1]->Z();
	double L  = sqrt(pow(x1-x0,2.0)+pow(y1-y0,2.0)+pow(z1-z0,2.0));
	double l = (x1-x0)/L;
	double m = (y1-y0)/L;
	double n = (z1-z0)/L;

	// Mount T Matrix
	T.Resize(2, _ndim*2);
	if (_ndim==2)
		T = l, m, 0, 0, 
	        0, 0, l, m;
	else 
		T = l, m, n, 0, 0, 0,
	        0, 0, 0, l, m, n;
}

inline void Spring::B_Matrix(Mat_t const & derivs, Mat_t const & J, Mat_t & B) const
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
	for(size_t i=0; i<Conn.Size(); i++)
		os << Conn[i]->GetID() << " ";
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
