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
#include "util/exception.h"
#include "fem/elems/vtkCellType.h"

namespace FEM
{

class Rod3 : public EquilibElem
{
public:
	// Constructor
	Rod3 () : _gam(0.0), _E(-1), _A(-1) { NNodes=3; Conn.Resize(NNodes); Conn.SetValues(NULL); SetIntPoints(2);}

	// Derived methods
	char const * Name() const { return "Rod3"; }

	// Derived methods
	void   Shape(double r, double s, double t, Vec_t & N) const;
	void   Derivs(double r, double s, double t, Mat_t & Derivs) const;
	void   LocalCoords(Mat_t & coords) const;
	void   SetIntPoints (int NumGaussPoints);
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
	int    VTKCellType  () const { return VTK_QUADRATIC_EDGE; }
	void   VTKConnect   (String & Nodes) const { Nodes.Printf("%d %d %d",Conn[0]->GetID(),Conn[2]->GetID(),Conn[1]->GetID()); }
	void   OutInfo(std::ostream & os) const;

	// Methods
	double N(double l) const; ///< Axial force (0 < l < 1) (Must be used after CalcDepVars())

private:
	// Data
	double _gam; ///< Specific weigth
	double _E; ///< Young modulus
	double _A; ///< Cross-sectional area

	// Depedent variables (calculated by CalcDepVars)
	mutable Vector<double> _Eps;   ///< Strain at nodal points
	mutable Vector<double> _EpsIP; ///< Strain at integration points

	// Private methods
	int  _geom                        () const { return 1; }              ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void _initialize                  ();                                 ///< Initialize the element
	void _calc_initial_internal_state ();                                 ///< Calculate initial internal state
	void _transf_mat                  (Mat_t & T) const; ///< Calculate transformation matrix

}; // class Rod3


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */
inline void  Rod3::Shape(double r, double s, double t, Vec_t & N) const 
{

	/*
	 *       -----o=========o=========o----->  r
	 *            0         1         2
	 *
	 */

	N.Resize(3);
	N(0) = 0.5*(r*r-r);
	N(1) = 1.0 - r*r;
	N(2) = 0.5*(r*r+r);
}

inline void Rod3::Derivs(double r, double s, double t, Mat_t & Derivs) const 
{
	Derivs.Resize(1,3);
	Derivs(0,0) =  r-0.5;
	Derivs(0,1) = -2.0*r;
	Derivs(0,2) =  r+0.5;
}

inline void Rod3::LocalCoords(Mat_t & coords) const 
{
	if (_ndim==2)
	{
		coords.Resize(3,3);
		coords =  -1.0,  0.0, 1.0,
				   0.0,  0.0, 1.0,
				   1.0,  0.0, 1.0;
	}
	else
	{
		coords.Resize(3,4);
		coords =  -1.0,  0.0, 0.0, 1.0,
				   0.0,  0.0, 0.0, 1.0,
				   1.0,  0.0, 0.0, 1.0;
	}
}

inline void Rod3::SetIntPoints(int NumGaussPoints)
{
	// Setup pointer to the array of Integration Points
	if      (NumGaussPoints==2) _a_int_pts = LIN_IP2;
	else if (NumGaussPoints==3) _a_int_pts = LIN_IP3;
	else throw new Fatal("Rod3::SetIntPoints: Error in number of integration points.");

	_n_int_pts      = NumGaussPoints;
	_n_face_int_pts = 0;
}

inline bool Rod3::CheckModel() const
{
	if (_E<0.0 || _A<0.0) return false;
	return true;
}

inline void Rod3::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// Check _ndim
	if (_ndim<1) throw new Fatal("Rod3::SetModel: The space dimension (SetDim) must be set before calling this method");
	if (CheckConnect()==false) throw new Fatal("Rod3::SetModel: Connectivity is not correct. Connectivity MUST be set before calling this method");

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
		else throw new Fatal("Rod3::SetModel: Parameter name (%s) is invalid",names[i].CStr());
	}
}

inline void Rod3::SetProps(char const * Properties)
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

inline void Rod3::UpdateState(double TimeInc, Vec_t const & dUglobal, Vec_t & dFint)
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

inline void Rod3::ApplyBodyForces() // TODO
{
	// Verify if element is active
	if (_is_active==false) return;

	// Weight
	double dx = Conn[1]->X()-Conn[0]->X();
	double dy = Conn[1]->Y()-Conn[0]->Y();
	double L  = sqrt(dx*dx+dy*dy);
	double W  = _A*L*_gam;

	// Set boundary conditions
	if (_ndim==1) throw new Fatal("Rod3::ApplyBodyForces: feature not available for NDim==1");
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

inline void Rod3::CalcDepVars() const
{
	// Calculate the strain at nodes
	_EpsIP.Resize(_n_int_pts);

	// Assemble (local/element) displacements vector
	Vec_t U(_nd*NNodes); 
	for (size_t i=0; i<NNodes; ++i) for (int j=0; j<_nd; ++j)
		U(i*_nd+j) = Conn[i]->Val(UD[_d][j]);

	// Calculate Strain
	Mat_t derivs; // size = NumLocalCoords(ex.: r,s,t) x NNodes
	Mat_t J;      // Jacobian matrix
	Mat_t B;      // strain-displacement matrix

	// Loop along integration points
	for (size_t i=0; i<_n_int_pts; ++i)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i].r;
		Derivs(r, 0, 0, derivs);      // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian(derivs, J);          // Calculate J (Jacobian) matrix for i Integration Point
		B_Matrix(derivs,J, B);        // Calculate B matrix for i Integration Point
		Vector<double> C;
		C = B*U;
		_EpsIP(i)=C(0);
	}
	Extrapolate(_EpsIP, _Eps);
}

inline double Rod3::Val(int iNodeLocal, char const * Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[_d][j])==0) return Conn[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[_d][j])==0) return Conn[iNodeLocal]->DOFVar(Name).NaturalVal;

	if (_Eps.Size()<1) throw new Fatal("Rod3::Val: Please, call CalcDepVars() before calling this method");

	     if (strcmp(Name,"Ea")==0) return _Eps(iNodeLocal);
	else if (strcmp(Name,"Sa")==0) return _Eps(iNodeLocal)*_E;
	else if (strcmp(Name,"N" )==0) return _Eps(iNodeLocal)*_E*_A;
	else throw new Fatal("Rod3::Val: This element does not have a Val named %s",Name);
}

inline double Rod3::Val(char const * Name) const
{
	throw new Fatal("Rod3::Val: Feature not available");
}

inline void Rod3::OutInfo(std::ostream & os) const
{
	CalcDepVars();
	for (size_t i=0; i<_n_int_pts; i++)
	{
		os << "IP # " << i << " Ea,Sa = " << _12_6 << _EpsIP(i) << _12_6 << _EpsIP(i)*_E;
	}
}


inline void Rod3::Order1Matrix(size_t Index, Mat_t & Ke) const
{
	// Resize Ke
	Ke.Resize(_ndim*NNodes, _ndim*NNodes); // sum(Bt*D*B*det(J)*w)
	Ke.SetValues(0.0);

	// Allocate entities used for every integration point
	Mat_t derivs; // size = NumLocalCoords(ex.: r,s,t) x NNodes
	Mat_t J;      // Jacobian matrix
	Mat_t B;      // strain-displacement matrix

	// Loop along integration points
	for (size_t i_ip=0; i_ip<_n_int_pts; ++i_ip)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i_ip].r;
		double w = _a_int_pts[i_ip].w;

		Derivs(r, 0, 0, derivs);      // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian(derivs, J);          // Calculate J (Jacobian) matrix for i_ip Integration Point
		B_Matrix(derivs,J, B);        // Calculate B matrix for i_ip Integration Point

		// Calculate Tangent Stiffness
		Ke += trn(B)*_E*B*det(J)*_A*w;
	}
}

inline void Rod3::B_Matrix(Mat_t const & derivs, Mat_t const & J, Mat_t & B) const
{
	Mat_t cart_derivs(_ndim, _ndim*NNodes); // cartesian derivatives
	double det_J = det(J);

	cart_derivs.SetValues(0.0);
	for(size_t i=0; i<NNodes; i++) for(int j=0; j<_ndim; j++) cart_derivs(j, i*_ndim+j) = derivs(0,i);

	B = 1.0/(det_J*det_J)*J*cart_derivs; // B matrix for a rod element (2D and 3D)
}

/* private */

inline void Rod3::_initialize()
{
	if (_ndim<1) throw new Fatal("Rod3::_initialize: For this element, _ndim must be greater than or equal to 1 (%d is invalid)",_ndim);
	_d  = _ndim-1;
	_nd = EquilibElem::ND[_d];
	_nl = EquilibElem::NL[_geom()-1];
}

inline void Rod3::_calc_initial_internal_state()
{
	throw new Fatal("Rod3::_calc_initial_internal_state: Feature not available");
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Rod3 element
Element * Rod3Maker()
{
	return new Rod3();
}

// Register a Rod3 element into ElementFactory array map
int Rod3Register()
{
	ElementFactory["Rod3"] = Rod3Maker;
	return 0;
}

// Execute the autoregistration
int __Rod3_dummy_int  = Rod3Register();

}; // namespace FEM

#endif // MECHSYS_FEM_ROD3_H
