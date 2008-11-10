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


/* __ Element for equilibrium simulations __

  Solves:
              dσ
              --:I + ρb = 0
              dx

  where:                      1 du    duT
             σ = D:ε       ε= -[-- + (--)]      t = σ . n
                              2 dx    dx

  Primary variale:     u == displacements
  Volumetric variable: b == body forces
  Secondary variable:  σ == stress
  Secondary variable:  ε == strain
  Secondary variable:  t == traction
  Secondary variable:  n == unit normal on boundary

*/

#ifndef MECHSYS_FEM_EQUILIB_H
#define MECHSYS_FEM_EQUILIB_H

// MechSys
#include "fem/element.h"
#include "models/equilibmodel.h"
#include "util/string.h"
#include "util/util.h"
#include "linalg/laexpr.h"

using Util::SQ2;

namespace FEM
{

class EquilibElem : public virtual Element
{
public:
	// Constants
	static const size_t ND [3];        ///< Number of DOFs to add to a node == ND[_ndim-1]
	static const size_t NDB[3];        ///< (Beam) Number of DOFs to add to a node == NDB[_ndim-1]
	static const char   UD [3][6][3];  ///< Essential DOF names == UD[_ndim-1][iDOF]
	static const char   FD [3][6][3];  ///< Natural DOF names
	static const size_t NL [5];        ///< Number of additional labels (exceeding ND)
	static const char   LB [5][18][4]; ///< Additional labels

	// Constructor
	EquilibElem () : _body_force(0.0), _has_body_force(false), _d(-1), _nd(-1) {}

	// Destructor
	virtual ~EquilibElem() {}

	// Derived methods
	virtual bool CheckModel      () const;
	bool         IsEssential     (char const * Name) const;
	virtual void SetModel        (char const * ModelName, char const * Prms, char const * Inis);
	void         SetProps        (Array<double> const & P);
	Element    * Connect         (int iNodeLocal, FEM::Node * ptNode);
	virtual void UpdateState     (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
	bool         HasVolForces    () const { return _has_body_force; }
	void         AddVolForces    (LinAlg::Vector<double> & FVol) const;
	virtual void BackupState     ();
	virtual void RestoreState    ();
	void         GetLabels       (Array<String> & Labels) const;
	void         Deactivate      ();
	char const * ModelName       () const { return (_a_model.Size()>0 ? _a_model[0]->Name() : "__no_model__"); }

	// Derived methods to assemble DAS matrices
	size_t       nOrder1Matrices () const { return 1; }
	void         Order1MatMap    (size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const;
	virtual void Order1Matrix    (size_t Index, LinAlg::Matrix<double> & Ke) const; ///< Stiffness

	// Methods
	virtual void B_Matrix (LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const;

	// Access methods
	virtual void   CalcDepVars () const;                                  ///< Calculate dependent variables (to be called before Val() or OutNodes() for example). Necessary for output of principal stresses, for example.
	virtual double Val         (int iNodeLocal, char const * Name) const; ///< Return values at nodes
	virtual double Val         (                char const * Name) const; ///< Return values at the CG of the element

protected:
	// Data
	Array<EquilibModel*> _a_model;        ///< Array of pointers to constitutive models
	Tensors::Tensor1     _body_force;     ///< Body force
	bool                 _has_body_force; ///< Has body force?
	int                  _d;              ///< Dimension index == _ndim-1
	int                  _nd;             ///< Number of DOFs == ND[_d] or NDB[_d]

	// Private methods
	virtual void _calc_initial_internal_state (); ///< Calculate initial internal state

	// Private methods that MUST be derived
	virtual int  _geom() const =0;               ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)

private:
	void _equations_map(Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const;

}; // class EquilibElem

// UD[_ndim-1][iDOF]                       1D                         2D                          3D
const size_t EquilibElem::ND [3]       = { 1,                         2,                          3};
const size_t EquilibElem::NDB[3]       = { 2,                         3,                          6};
const char   EquilibElem::UD [3][6][3] = {{"ux","wz","","","",""},  {"ux","uy","wz","","",""},  {"ux","uy","uz","wx","wy","wz"}};
const char   EquilibElem::FD [3][6][3] = {{"fx","mz","","","",""},  {"fx","fy","mz","","",""},  {"fx","fy","fz","mx","my","mz"}};

// LB[_geom-1][iLabel]
const size_t EquilibElem::NL[5]        = { 2, 12, 18, 10, 18 };
const char   EquilibElem::LB[5][18][4] = {
	{"Ea", "Sa", ""  ,  ""   , ""   , ""   , ""  , ""   , ""  , ""   , ""   , ""   , ""  , ""  , ""  , ""  , ""  , ""  }, // 1D
	{"Ex", "Ey", "Ez",  "Exy", "Sx" , "Sy" , "Sz", "Sxy", "E1", "E2" , "S1" , "S2" , ""  , ""  , ""  , ""  , ""  , ""  }, // 2D (plane-strain)
	{"Ex", "Ey", "Ez",  "Exy", "Eyz", "Ezx", "Sx", "Sy" , "Sz", "Sxy", "Syz", "Szx", "E1", "E2", "E3", "S1", "S2", "S3"}, // 3D
	{"Ex", "Ey", "Exy", "Sx" , "Sy" , "Sxy", "E1", "E2" , "S1", "S2" , ""   , ""   , ""  , ""  , ""  , ""  , ""  , ""  }, // 2D (plane-stress)
	{"Ex", "Ey", "Ez",  "Exy", "Eyz", "Ezx", "Sx", "Sy" , "Sz", "Sxy", "Syz", "Szx", "E1", "E2", "E3", "S1", "S2", "S3"}  // 2D (axis-symmetric)
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

// Derived methods

inline bool EquilibElem::CheckModel() const
{
	if (_a_model.Size()!=_a_int_pts.Size()) return false;
	for (size_t i=0; i<_a_int_pts.Size(); ++i) if (_a_model[i]==NULL) return false;
	return true;
}

inline bool EquilibElem::IsEssential(char const * Name) const
{
	for (int i=0; i<_nd; ++i) if (strcmp(Name,UD[_d][i])==0) return true;
	return false;
}

inline void EquilibElem::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// Check _ndim
	if (_ndim<1) throw new Fatal("EquilibElem::SetModel: The space dimension (SetDim) must be set before calling this method");
	if (CheckConnect()==false) throw new Fatal("EquilibElem::SetModel: Connectivity is not correct. Connectivity MUST be set before calling this method");

	// If pointers to model was not already defined => No model was allocated
	if (_a_model.Size()==0)
	{
		_a_model.Resize(_a_int_pts.Size());
		for (size_t i=0; i<_a_int_pts.Size(); ++i)
		{
			_a_model[i] = static_cast<EquilibModel*>(AllocModel(ModelName));
			_a_model[i]->SetGeom (_geom());
			_a_model[i]->SetPrms (Prms);
			_a_model[i]->SetInis (Inis);
		}
		_calc_initial_internal_state ();
	}
	else throw new Fatal("EquilibElem::SetModel: Feature not implemented.");
}

inline void EquilibElem::SetProps(Array<double> const & P)
{
	_body_force = 0.0;
	     if (P.Size()==1 && _ndim==1) _body_force(0) = P[0];
	else if (P.Size()==2 && _ndim==2) _body_force    = P[0], P[1];
	else if (P.Size()==3 && _ndim==3) _body_force    = P[0], P[1], P[2];
	else throw new Fatal("EquilibElem::SetProps: ElemProp(P)==BodyForce must have the same size as the number of dimensions (=%d)",_ndim);
}

inline Element * EquilibElem::Connect(int iNodeLocal, FEM::Node * ptNode)
{
	// Check
	if (_ndim<0 || _d<0 || _nd<0) throw new Fatal("EquilibElem::Connect: __Internal Error__: There is a problem with _ndim=%d, _d=%d, or _nd=%d\n (_ndim=space dimension, _d=dimension index==_ndim-1, and _nd=number of degrees of freedom)",_ndim,_d,_nd);

	// Connectivity
	_connects[iNodeLocal] = ptNode;

	// Add Degree of Freedom to a node (Essential, Natural)
	for (int i=0; i<_nd; ++i) _connects[iNodeLocal]->AddDOF (UD[_d][i], FD[_d][i]);

	// Set shared
	_connects[iNodeLocal]->SetSharedBy(_my_id);

	return this;
}

inline void EquilibElem::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
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

	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs;  // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;       // Jacobian matrix
	LinAlg::Matrix<double> B;       // strain-displacement matrix
	LinAlg::Vector<double> deps;    // delta strain vector
	LinAlg::Vector<double> dsig;    // delta stress vector

	// Update model and calculate internal force vector;
	for (size_t i=0; i<_a_int_pts.Size(); ++i)
	{
		double r = _a_int_pts[i].r;
		double s = _a_int_pts[i].s;
		double t = _a_int_pts[i].t; // only for 3D cases
		double w = _a_int_pts[i].w;

		Derivs   (r,s,t, derivs);  // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);      // Calculate J (Jacobian) matrix for i Integration Point
		B_Matrix (derivs, J, B);   // Calculate B matrix for i Integration Point

		deps = B*du;
		_a_model[i]->StateUpdate(deps, dsig);
		df += trn(B)*dsig*det(J)*w;
	}

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_n_nodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		dFint(_connects[i]->DOFVar(UD[_d][j]).EqID) += df(i*_ndim+j);
}

inline void EquilibElem::AddVolForces(LinAlg::Vector<double> & FVol) const
{
	if (_has_body_force)
	{
		// Allocate (local/element) external volume force vector
		LinAlg::Vector<double> fvol(_nd*_n_nodes);
		fvol.SetValues(0.0);

		// Allocate entities used for every integration point
		LinAlg::Vector<double> shape;
		LinAlg::Matrix<double> derivs;
		LinAlg::Matrix<double> J;

		// Calculate local external volume force
		for (size_t i=0; i<_a_int_pts.Size(); ++i)
		{
			double r = _a_int_pts[i].r;
			double s = _a_int_pts[i].s;
			double t = _a_int_pts[i].t;
			double w = _a_int_pts[i].w;

			Shape    (r,s,t, shape);   // Calculate shape functions for i IP
			Derivs   (r,s,t, derivs);  // Calculate Derivatives of Shape functions w.r.t local coordinate system
			Jacobian (derivs, J);      // Calculate J (Jacobian) matrix for i Integration Point

			for (size_t j=0; j<_n_nodes; ++j)
			for (size_t k=0; k<3;        ++k)
				fvol(j*_ndim+k) += _body_force(k)*shape(j)*det(J)*w;
		}

		// Sum up contribution to internal forces vector
		for (size_t i=0; i<_n_nodes; ++i)
		for (size_t j=0; j<3;        ++j)
			FVol(_connects[i]->DOFVar(UD[_d][j]).EqID) += fvol(i*_ndim+j);
	}
	else throw new Fatal("EquilibElem::AddVolForces: This element (%s # %d) does not have volumetric forces.",Name(),_my_id);
}

inline void EquilibElem::BackupState()
{
	for (size_t i=0; i<_a_int_pts.Size(); ++i) _a_model[i]->BackupState();
}

inline void EquilibElem::RestoreState()
{
	for (size_t i=0; i<_a_int_pts.Size(); ++i) _a_model[i]->RestoreState();
}

inline void EquilibElem::GetLabels(Array<String> & Labels) const
{
	const int p  = _geom()-1;
	const int q  = NL[p];   // number of additional labels
	const int nl = 2*_nd+q; // total number of labels
	Labels.Resize(nl);
	size_t k = 0;
	for (int i=0; i<_nd; ++i)
	{
		Labels[k] = UD[_d][i];  k++;
		Labels[k] = FD[_d][i];  k++;
	}
	for (int i=0; i<q; ++i)
	{
		Labels[k] = LB[p][i];  k++;
	}
}

inline void EquilibElem::CalcDepVars() const
{
	if (_a_model.Size()==_a_int_pts.Size()) for (size_t i=0; i<_a_int_pts.Size(); i++) _a_model[i]->CalcDepVars();
	else throw new Fatal("EquilibElem::CalcDepVars: Constitutive models for this element (ID==%d) were not set yet", _my_id);
}

inline double EquilibElem::Val(int iNodeLocal, char const * Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).NaturalVal;

	// Stress, strains, internal values, etc.
	LinAlg::Vector<double>    ip_values (_a_int_pts.Size()); // Vectors for extrapolation
	LinAlg::Vector<double> nodal_values (_n_nodes);

	// Get integration point values
	if (_a_model.Size()==_a_int_pts.Size())
		for (size_t i=0; i<_a_int_pts.Size(); i++)
			ip_values(i) = _a_model[i]->Val(Name);
	else throw new Fatal("EquilibElem::Val: Constitutive models for this element (ID==%d) were not set yet", _my_id);

	Extrapolate (ip_values, nodal_values);
	return nodal_values (iNodeLocal);
}

inline double EquilibElem::Val(char const * Name) const
{
	// Get integration point values
	double sum = 0.0;
	if (_a_model.Size()==_a_int_pts.Size())
		for (size_t i=0; i<_a_int_pts.Size(); i++)
			sum += _a_model[i]->Val(Name);
	else throw new Fatal("EquilibElem::Val: Constitutive models for this element (ID==%d) were not set yet", _my_id);

	// Output single value at CG
	return sum/_a_int_pts.Size();
}

inline void EquilibElem::Deactivate()
{
	throw new Fatal("EquilibElem::Deactivate: Feature not implemented yet");
}


// Derived methods to assemble DAS matrices

inline void EquilibElem::Order1MatMap(size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const
{
	_equations_map(RowsMap, ColsMap, RowsEssenPresc, ColsEssenPresc);
}

inline void EquilibElem::Order1Matrix(size_t index, LinAlg::Matrix<double> & Ke) const
{
	/* Stiffness:
	   ==========

	                 /    T
	        [Ke]  =  | [B]  * [D] * [B]  * dV
	                 /
	*/

	// Resize Ke
	Ke.Resize(_nd*_n_nodes, _nd*_n_nodes);
	Ke.SetValues(0.0);

	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs; // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;      // Jacobian matrix
	LinAlg::Matrix<double> B;      // strain-displacement matrix
	LinAlg::Matrix<double> D;      // Constitutive matrix

	// Calculate Tangent Stiffness
	for (size_t i=0; i<_a_int_pts.Size(); ++i)
	{
		double r = _a_int_pts[i].r;
		double s = _a_int_pts[i].s;
		double t = _a_int_pts[i].t;
		double w = _a_int_pts[i].w;

		Derivs   (r,s,t, derivs); // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);     // Calculate J (Jacobian) matrix for i Integration Point
		B_Matrix (derivs,J, B);   // Calculate B matrix for i Integration Point

		_a_model[i]->TgStiffness(D);
		Ke += trn(B)*D*B*det(J)*w;
	}
}


// Methods

inline void EquilibElem::B_Matrix(LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const
{
	/* OBS.:
	 *          This B matrix considers Soil Mechanics sign convention of stress and strains
	 *          Ex.: Compressive stresses/strains are positive
	 *          The B Matrix returns strains in Mandel notation
	 */

	// geometry type: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	switch (_geom())
	{
		case 1: // 1D
		{
			// Derivatives and determinand of Jacobian
			LinAlg::Matrix<double> nat_derivs(_ndim, _nd*_n_nodes);
			nat_derivs.SetValues(0.0);
			double det_J = det(J);
			for (size_t i=0; i<_n_nodes; i++)
			for (int    j=0; j<_nd;      j++)
				nat_derivs(j, i*_ndim+j) = derivs(0,i);
			// Assemble B matrix
			B = -1.0/(det_J*det_J)*J*nat_derivs; // B matrix for a linear element in 1D, 2D and 3D.
			return;
		}
		case 2: // 2D(plane-strain)
		{
			// Cartesian derivatives
			LinAlg::Matrix<double> cart_derivs;
			cart_derivs = inv(J)*derivs;
			// Resize B matrix
			const int n_scomps = 4; // number of stress compoments
			B.Resize(n_scomps,_nd*_n_nodes);
			// Loop along all nodes of the element
			double dNdX,dNdY;
			int  j=0; // j column of B
			for (size_t i=0; i<_n_nodes; ++i) // i row of B
			{
				// Assemble B matrix
				j = i*_ndim;
				dNdX=-cart_derivs(0,i);  dNdY=-cart_derivs(1,i);  // Negative values => Soil mechanics convention
				B(0,0+j) =      dNdX;   B(0,1+j) =  0.0;
				B(1,0+j) =       0.0;   B(1,1+j) = dNdY;
				B(2,0+j) =       0.0;   B(2,1+j) =  0.0;
				B(3,0+j) =  dNdY/SQ2;   B(3,1+j) = dNdX/SQ2;  // SQ2 => Mandel representation
			}
			return;
		}
		case 3: // 3D
		{
			// Cartesian derivatives
			LinAlg::Matrix<double> cart_derivs;
			cart_derivs = inv(J)*derivs;
			// Resize B matrix
			const int n_scomps = 6; // number of stress compoments
			B.Resize(n_scomps,_nd*_n_nodes);
			// Loop along all nodes of the element
			double dNdX,dNdY,dNdZ;
			int  j=0; // j column of B
			for (size_t i=0; i<_n_nodes; ++i) // i row of B
			{
				// Assemble B matrix
				j = i*_ndim;
				dNdX=-cart_derivs(0,i);  dNdY=-cart_derivs(1,i);  dNdZ=-cart_derivs(2,i); // Negative values => Soil mechanics convention
				B(0,0+j) =     dNdX;     B(0,1+j) =      0.0;     B(0,2+j) =      0.0;
				B(1,0+j) =      0.0;     B(1,1+j) =     dNdY;     B(1,2+j) =      0.0;
				B(2,0+j) =      0.0;     B(2,1+j) =      0.0;     B(2,2+j) =     dNdZ;
				B(3,0+j) = dNdY/SQ2;     B(3,1+j) = dNdX/SQ2;     B(3,2+j) =      0.0; // SQ2 => Mandel representation
				B(4,0+j) =      0.0;     B(4,1+j) = dNdZ/SQ2;     B(4,2+j) = dNdY/SQ2; // SQ2 => Mandel representation
				B(5,0+j) = dNdZ/SQ2;     B(5,1+j) =      0.0;     B(5,2+j) = dNdX/SQ2; // SQ2 => Mandel representation
			}
			return;
		}
		case 5: // 2D(plane-stress)
		{
			// Cartesian derivatives
			LinAlg::Matrix<double> cart_derivs;
			cart_derivs = inv(J)*derivs;
			// Resize B matrix
			const int n_scomps = 3; // number of stress compoments
			B.Resize(n_scomps,_nd*_n_nodes);
			// Loop along all nodes of the element
			double dNdX,dNdY;
			int  j=0; // j column of B
			for (size_t i=0; i<_n_nodes; ++i) // i row of B
			{
				// Assemble B matrix
				j = i*_ndim;
				dNdX=-cart_derivs(0,i);  dNdY=-cart_derivs(1,i);  // Negative values => Soil mechanics convention
				B(0,0+j) =      dNdX;   B(0,1+j) =  0.0;
				B(1,0+j) =       0.0;   B(1,1+j) = dNdY;
				B(2,0+j) =  dNdY/SQ2;   B(2,1+j) = dNdX/SQ2;  // SQ2 => Mandel representation
			}
			return;
		}
		case 4: // 2D(axis-symmetric)
		default:
			throw new Fatal("EquilibElem::B_Matrix: GeometryType==%d is not implemented yet",_geom());
	}
}


/* private */

inline void EquilibElem::_calc_initial_internal_state()
{
	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> f(_nd*_n_nodes);
	f.SetValues(0.0);

	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs;  // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;       // Jacobian matrix
	LinAlg::Matrix<double> B;       // strain-displacement matrix
	LinAlg::Vector<double> sig;     // Stress vector in Mandel's notation

	// Calculate internal force vector;
	for (size_t i=0; i<_a_int_pts.Size(); ++i)
	{
		double r = _a_int_pts[i].r;
		double s = _a_int_pts[i].s;
		double t = _a_int_pts[i].t; // only for 3D
		double w = _a_int_pts[i].w;

		Derivs   (r,s,t, derivs); // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);     // Calculate J (Jacobian) matrix for i Integration Point
		B_Matrix (derivs, J, B);  // Calculate B matrix for i Integration Point

		_a_model[i]->Sig(sig);
		f += trn(B)*sig*det(J)*w;
	}

	// Assemble (local/element) displacements vector.
	for (size_t i=0; i<_n_nodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		_connects[i]->DOFVar(UD[_d][j]).NaturalVal += f(i*_ndim+j);
}

inline void EquilibElem::_equations_map(Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const
{
	// Size of Stiffness/Mass Matrix
	int n_rows = _nd*_n_nodes; // == n_cols

	// Mounting a map of positions from Me to Global
	RowsMap       .Resize(n_rows);
	RowsEssenPresc.Resize(n_rows);

	int p = 0; // position inside matrix
	for (size_t i=0; i<_n_nodes; ++i)
	{
		for (int j=0; j<_nd; ++j)
		{
			RowsMap        [p] = _connects[i]->DOFVar(UD[_d][j]).EqID;
			RowsEssenPresc [p] = _connects[i]->DOFVar(UD[_d][j]).IsEssenPresc;
			p++;
		}
	}
	ColsMap        = RowsMap;
	ColsEssenPresc = RowsEssenPresc;
}

}; // namespace FEM

#endif // MECHSYS_FEM_EQUILIB_H
