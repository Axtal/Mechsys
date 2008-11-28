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
#include "util/numstreams.h"
#include "linalg/laexpr.h"

using Util::SQ2;
using Util::_12_6;

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
	EquilibElem () : _gam(0.0), _d(-1), _nd(-1) {}

	// Destructor
	virtual ~EquilibElem();

	// Derived methods
	virtual bool CheckModel          () const;
	bool         IsEssential         (char const * Name) const;
	virtual void SetModel            (char const * ModelName, char const * Prms, char const * Inis);
	void         SetProps            (char const * Properties);
	Element    * Connect             (int iNodeLocal, FEM::Node * ptNode);
	virtual void UpdateState         (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
	void         ApplyBodyForces     ();
	virtual void BackupState         ();
	virtual void RestoreState        ();
	void         GetLabels           (Array<String> & Labels) const;
	char const * ModelName           () const { return (_a_model.Size()>0 ? _a_model[0]->Name() : "__no_model__"); }
	virtual void ClearDispAndStrains ();
	void         SetActive           (bool Active);
	virtual void OutInfo             (std::ostream & os) const;

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
	Array<EquilibModel*> _a_model;  ///< Array of pointers to constitutive models
	double               _gam;      ///< Specific weigth
	int                  _d;        ///< Dimension index == _ndim-1
	int                  _nd;       ///< Number of DOFs == ND[_d] or NDB[_d]

	// Private methods
	virtual void _calc_initial_internal_state (); ///< Calculate initial internal state

	// Private methods that MUST be derived
	virtual int _geom() const =0; ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)

private:
	void _equations_map(Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const;
	void _dist_to_face_nodes(char const * Key, double const FaceValue, Array<Node*> const & FaceConnects) const;

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

inline EquilibElem::~EquilibElem()
{
	for (size_t i=0; i<_a_model.Size(); ++i) delete _a_model[i];
}

// Derived methods

inline bool EquilibElem::CheckModel() const
{
	if (_a_model.Size()!=_n_int_pts) return false;
	for (size_t i=0; i<_n_int_pts; ++i) if (_a_model[i]==NULL) return false;
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
		_a_model.Resize(_n_int_pts);
		for (size_t i=0; i<_n_int_pts; ++i)
		{
			_a_model[i] = static_cast<EquilibModel*>(AllocModel(ModelName));
			_a_model[i]->SetGeom (_geom());
			_a_model[i]->SetPrms (Prms);
			_a_model[i]->SetInis (Inis);
		}
		if (_is_active) _calc_initial_internal_state ();
	}
	else throw new Fatal("EquilibElem::SetModel: Feature not implemented.");
}

inline void EquilibElem::SetProps(char const * Properties)
{
	/* "gam=20.0 */
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

inline Element * EquilibElem::Connect(int iNodeLocal, FEM::Node * ptNode)
{
	// Check
	if (_n_nodes<1)               throw new Fatal("EquilibElem::Connect: __Internal Error__: There is a problem with the number of nodes: maybe derived elemet did not set _n_nodes");
	if (_connects.Size()<1)       throw new Fatal("EquilibElem::Connect: __Internal Error__: There is a problem with connectivity array: maybe derived elemet did not allocate _connect");
	if (_ndim<0 || _d<0 || _nd<0) throw new Fatal("EquilibElem::Connect: __Internal Error__: There is a problem with _ndim=%d, _d=%d, or _nd=%d\n (_ndim=space dimension, _d=dimension index==_ndim-1, and _nd=number of degrees of freedom)",_ndim,_d,_nd);

	// Connectivity
	_connects[iNodeLocal] = ptNode;

	if (_is_active)
	{
		// Add Degree of Freedom to a node (Essential, Natural)
		for (int i=0; i<_nd; ++i) _connects[iNodeLocal]->AddDOF (UD[_d][i], FD[_d][i]);

		// Set shared
		_connects[iNodeLocal]->SetSharedBy(_my_id);
	}

	return this;
}

inline void EquilibElem::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
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

	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs;  // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;       // Jacobian matrix
	LinAlg::Matrix<double> B;       // strain-displacement matrix
	LinAlg::Vector<double> deps;    // delta strain vector
	LinAlg::Vector<double> dsig;    // delta stress vector

	// Update model and calculate internal force vector;
	for (size_t i=0; i<_n_int_pts; ++i)
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
		dFint(_connects[i]->DOFVar(UD[_d][j]).EqID) += df(i*_nd+j);
}

inline void EquilibElem::ApplyBodyForces() 
{
	// Verify if element is active
	if (_is_active==false) return;
	
	// Allocate (local/element) external volume force vector
	LinAlg::Matrix<double> fvol(_n_nodes, _ndim);
	fvol.SetValues(0.0);

	// Allocate entities used for every integration point
	LinAlg::Vector<double> shape;
	LinAlg::Matrix<double> derivs;
	LinAlg::Matrix<double> J;
	LinAlg::Vector<double> b;

	// Mounting the body force vector
	     if (_ndim==3) { b.Resize(3); b = 0.0, 0.0, -_gam; }
	else if (_ndim==2) { b.Resize(2); b = 0.0, -_gam; }
	else if (_ndim==1) { b.Resize(1); b = -_gam; }

	// Calculate local external volume force
	for (size_t i=0; i<_n_int_pts; ++i)
	{
		double r = _a_int_pts[i].r;
		double s = _a_int_pts[i].s;
		double t = _a_int_pts[i].t;
		double w = _a_int_pts[i].w;

		Shape    (r,s,t, shape);   // Calculate shape functions for i IP
		Derivs   (r,s,t, derivs);  // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);      // Calculate J (Jacobian) matrix for i Integration Point

		fvol += shape*trn(b)*det(J)*w;
	}

	// Sum up contribution to external forces vector
	for (size_t i=0; i<_n_nodes; ++i)
	{
		              _connects[i]->Bry("fx",fvol(i,0));
		              _connects[i]->Bry("fy",fvol(i,1));
		if (_ndim==3) _connects[i]->Bry("fz",fvol(i,2));
	}
}

inline void EquilibElem::BackupState()
{
	for (size_t i=0; i<_n_int_pts; ++i) _a_model[i]->BackupState();
}

inline void EquilibElem::RestoreState()
{
	for (size_t i=0; i<_n_int_pts; ++i) _a_model[i]->RestoreState();
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
	if (_a_model.Size()==_n_int_pts) for (size_t i=0; i<_n_int_pts; i++) _a_model[i]->CalcDepVars();
	else throw new Fatal("EquilibElem::CalcDepVars: Constitutive models for this element (ID==%d) were not set yet", _my_id);
}

inline double EquilibElem::Val(int iNodeLocal, char const * Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).NaturalVal;

	// Stress, strains, internal values, etc.
	LinAlg::Vector<double>    ip_values (_n_int_pts); // Vectors for extrapolation
	LinAlg::Vector<double> nodal_values (_n_nodes);

	// Get integration point values
	if (_a_model.Size()==_n_int_pts)
		for (size_t i=0; i<_n_int_pts; i++)
			ip_values(i) = _a_model[i]->Val(Name);
	else throw new Fatal("EquilibElem::Val: Constitutive models for this element (ID==%d) were not set yet", _my_id);

	Extrapolate (ip_values, nodal_values);
	return nodal_values (iNodeLocal);
}

inline double EquilibElem::Val(char const * Name) const
{
	// Get integration point values
	double sum = 0.0;
	if (_a_model.Size()==_n_int_pts)
		for (size_t i=0; i<_n_int_pts; i++)
			sum += _a_model[i]->Val(Name);
	else throw new Fatal("EquilibElem::Val: Constitutive models for this element (ID==%d) were not set yet", _my_id);

	// Output single value at CG
	return sum/_n_int_pts;
}

inline void EquilibElem::ClearDispAndStrains()
{
	if (_is_active==false) return;

	// Clear displacements
	for (size_t i=0; i<_n_nodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		_connects[i]->DOFVar(UD[_d][j]).EssentialVal = 0.0;

	// Clear strains
	for (size_t i=0; i<_a_model.Size(); ++i) _a_model[i]->ClearStrain();
}

inline void EquilibElem::SetActive(bool Activate)
{
	if (_is_active==false && Activate)
	{
		// Set active
		_is_active = true;

		for (size_t i=0; i<_connects.Size(); ++i)
		{
			// Add Degree of Freedom to a node (Essential, Natural)
			for (int j=0; j<_nd; ++j) _connects[i]->AddDOF (UD[_d][j], FD[_d][j]);

			// Set SharedBy
			_connects[i]->SetSharedBy (_my_id);
		}

		// Apply body forces
		ApplyBodyForces ();
	}
	else throw new Fatal("EquilibElem::SetActive: Deactivation (excavation) is not available yet");
}

inline void EquilibElem::OutInfo(std::ostream & os) const
{
	for (size_t i=0; i<_n_int_pts; i++)
	{
		os << "IP # " << i << " Sx,Sy,Sz = " << _12_6 << _a_model[i]->Val("Sx") << _12_6 << _a_model[i]->Val("Sy") << _12_6 << _a_model[i]->Val("Sz");
		os <<                "  Ex,Ey,Ez = " << _12_6 << _a_model[i]->Val("Ex") << _12_6 << _a_model[i]->Val("Ey") << _12_6 << _a_model[i]->Val("Ez") << " ";
	}
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
	for (size_t i=0; i<_n_int_pts; ++i)
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
	 *          This B matrix considers Solid Mechanics sign convention of stress and strains
	 *          Ex.: Compressive stresses/strains are negative
	 *          The B Matrix returns strains in Mandel notation
	 *
	 *          Traction    => Positive
	 *          Compression => Negative
	 */

	// Cartesian derivatives
	LinAlg::Matrix<double> dN;
	dN = inv(J)*derivs;

	// geometry type: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	switch (_geom())
	{
		case 2: // 2D(plane-strain)
		{
			const int n_scomps = 4; // number of stress compoments
			B.Resize (n_scomps,_nd*_n_nodes);
			for (size_t i=0; i<_n_nodes; ++i) // i row of B
			{
				B(0,0+i*_nd) =     dN(0,i);  B(0,1+i*_nd) =         0.0;
				B(1,0+i*_nd) =         0.0;  B(1,1+i*_nd) =     dN(1,i);
				B(2,0+i*_nd) =         0.0;  B(2,1+i*_nd) =         0.0;
				B(3,0+i*_nd) = dN(1,i)/SQ2;  B(3,1+i*_nd) = dN(0,i)/SQ2; // SQ2 => Mandel representation
			}
			return;
		}
		case 3: // 3D
		{
			const int n_scomps = 6; // number of stress compoments
			B.Resize (n_scomps,_nd*_n_nodes);
			for (size_t i=0; i<_n_nodes; ++i) // i row of B
			{
				B(0,0+i*_nd) =     dN(0,i);  B(0,1+i*_nd) =         0.0;  B(0,2+i*_nd) =         0.0;
				B(1,0+i*_nd) =         0.0;  B(1,1+i*_nd) =     dN(1,i);  B(1,2+i*_nd) =         0.0;
				B(2,0+i*_nd) =         0.0;  B(2,1+i*_nd) =         0.0;  B(2,2+i*_nd) =     dN(2,i);
				B(3,0+i*_nd) = dN(1,i)/SQ2;  B(3,1+i*_nd) = dN(0,i)/SQ2;  B(3,2+i*_nd) =         0.0; // SQ2 => Mandel representation
				B(4,0+i*_nd) =         0.0;  B(4,1+i*_nd) = dN(2,i)/SQ2;  B(4,2+i*_nd) = dN(1,i)/SQ2; // SQ2 => Mandel representation
				B(5,0+i*_nd) = dN(2,i)/SQ2;  B(5,1+i*_nd) =         0.0;  B(5,2+i*_nd) = dN(0,i)/SQ2; // SQ2 => Mandel representation
			}
			return;
		}
		case 5: // 2D(plane-stress)
		{
			const int n_scomps = 3; // number of stress compoments
			B.Resize(n_scomps,_nd*_n_nodes);
			for (size_t i=0; i<_n_nodes; ++i) // i row of B
			{
				B(0,0+i*_nd) =      dN(0,i);   B(0,1+i*_nd) =         0.0;
				B(1,0+i*_nd) =          0.0;   B(1,1+i*_nd) =     dN(1,i);
				B(2,0+i*_nd) =  dN(1,i)/SQ2;   B(2,1+i*_nd) = dN(0,i)/SQ2; // SQ2 => Mandel representation
			}
			return;
		}
		case 1: // 1D
		case 4: // 2D(axis-symmetric)
		default: throw new Fatal("EquilibElem::B_Matrix: B_Matrix() method is not available for GeometryType==%d",_geom());
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
	for (size_t i=0; i<_n_int_pts; ++i)
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
		_connects[i]->DOFVar(UD[_d][j]).NaturalVal += f(i*_nd+j);
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

inline void EquilibElem::_dist_to_face_nodes(char const * Key, double const FaceValue, Array<Node*> const & FaceConnects) const
{
	// Compute face nodal values (integration along the face)

	// Conventional face boundary condition
	if (strcmp(Key,"Q")!=0)
	{
		Element::_dist_to_face_nodes(Key, FaceValue, FaceConnects);
		return;
	}

	// Check if the element is active
	if (_is_active==false) return;

	// Normal traction boundary condition
	LinAlg::Matrix<double> values;  values.Resize(_n_face_nodes, _ndim);  values.SetValues(0.0);
	LinAlg::Matrix<double> J;                         // Jacobian matrix. size = [1,2] x 3
	LinAlg::Vector<double> face_shape(_n_face_nodes); // Shape functions of a face/edge. size = _n_face_nodes
	LinAlg::Matrix<double> F;                         // Shape function matrix
	LinAlg::Vector<double> P;                         // Vector perpendicular to the face 
	for (size_t i=0; i<_n_face_int_pts; i++)
	{
		double r = _a_face_int_pts[i].r;
		double s = _a_face_int_pts[i].s;
		double w = _a_face_int_pts[i].w;
		FaceShape    (r, s, face_shape);
		F = trn(trn(face_shape)); // trick just to convert Vector face_shape to a col Matrix

		// Calculate perpendicular vector
		if (_ndim==3)
		{
			FaceJacobian (FaceConnects, r, s, J);
			LinAlg::Vector<double> V(3); V = J(0,0), J(0,1), J(0,2);
			LinAlg::Vector<double> W(3); W = J(1,0), J(1,1), J(1,2);
			P.Resize(3);
			P = V(1)*W(2) - V(2)*W(1),      // vectorial product
			    V(2)*W(0) - V(0)*W(2),
			    V(0)*W(1) - V(1)*W(0);
		}
		else
		{
			FaceJacobian (FaceConnects, r, J);
			P.Resize(2);
			P = J(0,1), -J(0,0);  
		}
		values += FaceValue*F*trn(P)*w;
	}

	// Set nodes Brys
	for (size_t i=0; i<_n_face_nodes; ++i)
	{
		              FaceConnects[i]->Bry("fx",values(i,0));
		              FaceConnects[i]->Bry("fy",values(i,1));
		if (_ndim==3) FaceConnects[i]->Bry("fz",values(i,2));
	}
}

}; // namespace FEM

#endif // MECHSYS_FEM_EQUILIB_H
