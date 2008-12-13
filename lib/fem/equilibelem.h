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

#ifndef MECHSYS_FEM_EQUILIBELEM_H
#define MECHSYS_FEM_EQUILIBELEM_H

// MechSys
#include "fem/probelem.h"
#include "fem/equilibvars.h"
#include "models/equilibmodel.h"
#include "util/string.h"
#include "util/util.h"
#include "util/numstreams.h"
#include "linalg/laexpr.h"

using Util::SQ2;
using Util::_12_6;

namespace FEM
{

class EquilibElem : public ProbElem
{
public:
	// Destructor
	virtual ~EquilibElem () { for (size_t i=0; i<_mdl.Size(); ++i) delete _mdl[i]; }

	// Methods related to PROBLEM
	void    AddVolForces ();
	void    ClearDisp    ();
	void    SetActive    (bool Activate, int ID);
	void    CalcDeps     () const;
	Str_t   ModelName    () const { return (_mdl.Size()>0 ? _mdl[0]->Name() : "__no_model__"); }
	double  Val          (int iNod, Str_t Key) const;
	double  Val          (          Str_t Key) const;
	void    Update       (double h, Vec_t const & dU, Vec_t & dFint);
	void    Backup       ();
	void    Restore      ();
	void    OutInfo      (std::ostream & os) const;
	size_t  NCMats       () const { return 1; }
	void    CMatrix      (size_t Idx, Mat_t & M) const;
	void    CMatMap      (size_t Idx, Array<size_t> & RMap, Array<size_t> & CMap, Array<bool> & RUPresc, Array<bool> & CUPresc) const;

	// Derived methods
	virtual bool CheckModel () const;

	// Method used when filling ProbElemFactory
	void __SetGeomIdx (int Idx) { _gi = Idx; }

protected:
	// Data
	Array<EquilibModel*> _mdl; ///< Array of pointers to constitutive models

	// Private methods that MAY be derived
	virtual void _initialize (Str_t Model); ///< Initialize this element
	virtual void _init_model (Str_t Inis);  ///< Initialize model
	virtual void _excavate   ();            ///< Excavate element

	// Private methods
	void _B_mat              (Mat_t const & dN, Mat_t const & J, Mat_t & B) const;      ///< Calculate B matrix
	void _dist_to_face_nodes (Str_t Key, double FaceValue, Array<Node*> const & FConn); ///< Distribute values from face/edges to nodes

private:
	void _init_internal_state (); ///< Initialize internal state
	void _equations_map       (Array<size_t> & RMap, Array<size_t> & CMap, Array<bool> & RUPresc, Array<bool> & CUPresc) const;

}; // class EquilibElem                                                                     


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void EquilibElem::AddVolForces()
{
	// Verify if element is active
	if (IsActive==false) return;
	
	// Allocate (local/element) external volume force vector
	Vec_t fvol(_ge->NNodes);
	fvol.SetValues (0.0);

	// Calculate local external volume force
	if (_props.Size()!=1) throw new Fatal("EquilibElem::AddVolForces: Array of properties must be set with gam==20, for example");
	double gam = _props[0];
	Vec_t N;
	Mat_t dN;
	Mat_t J;
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_ge->Shape    (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, N);
		_ge->Derivs   (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, dN);
		_ge->Jacobian (dN, J);
		fvol += -N*gam*det(J)*_ge->IPs[i].w;
	}

	// Sum up contribution to external forces vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	{
		if (_ge->NDim==2) _ge->Conn[i]->Bry("fy",fvol(i));
		if (_ge->NDim==3) _ge->Conn[i]->Bry("fz",fvol(i));
	}
}

inline void EquilibElem::ClearDisp()
{
	if (IsActive==false) return;

	// Clear displacements
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		_ge->Conn[i]->DOFVar(UD[j]).EssentialVal = 0.0;

	// Clear strains
	for (size_t i=0; i<_mdl.Size(); ++i) _mdl[i]->ClearStrain();
}

inline void EquilibElem::SetActive(bool Activate, int ID)
{
	if (IsActive==false && Activate)
	{
		// Set active
		IsActive = true;

		for (size_t i=0; i<_ge->Conn.Size(); ++i)
		{
			// Add Degree of Freedom to a node (Essential, Natural)
			for (int j=0; j<_nd; ++j) _ge->Conn[i]->AddDOF (UD[j], FD[j]);

			// Set SharedBy
			_ge->Conn[i]->SetSharedBy (ID);
		}

		// Apply body forces
		AddVolForces ();
	}
	if (IsActive && Activate==false)
	{
		// Set active
		IsActive = false;

		for (size_t i=0; i<_ge->Conn.Size(); ++i)
		{
			// Remove SharedBy
			_ge->Conn[i]->RemoveSharedBy (ID);

			// Remove Degree of Freedom to a node (Essential)
			if (_ge->Conn[i]->nSharedBy()==0) 
				for (int j=0; j<_nd; ++j) _ge->Conn[i]->RemoveDOF (UD[j]);
		}

		// Apply surface tensions to perform excavation
		_excavate();
	}
}

inline void EquilibElem::CalcDeps() const
{
	if (IsActive==false) throw new Fatal("EquilibElem::CalcDepVars: This element is inactive");
	if (_mdl.Size()==_ge->NIPs) for (size_t i=0; i<_ge->NIPs; i++) _mdl[i]->CalcDepVars();
	else throw new Fatal("EquilibElem::CalcDepVars: Constitutive models for this element were not set");
}

inline double EquilibElem::Val(int iNod, Str_t Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[j])==0) return _ge->Conn[iNod]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[j])==0) return _ge->Conn[iNod]->DOFVar(Name).NaturalVal;

	// Stress, strains, internal values, etc.
	Vec_t    ip_values (_ge->NIPs); // Vectors for extrapolation
	Vec_t nodal_values (_ge->NNodes);

	// Get integration point values
	if (_mdl.Size()==_ge->NIPs)
		for (size_t i=0; i<_ge->NIPs; i++) ip_values(i) = _mdl[i]->Val(Name);
	else throw new Fatal("EquilibElem::Val: Constitutive models for this element were not set yet");

	// Extrapolate
	_ge->Extrap (ip_values, nodal_values);
	return nodal_values (iNod);
}

inline double EquilibElem::Val(Str_t Name) const
{
	// Get integration point values
	double sum = 0.0;
	if (_mdl.Size()==_ge->NIPs)
		for (size_t i=0; i<_ge->NIPs; i++) sum += _mdl[i]->Val(Name);
	else throw new Fatal("EquilibElem::Val: Constitutive models for this element were not set yet");

	// Output single value at CG
	return sum/_ge->NIPs;
}

inline void EquilibElem::Update(double h, Vec_t const & dU, Vec_t & dFint)
{
	// Allocate (local/element) displacements vector
	Vec_t du(_nd*_ge->NNodes); // Delta disp. of this element

	// Assemble (local/element) displacements vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		du(i*_nd+j) = dU(_ge->Conn[i]->DOFVar(UD[j]).EqID);

	// Allocate (local/element) internal force vector
	Vec_t df(_nd*_ge->NNodes); // Delta internal force of this element
	df.SetValues (0.0);

	// Update model and calculate internal force vector;
	Mat_t dN;   // Shape derivatives
	Mat_t J;    // Jacobian
	Mat_t B;    // B matrix
	Vec_t deps; // Delta Strain
	Vec_t dsig; // Delta Stress
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_ge->Derivs   (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, dN);
		_ge->Jacobian (dN, J);
		_B_mat        (dN, J, B);
		deps = B*du;
		_mdl[i]->StateUpdate (deps, dsig);
		df += trn(B)*dsig*det(J)*_ge->IPs[i].w;
	}

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		dFint(_ge->Conn[i]->DOFVar(UD[j]).EqID) += df(i*_nd+j);
}

inline void EquilibElem::Backup()
{
	for (size_t i=0; i<_ge->NIPs; ++i) _mdl[i]->BackupState();
}

inline void EquilibElem::Restore()
{
	for (size_t i=0; i<_ge->NIPs; ++i) _mdl[i]->RestoreState();
}

inline void EquilibElem::OutInfo(std::ostream & os) const
{
	for (size_t i=0; i<_ge->NIPs; i++)
	{
		os << "IP # " << i << " Sx,Sy,Sz = " << _12_6 << _mdl[i]->Val("Sx") << _12_6 << _mdl[i]->Val("Sy") << _12_6 << _mdl[i]->Val("Sz");
		os <<                "  Ex,Ey,Ez = " << _12_6 << _mdl[i]->Val("Ex") << _12_6 << _mdl[i]->Val("Ey") << _12_6 << _mdl[i]->Val("Ez") << " ";
	}
}

inline void EquilibElem::CMatrix(size_t Idx, Mat_t & Ke) const
{
	/* Stiffness:
	   ==========
	                 /    T
	        [Ke]  =  | [B]  * [D] * [B]  * dV
	                 /
	*/

	Ke.Resize    (_nd*_ge->NNodes, _nd*_ge->NNodes);
	Ke.SetValues (0.0);
	Mat_t dN,J,B,D;
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_ge->Derivs   (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, dN);
		_ge->Jacobian (dN, J);
		_B_mat        (dN, J, B);
		_mdl[i]->TgStiffness (D);
		Ke += trn(B)*D*B*det(J)*_ge->IPs[i].w;
	}
}

inline void EquilibElem::CMatMap(size_t Idx, Array<size_t> & RMap, Array<size_t> & CMap, Array<bool> & RUPresc, Array<bool> & CUPresc) const
{
	_equations_map (RMap, CMap, RUPresc, CUPresc);
}

inline bool EquilibElem::CheckModel() const
{
	if (_prms.Size()<1)                                   return false;
	if (_mdl.Size()!=_ge->NIPs)                           return false;
	for (size_t i=0; i<_ge->NIPs; ++i) if (_mdl[i]==NULL) return false;
	return true;
}


/* private */

inline void EquilibElem::_initialize(Str_t Model)
{
	// Check
	if (_ge->NDim<1) throw new Fatal("EquilibElem::_initialize: The space dimension (SetDim) must be set before calling this method");

	// Set number of DOFs (_nd), number of labels (_nl), and arrays of essential UD, natural FD, and labels LB
	if (_gi==0)  // 3D
	{
		_nd = ND_EQUILIB_3D;
		UD  = UD_EQUILIB_3D;
		FD  = FD_EQUILIB_3D;
		_nl = NL_EQUILIB_3D; 
		LB  = LB_EQUILIB_3D;
	}
	else if (_gi==1)  // PlaneStrain
	{
		_nd = ND_EQUILIB_2D;
		UD  = UD_EQUILIB_2D;
		FD  = FD_EQUILIB_2D;
		_nl = NL_PSTRAIN; 
		LB  = LB_PSTRAIN;
	}
	else if (_gi==2)  // PlaneStress
	{
		_nd = ND_EQUILIB_2D;
		UD  = UD_EQUILIB_2D;
		FD  = FD_EQUILIB_2D;
		_nl = NL_PSTRESS; 
		LB  = LB_PSTRESS;
	}
	else throw new Fatal("EquilibElem::_initialize: GeometryIndex _gi==%d is invalid",_gi);

	// Model and parameters
	_mdl.Resize (_ge->NIPs);  _mdl.SetValues(NULL);
	for (size_t i=0; i<_ge->NIPs; ++i) _mdl[i] = static_cast<EquilibModel*>(AllocModel(Model));
	_prms.Resize (_mdl[0]->NPrms());
	PRMS = _mdl[0]->GetPrmNames();

	// Properties
	_props.Resize (1); // just "gam"
	PROP = EQUILIB_PROP;
}

inline void EquilibElem::_init_model(Str_t Inis)
{
	// Set model initial state
	if (_mdl.Size()!=_ge->NIPs)  throw new Fatal("EquilibElem::_init_state: Array with pointers to model was not properly initialized");
	for (size_t i=0; i<_ge->NIPs; ++i) _mdl[i]->Initialize (_gi, &_prms, Inis);

	// Initialize internal state
	if (IsActive) _init_internal_state();
}

inline void EquilibElem::_excavate()
{
	// Verify if element is on the boundary of excavation
	bool on_boundary = false;
	for (size_t i=0; i<_ge->NNodes; i++)
	{
		if (_ge->Conn[i]->nSharedBy()>0)
		{
			on_boundary=true;
			break;
		}
	}

	// Calculate internal forces for element on boundary
	if (on_boundary)
	{
		if (_props.Size()!=1) throw new Fatal("EquilibElem::_excavate: Array of properties must be set with gam==20, for example");
		double gam = _props[0];
		Vec_t F(_ge->NDim*_ge->NNodes);  F.SetValues(0.0);
		Vec_t S(_ge->NDim*_ge->NNodes);  S.SetValues(0.0);
		Vec_t N;      // Shape functions
		Mat_t dN;     // Derivative of shape functions
		Mat_t J;      // Jacobian matrix
		Mat_t B;      // strain-displacement matrix
		Vec_t sig;    // stress tensor
		for (size_t i=0; i<_ge->NIPs; ++i)
		{
			_ge->Shape    (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, N);
			_ge->Derivs   (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, dN);
			_ge->Jacobian (dN, J);
			_B_mat        (dN, J, B);

			// Mount S
			for (size_t j=0; j<_ge->NNodes; ++j) S(_ge->NDim*j+_ge->NDim-1) = N(j);

			// Get tensor
			_mdl[i]->Sig (sig);

			// Calculate internal force vector
			F += trn(B)*sig*det(J)*_ge->IPs[i].w + S*gam*det(J)*_ge->IPs[i].w;
		}

		// Add to boundary values of node if nSharedBy()>0
		for (size_t i=0; i<_ge->NNodes; ++i)
		{
			if (_ge->Conn[i]->nSharedBy()>0)
			{
				// Apply forces (only for equilibrium type neighbor elements)
				for (size_t j=0; j<_ge->NDim; ++j)
					_ge->Conn[i]->Bry (FD[j], F(_ge->NDim*i+j));
			}
		}
	}
}

inline void EquilibElem::_B_mat(Mat_t const & dN, Mat_t const & J, Mat_t & B) const
{
	/* OBS.:
	 *       1) This B matrix considers Solid Mechanics sign convention of stress and strains
	 *          Ex.: Compressive stresses/strains are negative
	 *          The B Matrix returns strains in Mandel notation
	 *
	 *          Traction    => Positive
	 *          Compression => Negative
	 *
	 *       2) This works for Equilib and Biot elements, but not for Beams and Rods
	 */

	// Cartesian derivatives
	Mat_t dC;
	dC = inv(J)*dN;

	int dim = _ge->NDim;
	switch (_gi)
	{
		case 0: // 3D
		{
			const int n_scomps = 6; // number of stress compoments
			B.Resize (n_scomps,dim*_ge->NNodes);
			for (size_t i=0; i<_ge->NNodes; ++i) // i row of B
			{
				B(0,0+i*dim) =     dC(0,i);  B(0,1+i*dim) =         0.0;  B(0,2+i*dim) =         0.0;
				B(1,0+i*dim) =         0.0;  B(1,1+i*dim) =     dC(1,i);  B(1,2+i*dim) =         0.0;
				B(2,0+i*dim) =         0.0;  B(2,1+i*dim) =         0.0;  B(2,2+i*dim) =     dC(2,i);
				B(3,0+i*dim) = dC(1,i)/SQ2;  B(3,1+i*dim) = dC(0,i)/SQ2;  B(3,2+i*dim) =         0.0; // SQ2 => Mandel representation
				B(4,0+i*dim) =         0.0;  B(4,1+i*dim) = dC(2,i)/SQ2;  B(4,2+i*dim) = dC(1,i)/SQ2; // SQ2 => Mandel representation
				B(5,0+i*dim) = dC(2,i)/SQ2;  B(5,1+i*dim) =         0.0;  B(5,2+i*dim) = dC(0,i)/SQ2; // SQ2 => Mandel representation
			}
			return;
		}
		case 1: // 2D(plane-strain)
		{
			const int n_scomps = 4; // number of stress compoments
			B.Resize (n_scomps,dim*_ge->NNodes);
			for (size_t i=0; i<_ge->NNodes; ++i) // i row of B
			{
				B(0,0+i*dim) =     dC(0,i);  B(0,1+i*dim) =         0.0;
				B(1,0+i*dim) =         0.0;  B(1,1+i*dim) =     dC(1,i);
				B(2,0+i*dim) =         0.0;  B(2,1+i*dim) =         0.0;
				B(3,0+i*dim) = dC(1,i)/SQ2;  B(3,1+i*dim) = dC(0,i)/SQ2; // SQ2 => Mandel representation
			}
			return;
		}
		case 2: // 2D(plane-stress)
		{
			const int n_scomps = 3; // number of stress compoments
			B.Resize(n_scomps,dim*_ge->NNodes);
			for (size_t i=0; i<_ge->NNodes; ++i) // i row of B
			{
				B(0,0+i*dim) =      dC(0,i);   B(0,1+i*dim) =         0.0;
				B(1,0+i*dim) =          0.0;   B(1,1+i*dim) =     dC(1,i);
				B(2,0+i*dim) =  dC(1,i)/SQ2;   B(2,1+i*dim) = dC(0,i)/SQ2; // SQ2 => Mandel representation
			}
			return;
		}
		case 3: // 2D(axis-symmetric)
		default: throw new Fatal("EquilibElem::_B_mat: _B_mat() method is not available for GeometryIndex(gi)==%d",_gi);
	}
}

inline void EquilibElem::_equations_map(Array<size_t> & RMap, Array<size_t> & CMap, Array<bool> & RUPresc, Array<bool> & CUPresc) const
{
	// Map of positions from Me to Global
	RMap   .Resize(_nd*_ge->NNodes);
	RUPresc.Resize(_nd*_ge->NNodes);
	int p = 0; // position inside matrix
	for (size_t i=0; i<_ge->NNodes; ++i)
	{
		for (int j=0; j<_nd; ++j)
		{
			RMap    [p] = _ge->Conn[i]->DOFVar(UD[j]).EqID;
			RUPresc [p] = _ge->Conn[i]->DOFVar(UD[j]).IsEssenPresc;
			p++;
		}
	}
	CMap    = RMap;
	CUPresc = RUPresc;
}

inline void EquilibElem::_dist_to_face_nodes(Str_t Key, double const FaceValue, Array<Node*> const & FConn)
{
	// Key=="Q" => Normal traction boundary condition
	if (strcmp(Key,"Q")==0)
	{
		// Check if the element is active
		if (IsActive==false) return;

		Mat_t values;  values.Resize (_ge->NFNodes, _ge->NDim);  values.SetValues(0.0);
		Mat_t J;
		Vec_t FN(_ge->NFNodes); // Face shape
		Mat_t FNmat;            // Shape function matrix
		Vec_t P;                // Vector perpendicular to the face
		for (size_t i=0; i<_ge->NFIPs; i++)
		{
			_ge->FaceShape (_ge->FIPs[i].r, _ge->FIPs[i].s, FN);
			FNmat = trn(trn(FN)); // trick just to convert Vector FN to a col Matrix

			// Calculate perpendicular vector
			if (_ge->NDim==3)
			{
				_ge->FaceJacob (FConn, _ge->FIPs[i].r, _ge->FIPs[i].s, J);
				Vec_t V(3); V = J(0,0), J(0,1), J(0,2);
				Vec_t W(3); W = J(1,0), J(1,1), J(1,2);
				P.Resize(3);
				P = V(1)*W(2) - V(2)*W(1),      // vectorial product
					V(2)*W(0) - V(0)*W(2),
					V(0)*W(1) - V(1)*W(0);
			}
			else
			{
				_ge->FaceJacob (FConn, _ge->FIPs[i].r, _ge->FIPs[i].s, J);
				P.Resize(2);
				P = J(0,1), -J(0,0);
			}
			values += FaceValue*FNmat*trn(P)*_ge->FIPs[i].w;
		}

		// Set nodes Brys
		for (size_t i=0; i<_ge->NFNodes; ++i)
		for (size_t j=0; j<_ge->NDim;    ++j)
			FConn[i]->Bry (FD[j], values(i,j));
	}
	else ProbElem::_dist_to_face_nodes (Key,FaceValue,FConn);
}

inline void EquilibElem::_init_internal_state()
{
	// Allocate (local/element) internal force vector
	Vec_t f(_nd*_ge->NNodes);
	f.SetValues (0.0);

	// Calculate internal force vector;
	Mat_t dN;  // Shape Derivs
	Mat_t J;   // Jacobian matrix
	Mat_t B;   // strain-displacement matrix
	Vec_t sig; // Stress vector in Mandel's notation
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_ge->Derivs   (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, dN);
		_ge->Jacobian (dN, J);
		_B_mat        (dN, J, B);
		_mdl[i]->Sig (sig);
		f += trn(B)*sig*det(J)*_ge->IPs[i].w;
	}

	// Assemble (local/element) displacements vector.
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		_ge->Conn[i]->DOFVar(UD[j]).NaturalVal += f(i*_nd+j);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new 3D Equilib element:
ProbElem * EquilibMaker() 
{ 
	EquilibElem * Ptr = new EquilibElem;
	Ptr->__SetGeomIdx(0);
	return Ptr; 
}
// Register element
int EquilibRegister() { ProbElemFactory["Equilib"]=EquilibMaker;  return 0; }
// Call register
int __Equilib_dummy_int  = EquilibRegister();


// Allocate a new PStrain element:
ProbElem * PStrainMaker() 
{ 
	EquilibElem * Ptr = new EquilibElem;
	Ptr->__SetGeomIdx(1);
	return Ptr; 
}
// Register element
int PStrainRegister() { ProbElemFactory["PStrain"]=PStrainMaker;  return 0; }
// Call register
int __PStrain_dummy_int  = PStrainRegister();


// Allocate a new PStress element:
ProbElem * PStressMaker() 
{ 
	EquilibElem * Ptr = new EquilibElem;
	Ptr->__SetGeomIdx(2);
	return Ptr; 
}
// Register element
int PStressRegister() { ProbElemFactory["PStress"]=PStressMaker;  return 0; }
// Call register
int __PStress_dummy_int  = PStressRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_EQUILIBELEM_H
