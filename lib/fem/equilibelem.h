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

  Primary variable:    u == displacements
  Natural variable:    f == force
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
#include "models/equilibmodel.h"
#include "util/string.h"
#include "util/util.h"
#include "util/numstreams.h"
#include "linalg/laexpr.h"
#include "tensors/tensors.h"

using Util::SQ2;
using Util::_12_6;
using Tensors::Tensor1;
using Tensors::Tensor2;

namespace FEM
{

class EquilibElem : public ProbElem
{
public:
	// Typedefs
	typedef Array<double>           IntVals; ///< Internal values (specific volume, yield surface size, etc.)
    typedef std::map<String,double> Ini_t;   ///< Initial values. Ex.: Sx=0.0

	//{ Constants
	static const size_t ND_EQUILIB_3D;        ///< Number of DOFs Eq3D
	static const char   UD_EQUILIB_3D[ 3][4]; ///< Essential DOF vars 3D
	static const char   FD_EQUILIB_3D[ 3][4]; ///< Natural DOF vars 3D
	static const size_t ND_EQUILIB_2D;        ///< Number of DOFs 2D
	static const char   UD_EQUILIB_2D[ 2][4]; ///< Essential DOF vars 2D
	static const char   FD_EQUILIB_2D[ 2][4]; ///< Natural DOF vars 2D
	static const size_t NL_EQUILIB_3D;        ///< Number of labels 3D
	static const char   LB_EQUILIB_3D[16][4]; ///< Name of labels 3D
	static const size_t NL_PSTRAIN;           ///< Number of labels Plane-strain
	static const char   LB_PSTRAIN   [12][4]; ///< Name of labels Plane-strain
	static const size_t NL_PSTRESS;           ///< Number of labels Plane-stress
	static const char   LB_PSTRESS   [ 6][4]; ///< Name of labels Plane-stress
	static const char   EQUILIB_PROP [ 1][8]; ///< Properties
	//}

	// Destructor
	virtual ~EquilibElem () {}

	// Methods related to PROBLEM
	int         InitCtes     (int nDim);
	int         NProps       () const { return 1; } ///< just "gam"
	ProName_t * Props        () const { return EQUILIB_PROP; }
	void        AddVolForces ();
	void        ClearDisp    ();
	void        SetActive    (bool Activate, int ID);
	void        CalcDeps     () const;
	double      Val          (int iNod, Str_t Key) const;
	double      Val          (          Str_t Key) const;
	void        Update       (double h, Vec_t const & dU, Vec_t & dFint);
	void        Backup       ();
	void        Restore      ();
	void        OutInfo      (std::ostream & os) const;
	size_t      NCMats       () const { return 1; }
	void        CMatrix      (size_t Idx, Mat_t & M) const;
	void        CMatMap      (size_t Idx, Array<size_t> & RMap, Array<size_t> & CMap, Array<bool> & RUPresc, Array<bool> & CUPresc) const;

	// Method used when filling ProbElemFactory
	void __SetGeomIdx (int Idx) { _gi = Idx; }

protected:
	// Data at each Integration Point (IP)
	Array<Tensor2> _sig;      ///< Stress (or axial force for linear elements)
	Array<Tensor2> _eps;      ///< Strain
	Array<IntVals> _ivs;      ///< Internal values
	Array<Tensor2> _sig_bkp;  ///< Backup stress
	Array<Tensor2> _eps_bkp;  ///< Backup strain
	Array<IntVals> _ivs_bkp;  ///< Backup internal values

	// Private methods that MAY be derived
	virtual void _initialize (Str_t Inis); ///< Initialize this element
	virtual void _excavate   ();           ///< Excavate element

	// Private methods
	void _B_mat              (Mat_t const & dN, Mat_t const & J, Mat_t & B) const;      ///< Calculate B matrix
	void _dist_to_face_nodes (Str_t Key, double FaceValue, Array<Node*> const & FConn); ///< Distribute values from face/edges to nodes

private:
	void _init_internal_state (); ///< Initialize internal state

}; // class EquilibElem                                                                     

//{ Constants
const size_t EquilibElem::ND_EQUILIB_3D        = 3;
const char   EquilibElem::UD_EQUILIB_3D[ 3][4] = {"ux", "uy", "uz"};
const char   EquilibElem::FD_EQUILIB_3D[ 3][4] = {"fx", "fy", "fz"};
const size_t EquilibElem::ND_EQUILIB_2D        = 2;
const char   EquilibElem::UD_EQUILIB_2D[ 2][4] = {"ux", "uy"};
const char   EquilibElem::FD_EQUILIB_2D[ 2][4] = {"fx", "fy"};
const size_t EquilibElem::NL_EQUILIB_3D        = 16;
const char   EquilibElem::LB_EQUILIB_3D[16][4] = {"Ex", "Ey", "Ez", "Exy", "Eyz", "Ezx", "Sx", "Sy", "Sz", "Sxy", "Syz", "Szx", "p", "q", "Ev", "Ed"};
const size_t EquilibElem::NL_PSTRAIN           = 12;
const char   EquilibElem::LB_PSTRAIN   [12][4] = {"Ex", "Ey", "Ez", "Exy", "Sx", "Sy", "Sz", "Sxy", "p", "q", "Ev", "Ed"};
const size_t EquilibElem::NL_PSTRESS           = 6;
const char   EquilibElem::LB_PSTRESS   [ 6][4] = {"Ex", "Ey", "Exy", "Sx", "Sy", "Sxy"};
const char   EquilibElem::EQUILIB_PROP [ 1][8] = {"gam"};
//}


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline int EquilibElem::InitCtes(int nDim)
{
	// Check
	if (nDim<2 || nDim>3) throw new Fatal("EquilibElem::InitCtes: The space dimension must be 2 or 3. nDim==%d is invalid",nDim);

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
	else throw new Fatal("EquilibElem::InitCtes: GeometryIndex _gi==%d is invalid",_gi);

	// Return geometry index
	return _gi;
}

inline void EquilibElem::AddVolForces()
{
	// Verify if element is active
	if (IsActive==false) return;
	
	// Allocate (local/element) external volume force vector
	Vec_t fvol(_ge->NNodes);
	fvol.SetValues (0.0);

	// Calculate local external volume force
	double gam = Prop("gam");
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
	for (size_t i=0; i<_ge->NIPs; i++)
	{
		_eps    [i] = 0.0,0.0,0.0, 0.0,0.0,0.0;
		_eps_bkp[i] = 0.0,0.0,0.0, 0.0,0.0,0.0;
	}
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
	/*
	for (size_t i=0; i<_ge->NIPs; i++)
	{
		// Calculate principal values
		Tensors::Eigenvals (_sig[i], _sigp[i]);
		Tensors::Eigenvals (_eps[i], _epsp[i]);

		// Sort (increasing)
		Tensors::Sort (_sigp[i]); // S1,S2,S3 = _sigp[0], _sigp[1], _sigp[2]
		Tensors::Sort (_epsp[i]); // E1,E2,E3 = _epsp[0], _epsp[1], _epsp[2]
	}
	*/
}

inline double EquilibElem::Val(int iNod, Str_t Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[j])==0) return _ge->Conn[iNod]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[j])==0) return _ge->Conn[iNod]->DOFVar(Name).NaturalVal;

	// Stress, strains, internal values, etc.
	Vec_t    ip_vals (_ge->NIPs); // Vectors for extrapolation
	Vec_t nodal_vals (_ge->NNodes);

	// Get integration point values
	for (size_t i=0; i<_ge->NIPs; i++) ip_vals(i) = Tensors::Val (_sig[i], _eps[i], Name);

	// Extrapolate
	_ge->Extrap (ip_vals, nodal_vals);
	return nodal_vals (iNod);
}

inline double EquilibElem::Val(Str_t Name) const
{
	double ave = 0.0;
	for (size_t i=0; i<_ge->NIPs; i++) ave += Tensors::Val (_sig[i], _eps[i], Name);
	return ave/_ge->NIPs;
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
	Mat_t dN,J,B;
	Vec_t deps,dsig;
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_ge->Derivs       (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, dN);
		_ge->Jacobian     (dN, J);
		_B_mat            (dN, J, B);
		deps = B*du;
		_mdl->StateUpdate (deps, _sig[i], _eps[i], _ivs[i], dsig);
		df += trn(B)*dsig*det(J)*_ge->IPs[i].w;
	}

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		dFint(_ge->Conn[i]->DOFVar(UD[j]).EqID) += df(i*_nd+j);
}

inline void EquilibElem::Backup()
{
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_sig_bkp[i] = _sig[i];
		_eps_bkp[i] = _eps[i];
		_ivs_bkp[i] = _ivs[i];
	}
}

inline void EquilibElem::Restore()
{
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_sig[i] = _sig_bkp[i];
		_eps[i] = _eps_bkp[i];
		_ivs[i] = _ivs_bkp[i];
	}
}

inline void EquilibElem::OutInfo(std::ostream & os) const
{
	for (size_t i=0; i<_ge->NIPs; i++)
	{
		os << "IP # " << i << " Sx,Sy,Sz = " << _12_6 << _sig[i](0) << _12_6 << _sig[i](1) << _12_6 << _sig[i](2);
		os <<                "  Ex,Ey,Ez = " << _12_6 << _sig[i](0) << _12_6 << _sig[i](1) << _12_6 << _sig[i](2) << " ";
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
		_ge->Derivs       (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, dN);
		_ge->Jacobian     (dN, J);
		_B_mat            (dN, J, B);
		_mdl->TgStiffness (_sig[i], _eps[i], _ivs[i], D);
		Ke += trn(B)*D*B*det(J)*_ge->IPs[i].w;
	}
}

inline void EquilibElem::CMatMap(size_t Idx, Array<size_t> & RMap, Array<size_t> & CMap, Array<bool> & RUPresc, Array<bool> & CUPresc) const
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


/* private */

inline void EquilibElem::_initialize(Str_t Inis)
{
	// Resize IP data
	_sig    .Resize (_ge->NIPs);
	_eps    .Resize (_ge->NIPs);
	_ivs    .Resize (_ge->NIPs);
	_sig_bkp.Resize (_ge->NIPs);
	_eps_bkp.Resize (_ge->NIPs);
	_ivs_bkp.Resize (_ge->NIPs);

	// Parse values
	LineParser           lp(Inis);
    Ini_t                names_vals;
	lp.BreakExpressions (names_vals);

	// Set initial values
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		// Stress and strain
		_sig[i] = 0.0,0.0,0.0, 0.0,0.0,0.0;
		_eps[i] = 0.0,0.0,0.0, 0.0,0.0,0.0;
		for (Ini_t::const_iterator it=names_vals.begin(); it!=names_vals.end(); ++it)
			Tensors::SetVal (it->first.CStr(), it->second, _sig[i], /*WithError*/false);

		// Initialize internal values
		_mdl->InitIVS (names_vals, _sig[i], _eps[i], _ivs[i]);
	}

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
		double gam = Prop("gam");
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
			Tensor2ToVector (_gi,_sig[i], sig);

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
	Vec_t f(_ge->NDim*_ge->NNodes); // only forces (fx,fy,fz) (NDim) DOFs
	f.SetValues (0.0);

	// Calculate internal force vector;
	Mat_t dN;  // Shape Derivs
	Mat_t J;   // Jacobian matrix
	Mat_t B;   // strain-displacement matrix
	Vec_t sig; // Stress vector in Mandel's notation
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_ge->Derivs     (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, dN);
		_ge->Jacobian   (dN, J);
		_B_mat          (dN, J, B);
		Tensor2ToVector (_gi,_sig[i], sig);
		f += trn(B)*sig*det(J)*_ge->IPs[i].w;
	}

	// Assemble (local/element) displacements vector.
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (size_t j=0; j<_ge->NDim;   ++j)
		_ge->Conn[i]->DOFVar(UD[j]).NaturalVal += f(i*_ge->NDim+j); // only forces (fx,fy,fz) (NDim) DOFs
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
