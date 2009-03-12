/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo, Fernando Alonso                   *
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

/* __ Element for diffusion transport simulations __

  Solves:
             dv                       d   du
           - --:I + s = 0    ==    - --(k.--):I = s
             dx                      dx   dx

  where:
                              du
            v = -k.i      i = --      qn = v . n
                              dx

  Primary variable:    u  == temperature/total head              ~~ displacements
  Natural variable:    f  == flow(volume)                        ~~ force
  Volumetric variable: s  == heat source/wate recharge (pumping) ~~ body forces
  Secondary variable:  v  == heat flux/Darcy's velocity          ~~ stress
  Secondary variable:  i  == gradient                            ~~ strain
  Secondary variable:  qn == normal flow                         ~~ traction
  Secondary variable:  n  == unit normal on boundary

*/

#ifndef MECHSYS_FEM_DIFFUSIONELEM_H
#define MECHSYS_FEM_DIFFUSIONELEM_H

// MechSys
#include "fem/probelem.h"
#include "models/diffusionmodel.h"
#include "util/string.h"
#include "util/util.h"
#include "util/numstreams.h"
#include "linalg/laexpr.h"

using Util::SQ2;
using Util::_12_6;
using Tensors::Tensor1;

namespace FEM
{

class DiffusionElem : public ProbElem
{
public:
	// Typedefs
	typedef Array<double>                 IntVals; ///< Internal values (specific volume, yield surface size, etc.)
	typedef blitz::TinyVector<double,3>   Vec3_t;
	typedef blitz::TinyMatrix<double,3,3> Mat3_t;
    typedef std::map<String,double>       Ini_t;   ///< Initial values. Ex.: Sx=0.0

	//{ Constants
	static const size_t ND_DIFFUSION;          ///< Number of DOFs
	static const char   UD_DIFFUSION   [1][4]; ///< Essential DOF vars
	static const char   FD_DIFFUSION   [1][4]; ///< Natural DOF vars
	static const size_t NL_DIFFUSION_3D;       ///< Number of labels 3D
	static const char   LB_DIFFUSION_3D[6][4]; ///< Name of labels 3D
	static const size_t NL_DIFFUSION_2D;       ///< Number of labels 2D
	static const char   LB_DIFFUSION_2D[4][4]; ///< Name of labels 2D
	static const char   DIFFUSION_PROP [1][8]; ///< Properties
	//}

	// Destructor
	virtual ~DiffusionElem () {}

	// Methods related to PROBLEM
	int         InitCtes     (int nDim);
	int         NProps       () const { return 1; } ///< just "s" (source)
	ProName_t * Props        () const { return DIFFUSION_PROP; }
	void        AddVolForces ();
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

protected:
	// Data at each Integration Point (IP)
	Array<Vec3_t>  _vel;      ///< Stress (or axial force for linear elements)
	Array<Vec3_t>  _gra;      ///< Strain
	Array<IntVals> _ivs;      ///< Internal values
	Array<Vec3_t>  _vel_bkp;  ///< Backup stress
	Array<Vec3_t>  _gra_bkp;  ///< Backup strain
	Array<IntVals> _ivs_bkp;  ///< Backup internal values

	// Private methods that MAY be derived
	virtual void _initialize (Str_t Inis); ///< Initialize this element

	// Private methods
	void _B_mat              (Mat_t const & dN, Mat_t const & J, Mat_t & B) const;      ///< Calculate B matrix
	void _dist_to_face_nodes (Str_t Key, double FaceValue, Array<Node*> const & FConn); ///< Distribute values from face/edges to nodes

private:
	void _init_internal_state (); ///< Initialize internal state

}; // class DiffusionElem                                                                     

//{ Constants
const size_t DiffusionElem::ND_DIFFUSION           = 1;
const char   DiffusionElem::UD_DIFFUSION   [1][4] = {"u"};
const char   DiffusionElem::FD_DIFFUSION   [1][4] = {"f"};
const size_t DiffusionElem::NL_DIFFUSION_3D       = 6;
const char   DiffusionElem::LB_DIFFUSION_3D[6][4] = {"Vx", "Vy", "Vz", "Ix", "Iy", "Iz"};
const size_t DiffusionElem::NL_DIFFUSION_2D       = 4;
const char   DiffusionElem::LB_DIFFUSION_2D[4][4] = {"Vx", "Vy", "Ix", "Iy"};
const char   DiffusionElem::DIFFUSION_PROP [1][8] = {"s"};
//}


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline int DiffusionElem::InitCtes(int nDim)
{
	// Check
	if (nDim<2 || nDim>3) throw new Fatal("DiffusionElem::InitCtes: The space dimension must be 2 or 3. nDim==%d is invalid",nDim);
	_gi = (nDim==2 ? 1 : 0);

	// Essential/Natural
	_nd = ND_DIFFUSION;
	UD  = UD_DIFFUSION;
	FD  = FD_DIFFUSION;

	// Set number of DOFs (_nd), number of labels (_nl), and arrays of essential UD, natural FD, and labels LB
	if (_gi==0)  // 3D
	{
		_nl = NL_DIFFUSION_3D; 
		LB  = LB_DIFFUSION_3D;
	}
	else if (_gi==1)  // 2D
	{
		_nl = NL_DIFFUSION_2D; 
		LB  = LB_DIFFUSION_2D;
	}
	else throw new Fatal("DiffusionElem::InitCtes: GeometryIndex _gi==%d is invalid",_gi);

	// Return geometry index
	return _gi;
}

inline void DiffusionElem::AddVolForces()
{
	// Verify if element is active
	if (IsActive==false) return;
	
	// Allocate (local/element) external volume force vector
	Vec_t fvol(_ge->NNodes);
	fvol.SetValues (0.0);

	// Calculate local external volume force
	double s = Prop("s");
	Vec_t N;
	Mat_t dN;
	Mat_t J;
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_ge->Shape    (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, N);
		_ge->Derivs   (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, dN);
		_ge->Jacobian (dN, J);
		fvol += s*N*det(J)*_ge->IPs[i].w;
	}

	// Sum up contribution to external forces vector
	for (size_t i=0; i<_ge->NNodes; ++i)
		_ge->Conn[i]->Bry("f",fvol(i));
}

inline void DiffusionElem::SetActive(bool Activate, int ID)
{
	throw new Fatal("DiffusionElem::SetActive: Element # %d: Method not available yet",ID);
}

inline void DiffusionElem::CalcDeps() const
{
	if (IsActive==false) throw new Fatal("DiffusionElem::CalcDepVars: This element is inactive");
}

inline double DiffusionElem::Val(int iNod, Str_t Name) const
{
	// Essential
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[j])==0) return _ge->Conn[iNod]->DOFVar(Name).EssentialVal;

	// Natural
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[j])==0) return _ge->Conn[iNod]->DOFVar(Name).NaturalVal;

	// Veloc, grads, internal values, etc.
	Vec_t    ip_vals (_ge->NIPs); // Vectors for extrapolation
	Vec_t nodal_vals (_ge->NNodes);

	// Get integration point values
	//for (size_t i=0; i<_ge->NIPs; i++) ip_vals(i) = Tensors::Val (_vel[i], _gra[i], Name);

	// Extrapolate
	_ge->Extrap (ip_vals, nodal_vals);
	return nodal_vals (iNod);
}

inline double DiffusionElem::Val(Str_t Name) const
{
	double ave = 0.0;
	//for (size_t i=0; i<_ge->NIPs; i++) ave += Tensors::Val (_vel[i], _gra[i], Name);
	return ave/_ge->NIPs;
}

inline void DiffusionElem::Update(double h, Vec_t const & dU, Vec_t & dFint)
{
	// Allocate (local/element) displacements vector
	Vec_t du(_nd*_ge->NNodes); // Delta temp/head of this element

	// Assemble (local/element) displacements vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		du(i*_nd+j) = dU(_ge->Conn[i]->DOFVar(UD[j]).EqID);

	// Allocate (local/element) internal force vector
	Vec_t df(_nd*_ge->NNodes); // Delta internal force of this element
	df.SetValues (0.0);

	// Update model and calculate internal force vector;
	Mat_t dN,J,B;
	Vec_t dgrad,dvel;
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_ge->Derivs       (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, dN);
		_ge->Jacobian     (dN, J);
		_B_mat            (dN, J, B);
		dgrad = B*du;
		_mdl->StateUpdate (dgrad, _vel[i], _gra[i], _ivs[i], dvel);
		df += trn(B)*dvel*det(J)*_ge->IPs[i].w;
	}

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		dFint(_ge->Conn[i]->DOFVar(UD[j]).EqID) += df(i*_nd+j);
}

inline void DiffusionElem::Backup()
{
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_vel_bkp[i] = _vel[i];
		_gra_bkp[i] = _gra[i];
		_ivs_bkp[i] = _ivs[i];
	}
}

inline void DiffusionElem::Restore()
{
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_vel[i] = _vel_bkp[i];
		_gra[i] = _gra_bkp[i];
		_ivs[i] = _ivs_bkp[i];
	}
}

inline void DiffusionElem::OutInfo(std::ostream & os) const
{
	for (size_t i=0; i<_ge->NIPs; i++)
	{
		os << "IP # " << i << " Sx,Sy,Sz = " << _12_6 << _vel[i](0) << _12_6 << _vel[i](1) << _12_6 << _vel[i](2);
		os <<                "  Ex,Ey,Ez = " << _12_6 << _vel[i](0) << _12_6 << _vel[i](1) << _12_6 << _vel[i](2) << " ";
	}
}

inline void DiffusionElem::CMatrix(size_t Idx, Mat_t & Ke) const
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
		_ge->Derivs          (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, dN);
		_ge->Jacobian        (dN, J);
		_B_mat               (dN, J, B);
		_mdl->TgConductivity (_vel[i], _gra[i], _ivs[i], D);
		Ke += trn(B)*D*B*det(J)*_ge->IPs[i].w;
	}
}

inline void DiffusionElem::CMatMap(size_t Idx, Array<size_t> & RMap, Array<size_t> & CMap, Array<bool> & RUPresc, Array<bool> & CUPresc) const
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

inline void DiffusionElem::_initialize(Str_t Inis)
{
	// Resize IP data
	_vel    .Resize (_ge->NIPs);
	_gra    .Resize (_ge->NIPs);
	_ivs    .Resize (_ge->NIPs);
	_vel_bkp.Resize (_ge->NIPs);
	_gra_bkp.Resize (_ge->NIPs);
	_ivs_bkp.Resize (_ge->NIPs);

	// Parse values
	LineParser           lp(Inis);
    Ini_t                names_vals;
	lp.BreakExpressions (names_vals);

	// Set initial values
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		// Stress and strain
		_vel[i] = 0.0,0.0,0.0;
		_gra[i] = 0.0,0.0,0.0;

		// Initialize internal values
		_mdl->InitIVS (names_vals, _vel[i], _gra[i], _ivs[i]);
	}

	// Initialize internal state
	if (IsActive) _init_internal_state();
}

inline void DiffusionElem::_B_mat(Mat_t const & dN, Mat_t const & J, Mat_t & B) const
{
	B = inv(J)*dN;
}

inline void DiffusionElem::_dist_to_face_nodes(Str_t Key, double const FaceValue, Array<Node*> const & FConn)
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

inline void DiffusionElem::_init_internal_state()
{
	// Allocate (local/element) internal force vector
	Vec_t f(_ge->NNodes);
	f.SetValues (0.0);

	// Calculate internal force vector;
	Mat_t dN;  // Shape Derivs
	Mat_t J;   // Jacobian matrix
	Mat_t B;   // grad-head matrix
	Vec_t vel; // Velocity
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_ge->Derivs   (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, dN);
		_ge->Jacobian (dN, J);
		_B_mat        (dN, J, B);
		DiffusionModel::Vec3ToVec (_gi, _vel[i], vel);
		f += trn(B)*vel*det(J)*_ge->IPs[i].w;
	}

	// Assemble (local/element) displacements vector.
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		_ge->Conn[i]->DOFVar(UD[j]).NaturalVal += f(i*_nd+j);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new 3D Diffusion element:
ProbElem * DiffusionMaker() { return new DiffusionElem; }

// Register element
int DiffusionRegister() { ProbElemFactory["Diffusion"]=DiffusionMaker;  return 0; }

// Call register
int __Diffusion_dummy_int  = DiffusionRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_DIFFUSIONELEM_H
