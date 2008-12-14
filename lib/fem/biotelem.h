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

#ifndef MECHSYS_FEM_BIOTELEM_H
#define MECHSYS_FEM_BIOTELEM_H

// MechSys
#include "fem/equilibelem.h"
#include "util/string.h"
#include "util/util.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "linalg/laexpr.h"

using Util::SQ2;
using Util::_12_6;

namespace FEM
{


class BiotElem : public EquilibElem
{
public:
	//{ Constants
	static const size_t ND_BIOT_3D;        ///< Number of DOF vars 3D
	static const char   UD_BIOT_3D[ 4][4]; ///< Essential DOF vars 3D
	static const char   FD_BIOT_3D[ 4][4]; ///< Natural DOF vars 3D
	static const size_t ND_BIOT_2D;        ///< Number of DOF vars 2D
	static const char   UD_BIOT_2D[ 3][4]; ///< Essential DOF vars 2D
	static const char   FD_BIOT_2D[ 3][4]; ///< Natural DOF vars 2D
	static const size_t NL_BIOT_3D;        ///< Number of labels 3D
	static const char   LB_BIOT_3D[20][4]; ///< Name of labels 3D
	static const size_t NL_BIOT_2D;        ///< Number of labels 2D
	static const char   LB_BIOT_2D[15][4]; ///< Name of labels 2D
	static const char   BIOT_PROP [ 2][8]; ///< Properties
	//}
	
	// Methods related to PROBLEM
	int         InitCtes     (int nDim);
	int         NProps       () const { return 2; } ///< "gam" and "gw"
	ProName_t * Props        () const { return BIOT_PROP; }
	void        ClearDisp    ();
	void        CalcDeps     () const;
	double      Val          (int iNod, Str_t Key) const;
	double      Val          (          Str_t Key) const;
	void        Update       (double h, Vec_t const & dU, Vec_t & dFint);
	void        OutInfo      (std::ostream & os) const;
	size_t      NCMats       () const { return 3; }
	size_t      NHMats       () const { return 1; }
	size_t      NUVecs       () const { return 2; }
	void        CMatrix      (size_t Idx, Mat_t & M) const;
	void        HMatrix      (size_t Idx, Mat_t & M) const;
	void        UVector      (size_t Idx, Vec_t & V) const;
	void        CMatMap      (size_t Idx, Array<size_t> & RMap, Array<size_t> & CMap, Array<bool> & RUPresc, Array<bool> & CUPresc) const;
	void        HMatMap      (size_t Idx, Array<size_t> & RMap, Array<size_t> & CMap, Array<bool> & RUPresc, Array<bool> & CUPresc) const;
	void        UVecMap      (size_t Idx, Array<size_t> & RMap) const;

private:
	// Derived methods
	void _excavate () { throw new Fatal("BiotElem::_excavate: Method not available"); }

	// Private methods
	void   _Bp_mat     (Mat_t const & dN, Mat_t const & J, Mat_t & Bp) const { Bp = inv(J)*dN; }
	void   _equi_map   (Array<size_t> & RMap, Array<bool> & RUPresc) const;
	void   _flow_map   (Array<size_t> & RMap, Array<bool> & RUPresc) const;
	void   _compute_K  (Mat_t & Ke) const;
	void   _compute_C  (Mat_t & Ce) const;
	void   _compute_L  (Mat_t & Le) const;
	void   _compute_H  (Mat_t & He) const;
	void   _compute_Qb (Vec_t & Qb) const;
	double _val_ip     (size_t iIP, Str_t Name) const;

}; // class BiotElem

//{ Constants
const size_t BiotElem::ND_BIOT_3D        = 4;
const char   BiotElem::UD_BIOT_3D[ 4][4] = {"ux", "uy", "uz", "pwp"};
const char   BiotElem::FD_BIOT_3D[ 4][4] = {"fx", "fy", "fz", "vol"};
const size_t BiotElem::ND_BIOT_2D        = 3;
const char   BiotElem::UD_BIOT_2D[ 3][4] = {"ux", "uy", "pwp"};
const char   BiotElem::FD_BIOT_2D[ 3][4] = {"fx", "fy", "vol"};
const size_t BiotElem::NL_BIOT_3D        = 20;
const char   BiotElem::LB_BIOT_3D[20][4] = {"Ex", "Ey", "Ez", "Exy", "Eyz", "Ezx", "Sx", "Sy", "Sz", "Sxy", "Syz", "Szx", "p", "q", "Ev", "Ed", "Vx", "Vy", "Vz", "H"};
const size_t BiotElem::NL_BIOT_2D        = 15;
const char   BiotElem::LB_BIOT_2D[15][4] = {"Ex", "Ey", "Ez", "Exy", "Sx", "Sy", "Sz", "Sxy", "p", "q", "Ev", "Ed", "Vx", "Vy", "H"};
const char   BiotElem::BIOT_PROP [ 2][8] = {"gam", "gw"};
//}


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline int BiotElem::InitCtes(int nDim)
{
	// Check
	if (nDim<2 || nDim>3) throw new Fatal("BiotElem::InitCtes: The space dimension must be 2 or 3. nDim==%d is invalid",nDim);

	// Set number of DOFs (_nd), number of labels (_nl), and arrays of essential UD, natural FD, and labels LB
	_gi = (nDim==2 ? 1 : 0);
	if (_gi==0)  // 3D
	{
		_nd = ND_BIOT_3D;
		UD  = UD_BIOT_3D;
		FD  = FD_BIOT_3D;
		_nl = NL_BIOT_3D; 
		LB  = LB_BIOT_3D;
	}
	else if (_gi==1)  // PlaneStrain
	{
		_nd = ND_BIOT_2D;
		UD  = UD_BIOT_2D;
		FD  = FD_BIOT_2D;
		_nl = NL_BIOT_2D; 
		LB  = LB_BIOT_2D;
	}
	else throw new Fatal("BiotElem::InitCtes: GeometryIndex _gi==%d is invalid",_gi);

	// Return geometry index
	return _gi;
}

inline void BiotElem::ClearDisp()
{
	if (IsActive==false) return;

	// Clear displacements
	size_t nde = _nd-1; // nDOFs equilib
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (size_t j=0; j<nde;         ++j)
		_ge->Conn[i]->DOFVar(UD[j]).EssentialVal = 0.0;

	// Clear strains
	for (size_t i=0; i<_eps.Size(); ++i) _eps[i] = 0.0,0.0,0.0, 0.0,0.0,0.0;
}

inline void BiotElem::CalcDeps() const
{
	if (IsActive==false) throw new Fatal("BiotElem::CalcDepVars: This element is inactive");
}

inline double BiotElem::Val(int iNod, Str_t Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[j])==0) return _ge->Conn[iNod]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[j])==0) return _ge->Conn[iNod]->DOFVar(Name).NaturalVal;

	// Stress, strains, internal values, etc.
	Vec_t    ip_values (_ge->NIPs); // Vectors for extrapolation
	Vec_t nodal_values (_ge->NNodes);

	// Get integration point values
	for (size_t i=0; i<_ge->NIPs; i++) ip_values(i) = _val_ip(i,Name);

	// Extrapolate
	_ge->Extrap (ip_values, nodal_values);
	return nodal_values (iNod);
}

inline double BiotElem::Val(Str_t Name) const
{
	// Get integration point values
	double sum = 0.0;
	for (size_t i=0; i<_ge->NIPs; i++) sum += _val_ip(i,Name);

	// Output single value at CG
	return sum/_ge->NIPs;
}

inline void BiotElem::Update(double h, Vec_t const & dU, Vec_t & dFint)
{
	// nDOFs
	size_t nde = _nd-1; // nDOFs equilib
	size_t ndf =     1; // nDOFs flow

	// Allocate (local/element) displacements vector
	Vec_t du(nde*_ge->NNodes);
	Vec_t dp(ndf*_ge->NNodes);
	Vec_t  p(ndf*_ge->NNodes);

	// Assemble (local/element) displacements vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	{
		for (size_t j=0; j<nde; ++j) du(i*nde+j) = dU(_ge->Conn[i]->DOFVar(UD[j    ]).EqID);
		for (size_t j=0; j<ndf; ++j) dp(i*ndf+j) = dU(_ge->Conn[i]->DOFVar(UD[j+nde]).EqID);
		for (size_t j=0; j<ndf; ++j)  p(i*ndf+j) =    _ge->Conn[i]->DOFVar(UD[j+nde]).EssentialVal;
	}

	// Allocate (local/element) internal force vector
	Vec_t df  (nde*_ge->NNodes); // Delta internal force of this element
	Vec_t dvol(ndf*_ge->NNodes); // Delta internal volumes of this element
	Mat_t Ke;                    // Stiffness matrix
	Mat_t Ce;                    // Coupling matrix 1
	Mat_t Le;                    // Coupling matrix 2
	Mat_t He;                    // Permeability matrix
	CMatrix (0,Ke);
	CMatrix (1,Ce);
	CMatrix (2,Le);
	HMatrix (0,He);

	Vector<double> Qb;
	_compute_Qb(Qb);

	df   = Ke*du + Ce*dp;
	dvol = Le*du + h*He*p + h*Qb;

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
	}

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	{
		for (size_t j=0; j<nde; ++j) dFint(_ge->Conn[i]->DOFVar(UD[j    ]).EqID) += df  (i*nde+j);
		for (size_t j=0; j<ndf; ++j) dFint(_ge->Conn[i]->DOFVar(UD[j+nde]).EqID) += dvol(i*ndf+j);
	}
}

inline void BiotElem::OutInfo(std::ostream & os) const
{
	for (size_t i=0; i<_ge->NIPs; i++)
	{
		os << "IP # " << i << " Sx,Sy,Sz = " << _12_6 << _sig[i](0) << _12_6 << _sig[i](1) << _12_6 << _sig[i](2);
		os <<                "  Ex,Ey,Ez = " << _12_6 << _eps[i](0) << _12_6 << _eps[i](1) << _12_6 << _eps[i](2) << " ";
	}
}

inline void BiotElem::CMatrix(size_t Idx, Mat_t & M) const
{
	     if (Idx==0) _compute_K (M);
	else if (Idx==1) _compute_C (M);
	else if (Idx==2) _compute_L (M);
}

inline void BiotElem::HMatrix(size_t index, Mat_t & M) const
{
	_compute_H (M);
}

inline void BiotElem::UVector(size_t Idx, Vec_t & V) const
{
	if (Idx==0)
	{
		V.Resize(_ge->NNodes);
		size_t nde = _nd-1; // nDOFs equilib
		size_t ndf =     1; // nDOFs flow
		for (size_t i=0; i<_ge->NNodes; ++i)
		for (size_t j=0; j<ndf;      ++j)
			V(i*ndf+j) = _ge->Conn[i]->DOFVar (UD[j+nde]).EssentialVal;
		Mat_t H;
		_compute_H (H);
		V = H*V;
	}
	else if (Idx==1) _compute_Qb(V);
}

inline void BiotElem::CMatMap(size_t Idx, Array<size_t> & RMap, Array<size_t> & CMap, Array<bool> & RUPresc, Array<bool> & CUPresc) const
{
	if (Idx==0)
	{
		_equi_map (RMap, RUPresc);
		CMap    = RMap;
		CUPresc = RUPresc;
	}
	else if (Idx==1)
	{
		_equi_map (RMap, RUPresc);
		_flow_map (CMap, CUPresc);
	}
	else if (Idx==2)
	{
		_flow_map (RMap, RUPresc);
		_equi_map (CMap, CUPresc);
	}
	else throw new Fatal("BiotElem::CMatMap: Idx==%d is invalid",Idx);
}

inline void BiotElem::HMatMap(size_t Idx, Array<size_t> & RMap, Array<size_t> & CMap, Array<bool> & RUPresc, Array<bool> & CUPresc) const
{
	_flow_map (RMap, RUPresc);
	CMap    = RMap;
	CUPresc = RUPresc;
}

inline void BiotElem::UVecMap(size_t Idx, Array<size_t> & RMap) const
{
	Array<bool> rupresc;
	_flow_map (RMap, rupresc);
}


/* private */

inline void BiotElem::_equi_map(Array<size_t> & RMap, Array<bool> & RUPresc) const
{
	// Map of positions from Me to Global
	size_t nde = _nd-1; // nDOFs equilib
	RMap   .Resize(nde*_ge->NNodes);
	RUPresc.Resize(nde*_ge->NNodes);
	int p = 0; // position inside matrix
	for (size_t i=0; i<_ge->NNodes; ++i)
	{
		for (size_t j=0; j<nde; ++j)
		{
			RMap    [p] = _ge->Conn[i]->DOFVar (UD[j]).EqID;
			RUPresc [p] = _ge->Conn[i]->DOFVar (UD[j]).IsEssenPresc;
			p++;
		}
	}
}

inline void BiotElem::_flow_map(Array<size_t> & RMap, Array<bool> & RUPresc) const
{
	// Map of positions from Me to Global
	size_t nde = _nd-1; // nDOFs equilib
	size_t ndf =     1; // nDOFs flow
	RMap   .Resize(ndf*_ge->NNodes);
	RUPresc.Resize(ndf*_ge->NNodes);
	int p = 0; // position inside matrix
	for (size_t i=0; i<_ge->NNodes; ++i)
	{
		for (size_t j=0; j<ndf; ++j)
		{
			RMap    [p] = _ge->Conn[i]->DOFVar (UD[j+nde]).EqID;
			RUPresc [p] = _ge->Conn[i]->DOFVar (UD[j+nde]).IsEssenPresc;
			p++;
		}
	}
}

inline void BiotElem::_compute_K(Mat_t & Ke) const
{
	/* Stiffness:
	   ==========
	                 /    T
	        [Ke]  =  | [B]  * [D] * [B]  * dV
	                 /
	*/

	size_t nde = _nd-1; // nDOFs equilib
	Ke.Resize    (nde*_ge->NNodes, nde*_ge->NNodes);
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

inline void BiotElem::_compute_C(Mat_t & Ce) const
{
	/*  Coupling Matrix L1:
	    ===================
	
	                 /    T      T
	        [C]   =  | [B]  * {m} * {N}  * dV    // coupled saturated
	                 /    u            p
	
	    Note:   p => pore-pressure
	    =====   u => displacement
	
	               T
	            {m} = [ 1 1 1 0 0 0 ]
	*/
	
	size_t nde = _nd-1; // nDOFs equilib
	size_t ndf =     1; // nDOFs flow
	Ce.Resize    (nde*_ge->NNodes, ndf*_ge->NNodes); 
	Ce.SetValues (0.0);
	Vec_t N,m;
	Mat_t dN,J,B;
	     if (_ge->NDim==2) { m.Resize(4); m = 1.0,1.0,1.0, 0.0;         }
	else if (_ge->NDim==3) { m.Resize(6); m = 1.0,1.0,1.0, 0.0,0.0,0.0; }
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_ge->Shape    (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, N);
		_ge->Derivs   (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, dN);
		_ge->Jacobian (dN, J);
		_B_mat        (dN, J, B);
		Ce += trn(B)*m*trn(N)*det(J)*_ge->IPs[i].w;
	}
}

inline void BiotElem::_compute_L(Mat_t & Le) const
{
	Mat_t Ce;
	_compute_C (Ce);
	Le = trn(Ce);
}

inline void BiotElem::_compute_H(Mat_t & He) const
{
	/*  Permeability Matrix H:
	    ======================
	
	                   /    T
	         [H]   =  -| [Bp]  * [K] * [Bp]  * dV
	                   /
	
	    Note:   p => pore-pressure
	    =====   u => displacement
	
	
	               T
	            {m} = [ 1 1 1 0 0 0 ]
	*/
	
	double gw  = Prop("gw");
	size_t ndf = 1; // nDOFs flow
	He.Resize    (ndf*_ge->NNodes, ndf*_ge->NNodes);
	He.SetValues (0.0);
	Mat_t dN,J,Bp,K;
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_ge->Derivs          (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, dN);
		_ge->Jacobian        (dN, J);
		_Bp_mat              (dN, J, Bp);
		_mdl->TgPermeability (K);
		He += -trn(Bp)*K*Bp*det(J)*_ge->IPs[i].w/gw;
	}
}

inline void BiotElem::_compute_Qb(Vec_t & Qb) const
{
	/*	 Matrix Qb:
		 ============================
	
	                    1   /    T                   
	         [Qb]   =  ---  | [B]  * [K] * {b}  * dV
	                    gw  /    p            w  
	*/
	
	double gw  = Prop("gw");
	size_t ndf = 1; // nDOFs flow
	Qb.Resize    (ndf*_ge->NNodes); 
	Qb.SetValues (0.0);
	Vec_t b;
	Mat_t dN,J,Bp,K;
	     if (_ge->NDim==2) { b.Resize(2); b = 0.0, gw; }
	else if (_ge->NDim==3) { b.Resize(3); b = 0.0, 0.0, gw; }
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		_ge->Derivs          (_ge->IPs[i].r, _ge->IPs[i].s, _ge->IPs[i].t, dN);
		_ge->Jacobian        (dN, J);
		_Bp_mat              (dN, J, Bp);
		_mdl->TgPermeability (K);
		Qb += trn(Bp)*K*b*det(J)*_ge->IPs[i].w/gw;
	}
}

inline double BiotElem::_val_ip(size_t iIP, Str_t Name) const
{
	     if (strcmp(Name,"Sx" )==0)                          return _sig[iIP](0);
	else if (strcmp(Name,"Sy" )==0)                          return _sig[iIP](1);
	else if (strcmp(Name,"Sz" )==0)                          return _sig[iIP](2);
	else if (strcmp(Name,"Sxy")==0 || strcmp(Name,"Syx")==0) return _sig[iIP](3)/SQ2;
	else if (strcmp(Name,"Syz")==0 || strcmp(Name,"Szy")==0) return _sig[iIP](4)/SQ2;
	else if (strcmp(Name,"Szx")==0 || strcmp(Name,"Sxz")==0) return _sig[iIP](5)/SQ2;
	else if (strcmp(Name,"p"  )==0)                          return (_sig[iIP](0)+_sig[iIP](1)+_sig[iIP](2))/3.0;
	else if (strcmp(Name,"q"  )==0)                          return sqrt(((_sig[iIP](0)-_sig[iIP](1))*(_sig[iIP](0)-_sig[iIP](1)) + (_sig[iIP](1)-_sig[iIP](2))*(_sig[iIP](1)-_sig[iIP](2)) + (_sig[iIP](2)-_sig[iIP](0))*(_sig[iIP](2)-_sig[iIP](0)) + 3.0*(_sig[iIP](3)*_sig[iIP](3) + _sig[iIP](4)*_sig[iIP](4) + _sig[iIP](5)*_sig[iIP](5)))/2.0);
	else if (strcmp(Name,"Ex" )==0)                          return _eps[iIP](0);
	else if (strcmp(Name,"Ey" )==0)                          return _eps[iIP](1);
	else if (strcmp(Name,"Ez" )==0)                          return _eps[iIP](2);
	else if (strcmp(Name,"Exy")==0 || strcmp(Name,"Eyx")==0) return _eps[iIP](3)/SQ2;
	else if (strcmp(Name,"Eyz")==0 || strcmp(Name,"Ezy")==0) return _eps[iIP](4)/SQ2;
	else if (strcmp(Name,"Ezx")==0 || strcmp(Name,"Exz")==0) return _eps[iIP](5)/SQ2;
	else if (strcmp(Name,"Ev" )==0)                          return _eps[iIP](0)+_eps[iIP](1)+_eps[iIP](2); 
	else if (strcmp(Name,"Ed" )==0)                          return sqrt(2.0*((_eps[iIP](0)-_eps[iIP](1))*(_eps[iIP](0)-_eps[iIP](1)) + (_eps[iIP](1)-_eps[iIP](2))*(_eps[iIP](1)-_eps[iIP](2)) + (_eps[iIP](2)-_eps[iIP](0))*(_eps[iIP](2)-_eps[iIP](0)) + 3.0*(_eps[iIP](3)*_eps[iIP](3) + _eps[iIP](4)*_eps[iIP](4) + _eps[iIP](5)*_eps[iIP](5))))/3.0;

	// Flow velocities
	else if (strcmp(Name,"Vx" )==0) return 0.0;
	else if (strcmp(Name,"Vy" )==0) return 0.0;
	else if (strcmp(Name,"Vz" )==0) return 0.0;

	// Total head
	else if (strcmp(Name,"H" )==0) return 0.0;

	else throw new Fatal("BiotElem::_val_ip: Name==%s if not available for this element",Name);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new 3D Biot element:
ProbElem * BiotMaker() { return new BiotElem; }

// Register element
int BiotRegister() { ProbElemFactory["Biot"]=BiotMaker;  return 0; }

// Call register
int __Biot_dummy_int = BiotRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_BIOTELEM_H
