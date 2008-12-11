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
#include "fem/probelem.h"
#include "models/equilibmodel.h"
#include "util/string.h"
#include "util/util.h"
#include "util/lineparser.h"
#include "util/numstreams.h"
#include "linalg/laexpr.h"

using Util::SQ2;
using Util::_12_6;

namespace FEM
{

class EquilibElem : public ProbElem
{
public:
	// Constants
	static const char   UD [2][6][3];  ///< Essential DOF names == UD[_ge->NDim-1][iDOF]
	static const char   FD [2][6][3];  ///< Natural DOF names
	static const size_t NL [4];        ///< Number of additional labels (exceeding ND)
	static const char   LB [4][18][4]; ///< Additional labels

	// Constructor
	EquilibElem () : _gam(0.0), _di(-1), _gi(-1) {}

	// Destructor
	~EquilibElem () { for (size_t i=0; i<_mdl.Size(); ++i) delete _mdl[i]; }

	// Methods related to PROBLEM
	void    AddVolForces ();
	void    ClearDisp    ();
	void    SetActive    (bool Activate);
	void    EdgeBry      (Str_t Key, double Val, int iEdge);
	void    FaceBry      (Str_t Key, double Val, int iFace);
	void    CalcDeps     () const;
	Str_t   ModelName    () const { return (_mdl.Size()>0 ? _mdl[0]->Name() : "__no_model__"); }
	double  Val          (int iNod, Str_t Key) const;
	double  Val          (          Str_t Key) const;
	bool    IsEssen      (Str_t Key) const;
	void    SetProps     (Str_t Properties);
	void    SetModel     (Str_t ModelName, Str_t Prms, Str_t Inis);
	void    SetConn      (int iNod, FEM::Node * ptNode);
	void    Update       (double h, Vec_t const & dU, Vec_t & dFint);
	void    Backup       ();
	void    Restore      ();
	void    GetLbls      (Array<String> & Lbls) const;
	void    OutInfo      (std::ostream & os) const;
	size_t  NCMats       () const { return 1; }
	void    CMatrix      (size_t Idx, Mat_t & M) const;
	void    CMatMap      (size_t Idx, Array<size_t> & RMap, Array<size_t> & CMap, Array<bool> & RUPresc, Array<bool> & CUPresc) const;

	// Derived methods
	bool CheckModel () const;

private:
	// Data
	Array<EquilibModel*> _mdl; ///< Array of pointers to constitutive models
	double               _gam; ///< Specific weigth
	int                  _di;  ///< Dimension index == _ge->NDim-2
	int                  _gi;  ///< Geometry index: 3D=0, PStrain=1, PStress=2, Axis=3

	// Derived methods
	void _initialize (Str_t Type);

	// Private methods
	void _excavate           ();                                                   ///< Excavate element
	void _B_mat              (Mat_t const & dN, Mat_t const & J, Mat_t & B) const; ///< Calculate B matrix
	void _initial_state      ();                                                   ///< Calculate initial internal state
	void _equations_map      (Array<size_t> & RMap, Array<size_t> & CMap, Array<bool> & RUPresc, Array<bool> & CUPresc) const;
	void _dist_to_face_nodes (char const * Key, double FaceValue, Array<Node*> const & FaceConnects) const;

}; // class EquilibElem

// UD[_ge->NDim-1][iDOF]                   2D                          3D
const char   EquilibElem::UD [2][3][3] = {{"ux","uy",""},  {"ux","uy","uz"}};
const char   EquilibElem::FD [2][3][3] = {{"fx","fy",""},  {"fx","fy","fz"}};

// LB[_geom-1][iLabel]
const size_t EquilibElem::NL [4]  = {18,    // 3D                  
                                     16,    // 2D (plane-strain)
                                     10,    // 2D (plane-stress)
                                     18 };  // 2D (axis-symmetric)
const char   EquilibElem::LB [4][18][4] = {
	{"Ex", "Ey", "Ez",  "Exy", "Eyz", "Ezx", "Sx", "Sy" , "Sz", "Sxy", "Syz", "Szx", "E1", "E2", "E3", "S1", "S2", "S3" }, // 3D
	{"Ex", "Ey", "Ez",  "Exy", "Sx" , "Sy" , "Sz", "Sxy", "E1", "E2" , "S1" , "S2" , "p" , "q" , "Ev", "Ed", ""  , ""   }, // 2D (plane-strain)
	{"Ex", "Ey", "Exy", "Sx" , "Sy" , "Sxy", "E1", "E2" , "S1", "S2" , ""   , ""   , ""  , ""  , ""  , ""  , ""  , ""   }, // 2D (plane-stress)
	{"Ex", "Ey", "Ez",  "Exy", "Eyz", "Ezx", "Sx", "Sy" , "Sz", "Sxy", "Syz", "Szx", "E1", "E2", "E3", "S1", "S2", "S3" }  // 2D (axis-symmetric)
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void EquilibElem::AddVolForces() 
{
	// Verify if element is active
	if (_is_active==false) return;
	
	// Allocate (local/element) external volume force vector
	Vec_t fvol(_ge->NNodes);
	fvol.SetValues(0.0);

	// Allocate entities used for every integration point
	Vec_t N;
	Mat_t dN;
	Mat_t J;

	// Calculate local external volume force
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		double r = _ge->IPs[i].r;
		double s = _ge->IPs[i].s;
		double t = _ge->IPs[i].t;
		double w = _ge->IPs[i].w;

		Shape    (r,s,t, N);   // Calculate N functions for i IP
		Derivs   (r,s,t, dN);  // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (dN, J);      // Calculate J (Jacobian) matrix for i Integration Point

		fvol += -N*_gam*det(J)*w;
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
	if (_is_active==false) return;

	// Clear displacements
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		_ge->Conn[i]->DOFVar(UD[_d][j]).EssentialVal = 0.0;

	// Clear strains
	for (size_t i=0; i<_mdl.Size(); ++i) _mdl[i]->ClearStrain();
}

inline void EquilibElem::SetActive(bool Activate)
{
	if (_is_active==false && Activate)
	{
		// Set active
		_is_active = true;

		for (size_t i=0; i<_ge->Conn.Size(); ++i)
		{
			// Add Degree of Freedom to a node (Essential, Natural)
			for (int j=0; j<_nd; ++j) _ge->Conn[i]->AddDOF (UD[_d][j], FD[_d][j]);

			// Set SharedBy
			_ge->Conn[i]->SetSharedBy (_my_id);
		}

		// Apply body forces
		ApplyBodyForces ();
	}
	if (_is_active && Activate==false)
	{
		// Set active
		_is_active = false;

		for (size_t i=0; i<_ge->Conn.Size(); ++i)
		{
			// Remove SharedBy
			_ge->Conn[i]->RemoveSharedBy (_my_id); 

			// Remove Degree of Freedom to a node (Essential)
			if (_ge->Conn[i]->nSharedBy()==0) 
				for (int j=0; j<_nd; ++j) _ge->Conn[i]->RemoveDOF (UD[_d][j]);
		}

		// Apply surface tensions to perform excavation
		_excavate();

	}
}

inline Element * Element::EdgeBry(char const * Key, double Value, int iEdge)
{
	if (NDim<3) // For 1D/2D meshes, edges correspond to faces
	{

		// Skip if key is "Qb" TODO: Move to EquilibElem
		if (strcmp(Key,"Qb")==0) return this;

		Array<Node*> fnodes;
		GetFaceNodes        (iEdge, fnodes);
		_dist_to_face_nodes (Key, Value, fnodes);
	}
	else
	{
		throw new Fatal("Element::EdgeBry: Method not yet implemented for 3D meshes.");
	}
	return this;
}

inline Element * Element::FaceBry(char const * Key, double Value, int iFace)
{
	if (NDim==2) throw new Fatal("Element::FaceBry: This method must be called only for 3D meshes.");
	else
	{
		Array<Node*> fnodes;
		GetFaceNodes        (iFace, fnodes);
		_dist_to_face_nodes (Key, Value, fnodes);
	}
	return this;
}

inline void EquilibElem::CalcDeps() const
{
	if (_is_active==false) throw new Fatal("EquilibElem::CalcDepVars: This element is inactive (ID=%d, Tag=%d)",_my_id,_tag);
	if (_mdl.Size()==_ge->NIPs) for (size_t i=0; i<_ge->NIPs; i++) _mdl[i]->CalcDepVars();
	else throw new Fatal("EquilibElem::CalcDepVars: Constitutive models for this element (ID=%d, Tag=%d) were not set yet",_my_id,_tag);
}

inline double EquilibElem::Val(int iNod, char const * Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[_d][j])==0) return _ge->Conn[iNod]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[_d][j])==0) return _ge->Conn[iNod]->DOFVar(Name).NaturalVal;

	// Stress, strains, internal values, etc.
	Vec_t    ip_values (_ge->NIPs); // Vectors for extrapolation
	Vec_t nodal_values (_ge->NNodes);

	// Get integration point values
	if (_mdl.Size()==_ge->NIPs)
		for (size_t i=0; i<_ge->NIPs; i++)
			ip_values(i) = _mdl[i]->Val(Name);
	else throw new Fatal("EquilibElem::Val: Constitutive models for this element (ID==%d) were not set yet", _my_id);

	Extrapolate (ip_values, nodal_values);
	return nodal_values (iNod);
}

inline double EquilibElem::Val(char const * Name) const
{
	// Get integration point values
	double sum = 0.0;
	if (_mdl.Size()==_ge->NIPs)
		for (size_t i=0; i<_ge->NIPs; i++)
			sum += _mdl[i]->Val(Name);
	else throw new Fatal("EquilibElem::Val: Constitutive models for this element (ID==%d) were not set yet", _my_id);

	// Output single value at CG
	return sum/_ge->NIPs;
}

inline bool EquilibElem::IsEssen(char const * Name) const
{
	for (int i=0; i<_nd; ++i) if (strcmp(Name,UD[_d][i])==0) return true;
	return false;
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

inline void EquilibElem::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// Check _ge->NDim
	if (_ge->NDim<1) throw new Fatal("EquilibElem::SetModel: The space dimension (SetDim) must be set before calling this method");
	if (CheckConnect()==false) throw new Fatal("EquilibElem::SetModel: Connectivity is not correct. Connectivity MUST be set before calling this method");

	// If pointers to model was not already defined => No model was allocated
	if (_mdl.Size()==0)
	{
		_mdl.Resize(_ge->NIPs);
		for (size_t i=0; i<_ge->NIPs; ++i)
		{
			_mdl[i] = static_cast<EquilibModel*>(AllocModel(ModelName));
			_mdl[i]->SetGeom (_geom());
			_mdl[i]->SetPrms (Prms);
			_mdl[i]->SetInis (Inis);
		}
		if (_is_active) _initial_state ();
	}
	else throw new Fatal("EquilibElem::SetModel: Feature not implemented.");
}

inline void EquilibElem::SetConn(int iNod, FEM::Node * ptNode)
{
	// Check
	if (_ge->NNodes<1)               throw new Fatal("EquilibElem::Connect: __Internal Error__: There is a problem with the number of nodes: maybe derived elemet did not set _ge->NNodes");
	if (_ge->Conn.Size()<1)       throw new Fatal("EquilibElem::Connect: __Internal Error__: There is a problem with connectivity array: maybe derived elemet did not allocate _connect");
	if (_ge->NDim<0 || _d<0 || _nd<0) throw new Fatal("EquilibElem::Connect: __Internal Error__: There is a problem with _ge->NDim=%d, _d=%d, or _nd=%d\n (_ge->NDim=space dimension, _d=dimension index==_ge->NDim-1, and _nd=number of degrees of freedom)",_ge->NDim,_d,_nd);

	// Connectivity
	_ge->Conn[iNod] = ptNode;

	if (_is_active)
	{
		// Add Degree of Freedom to a node (Essential, Natural)
		for (int i=0; i<_nd; ++i) _ge->Conn[iNod]->AddDOF (UD[_d][i], FD[_d][i]);

		// Set shared
		_ge->Conn[iNod]->SetSharedBy(_my_id);
	}

	return this;
}

inline void EquilibElem::Update(double h, Vec_t const & dU, Vec_t & dFint)
{
	// Allocate (local/element) displacements vector
	Vec_t du(_nd*_ge->NNodes); // Delta disp. of this element

	// Assemble (local/element) displacements vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		du(i*_nd+j) = dU(_ge->Conn[i]->DOFVar(UD[_d][j]).EqID);

	// Allocate (local/element) internal force vector
	Vec_t df(_nd*_ge->NNodes); // Delta internal force of this element
	df.SetValues(0.0);

	// Allocate entities used for every integration point
	Mat_t dN;  // size = NumLocalCoords(ex.: r,s,t) x _ge->NNodes
	Mat_t J;       // Jacobian matrix
	Mat_t B;       // strain-displacement matrix
	Vec_t deps;    // delta strain vector
	Vec_t dsig;    // delta stress vector

	// Update model and calculate internal force vector;
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		double r = _ge->IPs[i].r;
		double s = _ge->IPs[i].s;
		double t = _ge->IPs[i].t; // only for 3D cases
		double w = _ge->IPs[i].w;

		Derivs   (r,s,t, dN);  // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (dN, J);      // Calculate J (Jacobian) matrix for i Integration Point
		_Bmat (dN, J, B);   // Calculate B matrix for i Integration Point

		deps = B*du;
		_mdl[i]->StateUpdate(deps, dsig);
		df += trn(B)*dsig*det(J)*w;
	}

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		dFint(_ge->Conn[i]->DOFVar(UD[_d][j]).EqID) += df(i*_nd+j);
}

inline void EquilibElem::Backup()
{
	for (size_t i=0; i<_ge->NIPs; ++i) _mdl[i]->BackupState();
}

inline void EquilibElem::Restore()
{
	for (size_t i=0; i<_ge->NIPs; ++i) _mdl[i]->RestoreState();
}

inline void EquilibElem::GetLbls(Array<String> & Lbls) const
{
	const int nl = 2*_nd+_nl; // total number of labels
	Lbls.Resize(nl);
	size_t k = 0;
	for (int i=0; i<_nd; ++i)
	{
		Lbls[k] = UD[_d][i];  k++;
		Lbls[k] = FD[_d][i];  k++;
	}
	for (int i=0; i<_nl; ++i)
	{
		Lbls[k] = LB[_geom()-1][i];  k++;
	}
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

	// Resize Ke
	Ke.Resize(_nd*_ge->NNodes, _nd*_ge->NNodes);
	Ke.SetValues(0.0);

	// Allocate entities used for every integration point
	Mat_t dN; // size = NumLocalCoords(ex.: r,s,t) x _ge->NNodes
	Mat_t J;      // Jacobian matrix
	Mat_t B;      // strain-displacement matrix
	Mat_t D;      // Constitutive matrix

	// Calculate Tangent Stiffness
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		double r = _ge->IPs[i].r;
		double s = _ge->IPs[i].s;
		double t = _ge->IPs[i].t;
		double w = _ge->IPs[i].w;

		Derivs   (r,s,t, dN); // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (dN, J);     // Calculate J (Jacobian) matrix for i Integration Point
		_Bmat (dN,J, B);   // Calculate B matrix for i Integration Point

		_mdl[i]->TgStiffness(D);
		Ke += trn(B)*D*B*det(J)*w;
	}
}

inline void EquilibElem::CMatMap(size_t Idx, Array<size_t> & RMap, Array<size_t> & CMap, Array<bool> & RUPresc, Array<bool> & CUPresc) const
{
	_equations_map(RMap, CMap, RUPresc, CUPresc);
}

inline bool EquilibElem::CheckModel() const
{
	if (_mdl.Size()!=_ge->NIPs) return false;
	for (size_t i=0; i<_ge->NIPs; ++i) if (_mdl[i]==NULL) return false;
	return true;
}


/* private */
inline void EquilibElem::_initialize(Str_t Type)
{
	     if (strcmp(Type,"")       ==0) _pt = (_ge->NDim==2 ? PStrain_T : Eq3D_T);
	else if (strcmp(Type,"PStrain")==0) _pt = PStrain_T;
	else if (strcmp(Type,"PStress")==0) _pt = PStress_T;
	else if (strcmp(Type,"Axis")   ==0) _pt = Axis_T;
	else if (strcmp(Type,"Rod")    ==0) _st = Rod_T;
	else if (strcmp(Type,"Beam")   ==0) _st = Beam_T;
	else if (strcmp(Type,"Spring") ==0) _st = Spring_T;
	_d  = _ge->NDim-1;
	_nd = EquilibElem::ND[_d];
	size_t idx = 0;
	     if (_st==PStrain_T) idx = 0;
	else if (_st==Eq3D_T   ) idx = 1;
	else if (_st==PStress_T) idx = 2;
	else if (_st==Axis_T   ) idx = 3;
	else if (_st==Rod_T    ) idx = 4;
	else if (_st==Beam_T   ) idx = 5;
	else if (_st==Spring_T ) idx = 6;
	_nl = EquilibElem::NL[idx];
}

inline void EquilibElem::_excavate()
{
	// Verify if element is in the boundary of excavation
	bool in_boundary = false;
	for (size_t i=0; i<_ge->NNodes; i++)
		if (_ge->Conn[i]->nSharedBy()>0) 
		{
			in_boundary=true;
			break;
		}

	// Calculate internal forces for element on boundary
	if (in_boundary)
	{
		Vec_t F(_ge->NDim*_ge->NNodes); F.SetValues(0.0);
		Vec_t S(_ge->NDim*_ge->NNodes); S.SetValues(0.0);

		// Allocate entities used for every integration point
		Vec_t N;
		Mat_t dN; // size = NumLocalCoords(ex.: r,s,t) x _ge->NNodes
		Mat_t J;      // Jacobian matrix
		Mat_t B;      // strain-displacement matrix
		Vec_t sig;    // stress tensor

		for (size_t i=0; i<_ge->NIPs; ++i)
		{
			// Temporary Integration Points
			double r = _ge->IPs[i].r;
			double s = _ge->IPs[i].s;
			double t = _ge->IPs[i].t;
			double w = _ge->IPs[i].w;

			Shape(r,s,t, N);       // Calculate N functions for i IP
			Derivs(r,s,t, dN);     // Calculate Derivatives of Shape functions w.r.t local coordinate system
			Jacobian(dN, J);       // Calculate J (Jacobian) matrix for i Integration Point
			_Bmat(dN,J, B);        // Calculate B matrix for i Integration Point

			// Mount S
			for (size_t j=0; j<_ge->NNodes; ++j) S(_ge->NDim*j+_ge->NDim-1)=N(j);

			// Get tensor for accumulate stress
			_mdl[i]->Sig(sig);

			// Calculate internal force vector
			F += trn(B)*sig*det(J)*w + S*_gam*det(J)*w;
		}	

		// Adding to boundary values of node if nSharedBy()>0
		for (size_t i=0; i<_ge->NNodes; ++i)
			if (_ge->Conn[i]->nSharedBy()>0) 
			{
				// Apply forces (only for equilibrium type neighbor elements)
							  _ge->Conn[i]->Bry("fx",F(_ge->NDim*i+0));
							  _ge->Conn[i]->Bry("fy",F(_ge->NDim*i+1));
				if (_ge->NDim==3) _ge->Conn[i]->Bry("fz",F(_ge->NDim*i+2));
			}
	}

}

inline void EquilibElem::_B_mat(Mat_t const & dN, Mat_t const & J, Mat_t & B) const
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
	Mat_t dN;
	dN = inv(J)*dN;

	// geometry type: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	switch (_geom())
	{
		case 2: // 2D(plane-strain)
		{
			const int n_scomps = 4; // number of stress compoments
			B.Resize (n_scomps,_nd*_ge->NNodes);
			for (size_t i=0; i<_ge->NNodes; ++i) // i row of B
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
			B.Resize (n_scomps,_nd*_ge->NNodes);
			for (size_t i=0; i<_ge->NNodes; ++i) // i row of B
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
			B.Resize(n_scomps,_nd*_ge->NNodes);
			for (size_t i=0; i<_ge->NNodes; ++i) // i row of B
			{
				B(0,0+i*_nd) =      dN(0,i);   B(0,1+i*_nd) =         0.0;
				B(1,0+i*_nd) =          0.0;   B(1,1+i*_nd) =     dN(1,i);
				B(2,0+i*_nd) =  dN(1,i)/SQ2;   B(2,1+i*_nd) = dN(0,i)/SQ2; // SQ2 => Mandel representation
			}
			return;
		}
		case 1: // 1D
		case 4: // 2D(axis-symmetric)
		default: throw new Fatal("EquilibElem::_Bmat: _Bmat() method is not available for GeometryType==%d",_geom());
	}
}

inline void EquilibElem::_initial_state()
{
	// Allocate (local/element) internal force vector
	Vec_t f(_nd*_ge->NNodes);
	f.SetValues(0.0);

	// Allocate entities used for every integration point
	Mat_t dN;  // size = NumLocalCoords(ex.: r,s,t) x _ge->NNodes
	Mat_t J;       // Jacobian matrix
	Mat_t B;       // strain-displacement matrix
	Vec_t sig;     // Stress vector in Mandel's notation

	// Calculate internal force vector;
	for (size_t i=0; i<_ge->NIPs; ++i)
	{
		double r = _ge->IPs[i].r;
		double s = _ge->IPs[i].s;
		double t = _ge->IPs[i].t; // only for 3D
		double w = _ge->IPs[i].w;

		Derivs   (r,s,t, dN); // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (dN, J);     // Calculate J (Jacobian) matrix for i Integration Point
		_Bmat (dN, J, B);  // Calculate B matrix for i Integration Point

		_mdl[i]->Sig(sig);
		f += trn(B)*sig*det(J)*w;
	}

	// Assemble (local/element) displacements vector.
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		_ge->Conn[i]->DOFVar(UD[_d][j]).NaturalVal += f(i*_nd+j);
}

inline void EquilibElem::_equations_map(Array<size_t> & RMap, Array<size_t> & CMap, Array<bool> & RUPresc, Array<bool> & CUPresc) const
{
	// Size of Stiffness/Mass Matrix
	int n_rows = _nd*_ge->NNodes; // == n_cols

	// Mounting a map of positions from Me to Global
	RMap       .Resize(n_rows);
	RUPresc.Resize(n_rows);

	int p = 0; // position inside matrix
	for (size_t i=0; i<_ge->NNodes; ++i)
	{
		for (int j=0; j<_nd; ++j)
		{
			RMap        [p] = _ge->Conn[i]->DOFVar(UD[_d][j]).EqID;
			RUPresc [p] = _ge->Conn[i]->DOFVar(UD[_d][j]).IsEssenPresc;
			p++;
		}
	}
	CMap        = RMap;
	CUPresc = RUPresc;
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
	Mat_t values;  values.Resize(_ge->NFNodes, _ge->NDim);  values.SetValues(0.0);
	Mat_t J;                         // Jacobian matrix. size = [1,2] x 3
	Vec_t face_shape(_ge->NFNodes); // Shape functions of a face/edge. size = _ge->NFNodes
	Mat_t F;                         // Shape function matrix
	Vec_t P;                         // Vector perpendicular to the face 
	for (size_t i=0; i<_ge->NFIps; i++)
	{
		double r = _ge->FIPs[i].r;
		double s = _ge->FIPs[i].s;
		double w = _ge->FIPs[i].w;
		FaceShape    (r, s, face_shape);
		F = trn(trn(face_shape)); // trick just to convert Vector face_shape to a col Matrix

		// Calculate perpendicular vector
		if (_ge->NDim==3)
		{
			FaceJacobian (FaceConnects, r, s, J);
			Vec_t V(3); V = J(0,0), J(0,1), J(0,2);
			Vec_t W(3); W = J(1,0), J(1,1), J(1,2);
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
	for (size_t i=0; i<_ge->NFNodes; ++i)
	{
		              FaceConnects[i]->Bry("fx",values(i,0));
		              FaceConnects[i]->Bry("fy",values(i,1));
		if (_ge->NDim==3) FaceConnects[i]->Bry("fz",values(i,2));
	}
}

}; // namespace FEM

#endif // MECHSYS_FEM_EQUILIB_H
