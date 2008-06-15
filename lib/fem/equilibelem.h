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


/* __ Equilibrium capable element __

ux = Nodal displacement increment in x direction
uy = Nodal displacement increment in y direction
uz = Nodal displacement increment in z direction
fx = Nodal force increment in x direction
fy = Nodal force increment in y direction
fz = Nodal force increment in z direction
tx = Traction increment in x direction on face
ty = Traction increment in y direction on face
tz = Traction increment in z direction on face

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
	// Destructor
	virtual ~EquilibElem() {}

	// Derived methods
	bool      IsReady         () const;
	bool      IsEssential     (char const * DOFName) const;
	void      SetModel        (char const * ModelName, char const * Prms, char const * Inis);
	Element * SetNode         (int iNodeLocal, int iNodeGlobal);
	void      UpdateState     (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
	void      BackupState     ();
	void      RestoreState    ();
	void      SetGeometryType (int Geom);  
	void      SetProperties   (Array<double> const & EleProps) { _unit_weight=EleProps[0]; }
	String    OutCenter       (bool PrintCaptionOnly) const;
	void      OutNodes        (LinAlg::Matrix<double> & Values, Array<String> & Labels) const;
	void      Deactivate      ();

	// Derived methods to assemble DAS matrices
	size_t nOrder1Matrices () const { return 1; }
	void   Order1MatMap    (size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const;
	void   Order1Matrix    (size_t Index, LinAlg::Matrix<double> & Ke) const; // Stiffness

	// Methods
	void B_Matrix (LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const;

	// Access methods
	double Val (int iNodeLocal, char const * Name) const; ///< Return computed values at the CG of the element. Ex.: Name="Sx", "Sxy", "Ex", etc.

private:
	// Data
	int                  _n_stress;
	Array<EquilibModel*> _a_model;
	double               _unit_weight;

	// Private methods
	void _calc_initial_internal_forces ();

}; // class EquilibElem


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */
	
// Derived methods

inline bool EquilibElem::IsReady() const
{
	if (_a_model.Size()==static_cast<size_t>(_n_int_pts) && _connects.Size()==static_cast<size_t>(_n_nodes)) return true;
	else return false;
}

inline void EquilibElem::SetGeometryType(int Geom)
{
	// Set geometry type: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	_geom = Geom;

	// Set the number of stresses associated with the geometry type
	switch (_geom)
	{
		case 2: { _n_stress = 4; return; } // 2D(plane-strain)
		case 3: { _n_stress = 6; return; } // 3D
		case 5: { _n_stress = 3; return; } // 2D(plane-stress)
		case 1: // 1D
		case 4: // 2D(axis-symmetric)
		default:
			throw new Fatal("EquilibElem::SetGeometryType: GeometryType==%d is not implemented yet",_geom);
	}
}

inline bool EquilibElem::IsEssential(char const * DOFName) const
{
	if (strcmp(DOFName,"ux")==0 || strcmp(DOFName,"uy")==0) return true;
	if (_n_dim==3               && strcmp(DOFName,"uz")==0) return true;
	return false;
}

inline void EquilibElem::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// If pointers to model was not already defined => No model was allocated
	if (_a_model.Size()==0)
	{
		// Resize the array of model pointers
		_a_model.Resize(_n_int_pts);

		// Loop along integration points
		for (int i=0; i<_n_int_pts; ++i)
		{
			// Allocate a new model and set parameters
			_a_model[i] = static_cast<EquilibModel*>(AllocModel(ModelName));
			_a_model[i]->SetGeom(_geom);
			_a_model[i]->SetPrms(Prms);
			_a_model[i]->SetInis(Inis);
		}

		// Calculate initial internal forces
		_calc_initial_internal_forces();
	}
	else throw new Fatal("EquilibElem::SetModel: Feature not implemented.");
}

inline Element * EquilibElem::SetNode(int iNodeLocal, int iNodeGlobal)
{
	// Connects
	_connects[iNodeLocal] = Nodes[iNodeGlobal];

	// Add Degree of Freedom to a node (Essential, Natural)
	Nodes[iNodeGlobal]->AddDOF("ux", "fx");
	Nodes[iNodeGlobal]->AddDOF("uy", "fy");
	if (_n_dim==3)
	Nodes[iNodeGlobal]->AddDOF("uz", "fz");

	// Shared
	Nodes[iNodeGlobal]->SetSharedBy(_my_id);

	return this;
}

inline void EquilibElem::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
{
	// Allocate (local/element) displacements vector
	LinAlg::Vector<double> dU(_n_dim*_n_nodes); // Delta disp. of this element
	
	// Assemble (local/element) displacements vector
	for (int i=0; i<_n_nodes; ++i)
	{
		dU(i*_n_dim  ) = dUglobal(_connects[i]->DOFVar("ux").EqID);
		dU(i*_n_dim+1) = dUglobal(_connects[i]->DOFVar("uy").EqID);
		if (_n_dim==3)
		dU(i*_n_dim+2) = dUglobal(_connects[i]->DOFVar("uz").EqID);
	}
	
	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> dF(_n_dim*_n_nodes); // Delta internal force of this element
	dF.SetValues(0.0);
	
	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs;  // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;       // Jacobian matrix
	LinAlg::Matrix<double> B;       // strain-displacement matrix
	LinAlg::Vector<double> DEps;    // Strain vector 
	LinAlg::Vector<double> DSig;    // Stress vector 

	// Loop along integration points
	for (int i_ip=0; i_ip<_n_int_pts; ++i_ip)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i_ip].r;
		double s = _a_int_pts[i_ip].s;
		double t = _a_int_pts[i_ip].t; // only for 3D cases
		double w = _a_int_pts[i_ip].w;

		Derivs   (r,s,0.0, derivs);  // Calculate Derivatives of Shape functions w.r.t local coordinate system
		if (_n_dim==3)
		Derivs   (r,s,t, derivs);  
		Jacobian (derivs, J);      // Calculate J (Jacobian) matrix for i_ip Integration Point
		B_Matrix (derivs, J, B);   // Calculate B matrix for i_ip Integration Point

		// Calculate a tensor for the increments of strain
		DEps = B*dU;
		
		// Update model
		_a_model[i_ip]->StressUpdate(DEps, DSig);

		// Calculate internal force vector;
		dF += trn(B)*DSig*det(J)*w;
	}

	// Return internal forces
	for (int i=0; i<_n_nodes; ++i)
	{
		// Sum up contribution to internal forces vector
		dFint(_connects[i]->DOFVar("fx").EqID) += dF(i*_n_dim  );
		dFint(_connects[i]->DOFVar("fy").EqID) += dF(i*_n_dim+1);
		if (_n_dim==3)
		dFint(_connects[i]->DOFVar("fz").EqID) += dF(i*_n_dim+2);
	}
}

inline void EquilibElem::BackupState()
{
	for (int i=0; i<_n_int_pts; ++i)
		_a_model[i]->BackupState();
}

inline void EquilibElem::RestoreState()
{
	for (int i=0; i<_n_int_pts; ++i)
		_a_model[i]->RestoreState();
}

inline String EquilibElem::OutCenter(bool PrintCaptionOnly=false) const
{
	// Auxiliar variables
	std::ostringstream oss;
	/*
	// Number of state values
	int n_int_state_vals = _a_model[0]->nInternalStateValues();
	
	// Print caption
	if (PrintCaptionOnly)
	{
		// Stress and strains
		oss << Util::_8s<< "p"  << Util::_8s<< "q"  << Util::_8s<< "sin3th" << Util::_8s<< "Ev"  << Util::_8s<< "Ed";
		oss << Util::_8s<< "Sx" << Util::_8s<< "Sy" << Util::_8s<< "Sz" << Util::_8s<< "Sxy" << Util::_8s<< "Syz" << Util::_8s<< "Sxz";
		oss << Util::_8s<< "Ex" << Util::_8s<< "Ey" << Util::_8s<< "Ez" << Util::_8s<< "Exy" << Util::_8s<< "Eyz" << Util::_8s<< "Exz";

		// Internal state values
		Array<String> str_state_names;   _a_model[0]->InternalStateNames(str_state_names);
		for (int i=0; i<n_int_state_vals; ++i)
			oss << Util::_8s<< str_state_names[i];
		oss << std::endl;
	}
	else
	{
		// Stress, strains and internal state values evaluated at the center of the element
		Tensors::Tensor2 sig_cen(0.0);
		Tensors::Tensor2 eps_cen(0.0);
		Array<double> int_state_vals_cen(n_int_state_vals);
		int_state_vals_cen = 0.0;

		// Loop over integration points
		for (int i_ip=0; i_ip<_n_int_pts; ++i_ip)
		{
			// Stress and strains
			sig_cen += _a_model[i_ip]->Sig();
			eps_cen += _a_model[i_ip]->Eps();

			// Internal state values
			Array<double> int_state_vals;    _a_model[i_ip]->InternalStateValues(int_state_vals);
			for (int j=0; j<n_int_state_vals; ++j)
				int_state_vals_cen[j] += int_state_vals[j];
		}
		
		// Average stress and strains
		sig_cen = sig_cen / _n_int_pts;
		eps_cen = eps_cen / _n_int_pts;
		
		// Average internal state values
		for (int j=0; j<n_int_state_vals; ++j)
			int_state_vals_cen[j] = int_state_vals_cen[j] / _n_int_pts;;

		// Calculate stress invariants
		double            p,q,sin3th;
		Tensors::Tensor2  S;
		Stress_p_q_S_sin3th(sig_cen,p,q,S,sin3th);

		// Calculate strain invariants
		double     Ev,Ed;
		Strain_Ev_Ed(eps_cen,Ev,Ed);

		// Output
		oss << Util::_8s<< p          << Util::_8s<< q          << Util::_8s<< sin3th << Util::_8s<< Ev*100.0       << Util::_8s<< Ed*100.0;
		oss << Util::_8s<< sig_cen(0) << Util::_8s<< sig_cen(1) << Util::_8s<< sig_cen(2)         << Util::_8s<< sig_cen(3)/SQ2 << Util::_8s<< sig_cen(4)/SQ2 << Util::_8s<< sig_cen(5)/SQ2;
		oss << Util::_8s<< eps_cen(0) << Util::_8s<< eps_cen(1) << Util::_8s<< eps_cen(2)         << Util::_8s<< eps_cen(3)/SQ2 << Util::_8s<< eps_cen(4)/SQ2 << Util::_8s<< eps_cen(5)/SQ2;
		for (int j=0; j<n_int_state_vals; ++j)
			oss << Util::_8s<< int_state_vals_cen[j];
		oss << std::endl;
	} */

	return oss.str(); 
}

inline void EquilibElem::OutNodes(LinAlg::Matrix<double> & Values, Array<String> & Labels) const
{
	if (_n_dim==2)
	{
		int const DATA_COMPS=10;
		Values.Resize(_n_nodes,DATA_COMPS);
		Labels.Resize(DATA_COMPS);
		Labels[ 0] = "ux"; Labels[ 1] = "uy";  
		Labels[ 2] = "fx"; Labels[ 3] = "fy"; 
		Labels[ 4] = "Ex"; Labels[ 5] = "Ey"; Labels[ 6] = "Exy";
		Labels[ 7] = "Sx"; Labels[ 8] = "Sy"; Labels[ 9] = "Sxy";

		for (int i_node=0; i_node<_n_nodes; i_node++)
		{
			Values(i_node,0) = _connects[i_node]->DOFVar("ux").EssentialVal;
			Values(i_node,1) = _connects[i_node]->DOFVar("uy").EssentialVal;
			Values(i_node,2) = _connects[i_node]->DOFVar("fx").NaturalVal;
			Values(i_node,3) = _connects[i_node]->DOFVar("fy").NaturalVal;
		}

		//Extrapolation
		LinAlg::Vector<double> ip_values(_n_int_pts);
		LinAlg::Vector<double> nodal_values(_n_nodes);

		// Stress & Strain
		for (int i_comp=0; i_comp<6; i_comp++) // six values: stress + strain
		{
			// Get integration point values
			for (int i_ip=0; i_ip<_n_int_pts; i_ip++)
				ip_values(i_ip) = _a_model[i_ip]->Val(Labels[i_comp+4].c_str());

			Extrapolate(ip_values, nodal_values);

			// Place nodal values
			for (int j_node=0; j_node<_n_nodes; j_node++)
				Values(j_node,i_comp+4) = nodal_values(j_node);
		}
	}
	else // _n_dim==3
	{
		int const DATA_COMPS=18;
		Values.Resize(_n_nodes,DATA_COMPS);
		Labels.Resize(DATA_COMPS);
		Labels[ 0] = "ux" ; Labels[ 1] = "uy" ; Labels[ 2] = "uz"; 
		Labels[ 3] = "fx" ; Labels[ 4] = "fy" ; Labels[ 5] = "fz";
		Labels[ 6] = "Ex"; Labels[ 7] = "Ey"; Labels[ 8] = "Ez"; Labels[ 9] = "Exy"; Labels[10] = "Eyz"; Labels[11] = "Exz";
		Labels[12] = "Sx"; Labels[13] = "Sy"; Labels[14] = "Sz"; Labels[15] = "Sxy"; Labels[16] = "Syz"; Labels[17] = "Sxz";

		for (int i_node=0; i_node<_n_nodes; i_node++)
		{
			Values(i_node,0) = _connects[i_node]->DOFVar("ux").EssentialVal;
			Values(i_node,1) = _connects[i_node]->DOFVar("uy").EssentialVal;
			Values(i_node,2) = _connects[i_node]->DOFVar("uz").EssentialVal;
			Values(i_node,3) = _connects[i_node]->DOFVar("fx").NaturalVal;
			Values(i_node,4) = _connects[i_node]->DOFVar("fy").NaturalVal;
			Values(i_node,5) = _connects[i_node]->DOFVar("fz").NaturalVal;
		}
		//Extrapolation
		LinAlg::Vector<double> ip_values(_n_int_pts);
		LinAlg::Vector<double> nodal_values(_n_nodes);
		
		// Stress & Strain
		for (int i_comp=0; i_comp<12; i_comp++) // six values: stress + strain
		{
			// Get integration point values
			for (int i_ip=0; i_ip<_n_int_pts; i_ip++)
				ip_values(i_ip) = _a_model[i_ip]->Val(Labels[i_comp+6].c_str());

			Extrapolate(ip_values, nodal_values);

			// Place nodal values
			for (int j_node=0; j_node<_n_nodes; j_node++)
				Values(j_node,i_comp+6) = nodal_values(j_node);
		}
	}
}

inline double EquilibElem::Val(int iNodeLocal, char const * Name) const
{
	LinAlg::Vector<double> ip_values(_n_int_pts);
	LinAlg::Vector<double> nodal_values(_n_nodes);

	// Get integration point values
	for (int i_ip=0; i_ip<_n_int_pts; i_ip++)
		ip_values(i_ip) = _a_model[i_ip]->Val(Name);

	// Extrapolation
	Extrapolate (ip_values, nodal_values);

	return nodal_values(iNodeLocal);
}

inline void EquilibElem::Deactivate()
{
	throw new Fatal("EquilibElem::Deactivate: Feature not implemented yet");
}

// Derived methods to assemble DAS matrices

inline void EquilibElem::Order1MatMap(size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const
{
	// Size of Ke
	int n_rows = _n_dim*_n_nodes; // == n_cols

	// Mounting a map of positions from Ke to Global
	int idx_Ke = 0;                // position (idx) inside Ke matrix
	RowsMap       .Resize(n_rows); // size=Ke.Rows()=Ke.Cols()
	RowsEssenPresc.Resize(n_rows); // size=Ke.Rows()=Ke.Cols()

	// Fill map of Ke position to K position of DOFs components
	for (int i_node=0; i_node<_n_nodes; ++i_node)
	{
		RowsMap        [idx_Ke] = _connects[i_node]->DOFVar("ux").EqID; 
		RowsEssenPresc [idx_Ke] = _connects[i_node]->DOFVar("ux").IsEssenPresc; 
		idx_Ke++;
		RowsMap        [idx_Ke] = _connects[i_node]->DOFVar("uy").EqID; 
		RowsEssenPresc [idx_Ke] = _connects[i_node]->DOFVar("uy").IsEssenPresc; 
		idx_Ke++;
		if (_n_dim==3)
		{
			RowsMap        [idx_Ke] = _connects[i_node]->DOFVar("uz").EqID; 
			RowsEssenPresc [idx_Ke] = _connects[i_node]->DOFVar("uz").IsEssenPresc; 
			idx_Ke++;
		}
	}
	ColsMap        = RowsMap;
	ColsEssenPresc = RowsEssenPresc;
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
	Ke.Resize(_n_dim*_n_nodes, _n_dim*_n_nodes); // sum(Bt*D*B*det(J)*w)
	Ke.SetValues(0.0);

	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs; // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;      // Jacobian matrix
	LinAlg::Matrix<double> B;      // strain-displacement matrix
	LinAlg::Matrix<double> D;      // Constitutive matrix

	// Loop along integration points
	for (int i_ip=0; i_ip<_n_int_pts; ++i_ip)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i_ip].r;
		double s = _a_int_pts[i_ip].s;
		double t = _a_int_pts[i_ip].t;
		double w = _a_int_pts[i_ip].w;

		Derivs(r,s,0.0, derivs);  // Calculate Derivatives of Shape functions w.r.t local coordinate system
		if (_n_dim==3)
		Derivs   (r,s,t, derivs);  
		Jacobian(derivs, J);          // Calculate J (Jacobian) matrix for i_ip Integration Point
		B_Matrix(derivs,J, B);        // Calculate B matrix for i_ip Integration Point

		// Constitutive tensor 
		_a_model[i_ip]->TgStiffness(D); 

		// Calculate Tangent Stiffness
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
	
	// Resize B matrix
	B.Resize(_n_stress,_n_dim*_n_nodes);

	// Inverse of Jacobian
	LinAlg::Matrix<double> inv_J(J.Rows(),J.Cols());
	J.Inv(inv_J);

	// Cartesian derivatives
	LinAlg::Matrix<double> cart_derivs;
	cart_derivs = inv(J)*derivs;

	// geometry type: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	switch (_geom)
	{
		case 2: // 2D(plane-strain)
		{
			// Loop along all nodes of the element
			double dNdX,dNdY;
			int  j=0; // j column of B
			for (int i=0; i<_n_nodes; ++i) // i row of B
			{
				// Assemble B matrix
				j = i*_n_dim;
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
			// Loop along all nodes of the element
			double dNdX,dNdY,dNdZ;
			int  j=0; // j column of B
			for (int i=0; i<_n_nodes; ++i) // i row of B
			{
				// Assemble B matrix
				j = i*_n_dim;
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
			// Loop along all nodes of the element
			double dNdX,dNdY;
			int  j=0; // j column of B
			for (int i=0; i<_n_nodes; ++i) // i row of B
			{
				// Assemble B matrix
				j = i*_n_dim;
				dNdX=-cart_derivs(0,i);  dNdY=-cart_derivs(1,i);  // Negative values => Soil mechanics convention
				B(0,0+j) =      dNdX;   B(0,1+j) =  0.0;
				B(1,0+j) =       0.0;   B(1,1+j) = dNdY;
				B(2,0+j) =  dNdY/SQ2;   B(2,1+j) = dNdX/SQ2;  // SQ2 => Mandel representation
			}
			return;
		}
		case 1: // 1D
		case 4: // 2D(axis-symmetric)
		default:
			throw new Fatal("EquilibElem::B_Matrix: GeometryType==%d is not implemented yet",_geom);
	}
}


/* private */

inline void EquilibElem::_calc_initial_internal_forces()
{
	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> F(_n_dim*_n_nodes);
	F.SetValues(0.0);

	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs;  // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;       // Jacobian matrix
	LinAlg::Matrix<double> B;       // strain-displacement matrix
	LinAlg::Vector<double> sig;     // Stress vector in Mandel's notation 

	// Loop along integration points
	for (int i_ip=0; i_ip<_n_int_pts; ++i_ip)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i_ip].r;
		double s = _a_int_pts[i_ip].s;
		double t = _a_int_pts[i_ip].t; // only for 3D
		double w = _a_int_pts[i_ip].w;

		Derivs(r,s,0.0, derivs);      // Calculate Derivatives of Shape functions w.r.t local coordinate system
		if (_n_dim==3)
		Derivs(r,s,t, derivs);
		Jacobian(derivs, J);          // Calculate J (Jacobian) matrix for i_ip Integration Point
		B_Matrix(derivs, J, B);       // Calculate B matrix for i_ip Integration Point

		_a_model[i_ip]->Sig(sig); 

		// Calculate internal force vector;
		F += trn(B)*sig*det(J)*w;
	}

	// Update nodal NaturVals
	for (int i_node=0; i_node<_n_nodes; ++i_node)
	{
		// Assemble (local/element) displacements vector.
		_connects[i_node]->DOFVar("fx").NaturalVal += F(i_node*_n_dim  ); // NaturalVal must be set to zero during AddDOF routine
		_connects[i_node]->DOFVar("fy").NaturalVal += F(i_node*_n_dim+1);
		if (_n_dim==3)
		_connects[i_node]->DOFVar("fz").NaturalVal += F(i_node*_n_dim+2);
	}
}


/////////////////////////////////////////////////////////////////////////////////////////// Map /////


// Register the DOF information of EquilibElement into DOFInfoMap
int EquilibDOFInfoRegister()
{
	// Temporary 
	DOFInfo D; 

	// Nodal

	// Insert into DOFInfoMap
	DOFInfoMap["Equilibrium"] = D;

	return 0;
}

// Execute the autoregistration
int __EquilibElemDOFInfo_dummy_int  = EquilibDOFInfoRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_EQUILIB_H
