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


/* __ Equilibrium capable element __ */


#ifndef MECHSYS_FEM_EQUI2D_H
#define MECHSYS_FEM_EQUI2D_H

// MechSys
#include "fem/element.h"
#include "models/equilibs/linelastic.h"
#include "tensors/tensors.h"
#include "tensors/functions.h"
#include "util/numstreams.h"
#include "util/string.h"
#include "linalg/laexpr.h"

namespace FEM
{

class EquilibElem2D : public virtual Element
{
public:
	// EquilibElem2D constants
	static int    NDIM;
	static int    NSTRESSCOMPS;
	static String DUX;
	static String DUY;
	static String DFX;
	static String DFY;
	static String DTX;
	static String DTY;
	
	// Destructor
	virtual ~EquilibElem2D() {}

	// Derived methods
	bool      IsEssential     (String const & DOFName) const;
	void      ReAllocateModel (String const & ModelName, String const & Prms, String const & Inis);
	Element * SetNode         (int iNodeLocal, int iNodeGlobal);
	void      UpdateState     (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
	void      BackupState     ();
	void      RestoreState    ();
	void      SetProperties   (Array<double> const & EleProps) { _unit_weight=EleProps[0]; }
	String    OutCenter       (bool PrintCaptionOnly) const;
	void      OutNodes        (LinAlg::Matrix<double> & Values, Array<String> & Labels) const;
	void      Deactivate      ();
	void      FaceNodalVals   (String const & FaceDOFName, double const FaceDOFValue, Array<FEM::Node*> const & APtrFaceNodes, String & NodalDOFName, LinAlg::Vector<double> & NodalValues) const;

	// Derived methods to assemble DAS matrices
	size_t nOrder1Matrices () const { return 1; }
	void   Order1MatMap    (size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const;
	void   Order1Matrix    (size_t Index, LinAlg::Matrix<double> & Ke) const; // Stiffness

	// Derived methods to output
	double OutScalar2 ()             const; ///< Ed
	void   OutTensor1 (String & Str) const; ///< Stress
	void   OutTensor2 (String & Str) const; ///< Strain

	// Methods
	void B_Matrix (LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const;

private:
	// Data
	Array<EquilibModel*> _a_model;
	double               _unit_weight;

	// Private methods
	void _calc_initial_internal_forces ();

}; // class EquilibElem2D

// EquilibElem2D constants
int    EquilibElem2D::NDIM         = 2;
int    EquilibElem2D::NSTRESSCOMPS = 3;
String EquilibElem2D::DUX          = _T("Dux");
String EquilibElem2D::DUY          = _T("Duy");
String EquilibElem2D::DFX          = _T("Dfx");
String EquilibElem2D::DFY          = _T("Dfy");
String EquilibElem2D::DTX          = _T("Dtx");
String EquilibElem2D::DTY          = _T("Dty");


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */
	
// Derived methods

inline bool EquilibElem2D::IsEssential(String const & DOFName) const
{
	if (DOFName==DUX || DOFName==DUY) return true;
	else return false;
}

inline void EquilibElem2D::ReAllocateModel(String const & ModelName, String const & Prms, String const & Inis)
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
			_a_model[i] = AllocEquilibModel(ModelName);
			_a_model[i]->SetGeom(2);
			_a_model[i]->SetPrms(Prms);
			_a_model[i]->SetInis(Inis);
		}

		// Calculate initial internal forces
		_calc_initial_internal_forces();
	}
	else throw new Fatal("EquilibElem2D::ReAllocateModel: Feature not implemented.");
}

inline Element * EquilibElem2D::SetNode(int iNodeLocal, int iNodeGlobal)
{
	// Connects
	_connects[iNodeLocal] = Nodes[iNodeGlobal];

	// Add Degree of Freedom to a node (Essential, Natural)
	Nodes[iNodeGlobal]->AddDOF(DUX, DFX);
	Nodes[iNodeGlobal]->AddDOF(DUY, DFY);

	// Shared
	Nodes[iNodeGlobal]->SetSharedBy(_my_id);

	return this;
}

inline void EquilibElem2D::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
{
	// Allocate (local/element) displacements vector
	LinAlg::Vector<double> dU(NDIM*_n_nodes); // Delta disp. of this element
	
	// Assemble (local/element) displacements vector
	for (int i=0; i<_n_nodes; ++i)
	{
		dU(i*NDIM  ) = dUglobal(_connects[i]->DOFVar(DUX).EqID);
		dU(i*NDIM+1) = dUglobal(_connects[i]->DOFVar(DUY).EqID);
	}
	
	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> dF(NDIM*_n_nodes); // Delta internal force of this element
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
		double w = _a_int_pts[i_ip].w;

		Derivs   (r,s,0.0, derivs);  // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);        // Calculate J (Jacobian) matrix for i_ip Integration Point
		B_Matrix (derivs, J, B);     // Calculate B matrix for i_ip Integration Point

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
		dFint(_connects[i]->DOFVar(DFX).EqID) += dF(i*NDIM  );
		dFint(_connects[i]->DOFVar(DFY).EqID) += dF(i*NDIM+1);
	}
}

inline void EquilibElem2D::BackupState()
{
	for (int i=0; i<_n_int_pts; ++i)
		_a_model[i]->BackupState();
}

inline void EquilibElem2D::RestoreState()
{
	for (int i=0; i<_n_int_pts; ++i)
		_a_model[i]->RestoreState();
}

inline String EquilibElem2D::OutCenter(bool PrintCaptionOnly=false) const
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
		oss << Util::_8s<< "Sx" << Util::_8s<< "Sy" << Util::_8s<< "Sxy";
		oss << Util::_8s<< "Ex" << Util::_8s<< "Ey" << Util::_8s<< "Exy";

		// Internal state values
		Array<String> str_state_names;   _a_model[0]->InternalStateNames(str_state_names);
		for (int i=0; i<n_int_state_vals; ++i)
			oss << Util::_8s<< str_state_names[i];
		oss << std::endl;
	}
	else
	{
		// Stress, strains and internal state values evaluated at the center of the element
		Vector<double> sig_cen;
		Vector<double> eps_cen;
		sig_cen.SetValues(0.0);
		eps_cen.SetValues(0.0);
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
		double sq2 = sqrt(2.0);
		oss << Util::_8s<< p          << Util::_8s<< q          << Util::_8s<< sin3th << Util::_8s<< Ev*100.0       << Util::_8s<< Ed*100.0;
		oss << Util::_8s<< sig_cen(0) << Util::_8s<< sig_cen(1) << Util::_8s<< sig_cen(2)         << Util::_8s<< sig_cen(3)/sq2 << Util::_8s<< sig_cen(4)/sq2 << Util::_8s<< sig_cen(5)/sq2;
		oss << Util::_8s<< eps_cen(0) << Util::_8s<< eps_cen(1) << Util::_8s<< eps_cen(2)         << Util::_8s<< eps_cen(3)/sq2 << Util::_8s<< eps_cen(4)/sq2 << Util::_8s<< eps_cen(5)/sq2;
		for (int j=0; j<n_int_state_vals; ++j)
			oss << Util::_8s<< int_state_vals_cen[j];
		oss << std::endl;
	} */

	return oss.str(); 
}

inline void EquilibElem2D::OutNodes(LinAlg::Matrix<double> & Values, Array<String> & Labels) const
{
	int const DATA_COMPS=10;
	Values.Resize(_n_nodes,DATA_COMPS);
	Labels.Resize(DATA_COMPS);
	Labels[ 0] = DUX ; Labels[ 1] = DUY ;  
	Labels[ 2] = DFX ; Labels[ 3] = DFY ; 
	Labels[ 4] = "Ex"; Labels[ 5] = "Ey"; Labels[ 6] = "Exy";
	Labels[ 7] = "Sx"; Labels[ 8] = "Sy"; Labels[ 9] = "Sxy";
	for (int i_node=0; i_node<_n_nodes; i_node++)
	{
		Values(i_node,0) = _connects[i_node]->DOFVar(DUX).EssentialVal;
		Values(i_node,1) = _connects[i_node]->DOFVar(DUY).EssentialVal;
		Values(i_node,2) = _connects[i_node]->DOFVar(DFX).NaturalVal;
		Values(i_node,3) = _connects[i_node]->DOFVar(DFY).NaturalVal;
	}

	//Extrapolation
	LinAlg::Vector<double> ip_values(_n_int_pts);
	LinAlg::Vector<double> nodal_values(_n_nodes);
	
	// Strains
	for (int i_comp=0; i_comp<NSTRESSCOMPS; i_comp++)
	{
		for (int j_ip=0; j_ip<_n_int_pts; j_ip++)
		{
			Vector<double> eps;
			_a_model[j_ip]->Eps(eps);
			ip_values(j_ip) = eps(i_comp); //getting IP values
		}
		Extrapolate(ip_values, nodal_values);
		for (int j_node=0; j_node<_n_nodes; j_node++)
			Values(j_node,i_comp+4) = nodal_values(j_node);
	}
	
	// Stresses
	for (int i_comp=0; i_comp<NSTRESSCOMPS; i_comp++)
	{
		for (int j_ip=0; j_ip<_n_int_pts; j_ip++)
		{
			Vector<double> sig;
			_a_model[j_ip]->Sig(sig);
			ip_values(j_ip) = sig(i_comp); //getting IP values
		}
		Extrapolate(ip_values, nodal_values);
		for (int j_node=0; j_node<_n_nodes; j_node++)
			Values(j_node,i_comp+7) = nodal_values(j_node);
	}

}

inline void EquilibElem2D::Deactivate()
{
	throw new Fatal("EquilibElem2D::Deactivate: Feature not implemented yet");
}

inline void EquilibElem2D::FaceNodalVals(String const & FaceDOFName, double const   FaceDOFValue, Array<FEM::Node*> const & APtrFaceNodes, String & NodalDOFName, LinAlg::Vector<double>& NodalValues) const
{
	if (FaceDOFName==DTX || FaceDOFName==DTY)
	{
		if (FaceDOFName==DTX) NodalDOFName=DFX;
		if (FaceDOFName==DTY) NodalDOFName=DFY;
		Dist2FaceNodes(APtrFaceNodes, FaceDOFValue, NodalValues);
	}
	else
	{
		std::ostringstream oss; oss << "Face nodes coordinates:\n";
		for (size_t i_node=0; i_node<APtrFaceNodes.Size(); ++i_node)
			oss << "X=" << APtrFaceNodes[i_node]->X() << ", Y=" << APtrFaceNodes[i_node]->Y() << std::endl;
		throw new Fatal(_("EquilibElem2D::CalcFaceNodalValues: This method must only be called for FaceDOFName< %s > equal to Dtx or Dty.\n %s"),
				FaceDOFName.c_str(), oss.str().c_str());
	}
}

// Derived methods to assemble DAS matrices

inline void EquilibElem2D::Order1MatMap(size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const
{
	// Size of Ke
	int n_rows = NDIM*_n_nodes; // == n_cols

	// Mounting a map of positions from Ke to Global
	int idx_Ke = 0;                // position (idx) inside Ke matrix
	RowsMap       .Resize(n_rows); // size=Ke.Rows()=Ke.Cols()
	RowsEssenPresc.Resize(n_rows); // size=Ke.Rows()=Ke.Cols()

	// Fill map of Ke position to K position of DOFs components
	for (int i_node=0; i_node<_n_nodes; ++i_node)
	{
		RowsMap        [idx_Ke] = _connects[i_node]->DOFVar(DUX).EqID; 
		RowsEssenPresc [idx_Ke] = _connects[i_node]->DOFVar(DUX).IsEssenPresc; 
		idx_Ke++;
		RowsMap        [idx_Ke] = _connects[i_node]->DOFVar(DUY).EqID; 
		RowsEssenPresc [idx_Ke] = _connects[i_node]->DOFVar(DUY).IsEssenPresc; 
		idx_Ke++;
	}
	ColsMap        = RowsMap;
	ColsEssenPresc = RowsEssenPresc;
}

inline void EquilibElem2D::Order1Matrix(size_t index, LinAlg::Matrix<double> & Ke) const
{
	/* Stiffness:
	   ==========
	
	                 /    T
	        [Ke]  =  | [B]  * [D] * [B]  * dV
	                 /
	*/

	// Resize Ke
	Ke.Resize(NDIM*_n_nodes, NDIM*_n_nodes); // sum(Bt*D*B*det(J)*w)
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
		double w = _a_int_pts[i_ip].w;

		Derivs(r,s,0.0, derivs);        // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian(derivs, J);          // Calculate J (Jacobian) matrix for i_ip Integration Point
		B_Matrix(derivs,J, B);        // Calculate B matrix for i_ip Integration Point

		// Constitutive tensor 
		_a_model[i_ip]->TgStiffness(D); 

		// Calculate Tangent Stiffness
		Ke += trn(B)*D*B*det(J)*w;
	}
}
	
// Derived methods to output

inline void EquilibElem2D::OutTensor1(String & Str) const
{ /*
	// Stress evaluated at the center of the element
	Tensors::Tensor2 s(0.0);

	// Loop over integration points
	for (int i_ip=0; i_ip<_n_int_pts; ++i_ip)
		s += _a_model[i_ip]->Sig();
	
	// Average stress
	s = s / _n_int_pts;

	// Output
	double sq2 = sqrt(2.0);
	Str.Printf(_(" %e %e %e  %e %e %e  %e %e %e "), s(0),s(3)/sq2,s(5)/sq2,  s(3)/sq2,s(1),s(4)/sq2,  s(5)/sq2,s(4)/sq2,s(2));
*/
}

inline void EquilibElem2D::OutTensor2(String & Str) const
{ /*
	// Strains evaluated at the center of the element
	Tensors::Tensor2 e(0.0);

	// Loop over integration points
	for (int i_ip=0; i_ip<_n_int_pts; ++i_ip)
		e += 100.0*_a_model[i_ip]->Eps();
	
	// Average strains
	e = e / _n_int_pts;

	// Output
	double sq2 = sqrt(2.0);
	Str.Printf(_(" %e %e %e  %e %e %e  %e %e %e "), e(0),e(3)/sq2,e(5)/sq2,  e(3)/sq2,e(1),e(4)/sq2,  e(5)/sq2,e(4)/sq2,e(2));
*/
}

inline double EquilibElem2D::OutScalar2() const
{ /*
	// Strains evaluated at the center of the element
	Tensors::Tensor2 e(0.0);

	// Loop over integration points
	for (int i_ip=0; i_ip<_n_int_pts; ++i_ip)
		e += 100.0*_a_model[i_ip]->Eps();
	
	// Average strains
	e = e / _n_int_pts;

	// Calculate strain invariants
	double   Ev,Ed;
	Strain_Ev_Ed(e,Ev,Ed);

	// Output
	return Ed;
*/
	return 0;
}

// Methods

inline void EquilibElem2D::B_Matrix(LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const
{
	/* OBS.:
	 *          This B matrix considers Soil Mechanics sign convention of stress and strains
	 *          Ex.: Compressive stresses/strains are positive
	 */
	
	// Resize B matrix
	B.Resize(NSTRESSCOMPS,NDIM*_n_nodes);

	// Cartesian derivatives
	LinAlg::Matrix<double> cart_derivs;
	cart_derivs = inv(J)*derivs;

	// Loop along all nodes of the element
	double dNdX,dNdY;
	int  j=0; // j column of B
	for (int i=0; i<_n_nodes; ++i) // i row of B
	{
		// Assemble B matrix
		j=i*NDIM;
		dNdX=-cart_derivs(0,i);  dNdY=-cart_derivs(1,i);  // Negative values => Soil mechanics convention
		B(0,0+j) =     dNdX;     B(0,1+j) =      0.0;     
		B(1,0+j) =      0.0;     B(1,1+j) =     dNdY;     
		B(2,0+j) =     dNdY;     B(2,1+j) =     dNdX;     
	}
}


/* private */

inline void EquilibElem2D::_calc_initial_internal_forces()
{
	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> F(NDIM*_n_nodes);
	F.SetValues(0.0);

	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs;  // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;       // Jacobian matrix
	LinAlg::Matrix<double> B;       // strain-displacement matrix

	// Loop along integration points
	for (int i_ip=0; i_ip<_n_int_pts; ++i_ip)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i_ip].r;
		double s = _a_int_pts[i_ip].s;
		double w = _a_int_pts[i_ip].w;

		Derivs(r,s,0.0, derivs);      // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian(derivs, J);          // Calculate J (Jacobian) matrix for i_ip Integration Point
		B_Matrix(derivs, J, B);       // Calculate B matrix for i_ip Integration Point

		Vector<double> sig;
		_a_model[i_ip]->Sig(sig); // Declare a reference to the actual (initial) stress tensor inside EquilibModel

		// Calculate internal force vector;
		F += trn(B)*sig*det(J)*w;
	}

	// Update nodal NaturVals
	for (int i_node=0; i_node<_n_nodes; ++i_node)
	{
		// Assemble (local/element) displacements vector. OBS.: NDIM=2
		_connects[i_node]->DOFVar(DFX).NaturalVal += F(i_node*NDIM  ); // NaturalVal must be set to zero during AddDOF routine
		_connects[i_node]->DOFVar(DFY).NaturalVal += F(i_node*NDIM+1);
	}
}


/////////////////////////////////////////////////////////////////////////////////////////// Map /////


// Register the DOF information of EquilibElement into DOFInfoMap
int Equilib2DDOFInfoRegister()
{
	// Temporary 
	DOFInfo D; 

	// Nodal
	D.NodeEssential.Push(EquilibElem2D::DUX + _("@Nodal displacement increment in x direction"));
	D.NodeEssential.Push(EquilibElem2D::DUY + _("@Nodal displacement increment in y direction"));
	D.NodeNatural  .Push(EquilibElem2D::DFX + _("@Nodal force increment in x direction"));
	D.NodeNatural  .Push(EquilibElem2D::DFY + _("@Nodal force increment in y direction"));

	// Face
	D.FaceEssential.Push(EquilibElem2D::DUX + _("@Displacement increment in x direction on face"));
	D.FaceEssential.Push(EquilibElem2D::DUY + _("@Displacement increment in y direction on face"));
	D.FaceNatural  .Push(EquilibElem2D::DTX + _("@Traction increment in x direction on face"));
	D.FaceNatural  .Push(EquilibElem2D::DTY + _("@Traction increment in y direction on face"));

	// Insert into DOFInfoMap
	DOFInfoMap["Equilibrium2D"] = D;

	return 0;
}

// Execute the autoregistration
int __Equilib2DElemDOFInfo_dummy_int  = Equilib2DDOFInfoRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_EQUI2D_H
