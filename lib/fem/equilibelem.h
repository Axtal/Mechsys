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
	Element * SetNode         (int iNodeLocal, FEM::Node * ptNode);
	void      UpdateState     (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
	void      BackupState     ();
	void      RestoreState    ();
	void      SetGeometryType (int Geom);  
	void      SetProperties   (Array<double> const & EleProps) { _unit_weight=EleProps[0]; }
	void      GetLabels       (Array<String> & Labels) const;
	void      Deactivate      ();

	// Derived methods to assemble DAS matrices
	size_t nOrder1Matrices () const { return 1; }
	void   Order1MatMap    (size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const;
	void   Order1Matrix    (size_t Index, LinAlg::Matrix<double> & Ke) const; // Stiffness

	// Methods
	void B_Matrix (LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const;

	// Access methods
	double Val (int iNodeLocal, char const * Name) const;
	double Val (                char const * Name) const;

private:
	// Data
	Array<EquilibModel*> _a_model;
	double               _unit_weight;

	// Private methods
	void _calc_initial_internal_forces ();

	// Private methods that MUST be derived
	virtual int _geom() const =0; ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)

}; // class EquilibElem


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */
	
// Derived methods

inline bool EquilibElem::IsReady() const
{
	if (_a_model.Size()==_n_int_pts && _connects.Size()==_n_nodes) return true;
	else return false;
}

inline bool EquilibElem::IsEssential(char const * DOFName) const
{
	if (strcmp(DOFName,"ux")==0 || strcmp(DOFName,"uy")==0) return true;
	if (_ndim==3                && strcmp(DOFName,"uz")==0) return true;
	return false;
}

inline void EquilibElem::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// Check _ndim
	if (_ndim<1) throw new Fatal("EquilibElem::SetModel: The space dimension (SetDim) must be set before calling this method");

	// If pointers to model was not already defined => No model was allocated
	if (_a_model.Size()==0)
	{
		// Resize the array of model pointers
		_a_model.Resize(_n_int_pts);

		// Loop along integration points
		for (size_t i=0; i<_n_int_pts; ++i)
		{
			// Allocate a new model and set parameters
			_a_model[i] = static_cast<EquilibModel*>(AllocModel(ModelName));
			_a_model[i]->SetGeom (_geom());
			_a_model[i]->SetPrms (Prms);
			_a_model[i]->SetInis (Inis);
		}

		// Calculate initial internal forces
		_calc_initial_internal_forces();
	}
	else throw new Fatal("EquilibElem::SetModel: Feature not implemented.");
}

inline Element * EquilibElem::SetNode(int iNodeLocal, FEM::Node * ptNode)
{
	// Connects
	_connects[iNodeLocal] = ptNode;

	// Add Degree of Freedom to a node (Essential, Natural)
	              _connects[iNodeLocal]->AddDOF("ux", "fx");
	              _connects[iNodeLocal]->AddDOF("uy", "fy");
	if (_ndim==3) _connects[iNodeLocal]->AddDOF("uz", "fz");

	// Shared
	_connects[iNodeLocal]->SetSharedBy(_my_id);

	return this;
}

inline void EquilibElem::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
{
	// Allocate (local/element) displacements vector
	LinAlg::Vector<double> dU(_ndim*_n_nodes); // Delta disp. of this element
	
	// Assemble (local/element) displacements vector
	for (size_t i=0; i<_n_nodes; ++i)
	{
		              dU(i*_ndim  ) = dUglobal(_connects[i]->DOFVar("ux").EqID);
		              dU(i*_ndim+1) = dUglobal(_connects[i]->DOFVar("uy").EqID);
		if (_ndim==3) dU(i*_ndim+2) = dUglobal(_connects[i]->DOFVar("uz").EqID);
	}
	
	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> dF(_ndim*_n_nodes); // Delta internal force of this element
	dF.SetValues(0.0);
	
	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs;  // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;       // Jacobian matrix
	LinAlg::Matrix<double> B;       // strain-displacement matrix
	LinAlg::Vector<double> DEps;    // Strain vector 
	LinAlg::Vector<double> DSig;    // Stress vector 

	// Loop along integration points
	for (size_t i=0; i<_n_int_pts; ++i)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i].r;
		double s = _a_int_pts[i].s;
		double t = _a_int_pts[i].t; // only for 3D cases
		double w = _a_int_pts[i].w;

		Derivs   (r,s,t, derivs);  // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);      // Calculate J (Jacobian) matrix for i Integration Point
		B_Matrix (derivs, J, B);   // Calculate B matrix for i Integration Point

		// Calculate a tensor for the increments of strain
		DEps = B*dU;
		
		// Update model
		_a_model[i]->StressUpdate(DEps, DSig);

		// Calculate internal force vector;
		dF += trn(B)*DSig*det(J)*w;
	}

	// Return internal forces
	for (size_t i=0; i<_n_nodes; ++i)
	{
		// Sum up contribution to internal forces vector
		              dFint(_connects[i]->DOFVar("fx").EqID) += dF(i*_ndim  );
		              dFint(_connects[i]->DOFVar("fy").EqID) += dF(i*_ndim+1);
		if (_ndim==3) dFint(_connects[i]->DOFVar("fz").EqID) += dF(i*_ndim+2);
	}
}

inline void EquilibElem::BackupState()
{
	for (size_t i=0; i<_n_int_pts; ++i)
		_a_model[i]->BackupState();
}

inline void EquilibElem::RestoreState()
{
	for (size_t i=0; i<_n_int_pts; ++i)
		_a_model[i]->RestoreState();
}

inline void EquilibElem::GetLabels(Array<String> & Labels) const
{
	// Get labels of all values to output
	switch (_geom()) // 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	{
		case 2: // 2D(plane-strain)
		{
			Labels.Resize(12); // 14 values to output
			Labels[0]="ux"; Labels[1]="uy";                                    // Displacements
			Labels[2]="fx"; Labels[3]="fy";                                    // Forces
			Labels[4]="Ex"; Labels[5]="Ey"; Labels[ 6]="Ez"; Labels[ 7]="Exy"; // Strains
			Labels[8]="Sx"; Labels[9]="Sy"; Labels[10]="Sz"; Labels[11]="Sxy"; // Stress
			return;
		}
		case 3: // 3D
		{
			Labels.Resize(18); // 18 values to output
			Labels[ 0]="ux"; Labels[ 1]="uy"; Labels[ 2]="uz";                                                       // Displacements
			Labels[ 3]="fx"; Labels[ 4]="fy"; Labels[ 5]="fz";                                                       // Forces
			Labels[ 6]="Ex"; Labels[ 7]="Ey"; Labels[ 8]="Ez"; Labels[ 9]="Exy"; Labels[10]="Eyz"; Labels[11]="Exz"; // Strains
			Labels[12]="Sx"; Labels[13]="Sy"; Labels[14]="Sz"; Labels[15]="Sxy"; Labels[16]="Syz"; Labels[17]="Sxz"; // Stress
			return;
		}
		case 5: // 2D(plane-stress)
		{
			Labels.Resize(10); // 10 values to output
			Labels[0]="ux"; Labels[1]="uy";                  // Displacements
			Labels[2]="fx"; Labels[3]="fy";                  // Forces
			Labels[4]="Ex"; Labels[5]="Ey"; Labels[6]="Exy"; // Strains
			Labels[7]="Sx"; Labels[8]="Sy"; Labels[9]="Sxy"; // Stress
			return;
		}
		case 1: // 1D
		case 4: // 2D(axis-symmetric)
		default:
			throw new Fatal("EquilibElem::GetLabels: GeometryType==%d is not implemented yet",_geom());
	}
}

inline double EquilibElem::Val(int iNodeLocal, char const * Name) const
{
	// Displacements
	if (strcmp(Name,"ux")==0 || strcmp(Name,"uy")==0 || strcmp(Name,"uz")==0)
		return _connects[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Forces
	else if (strcmp(Name,"fx")==0 || strcmp(Name,"fy")==0 || strcmp(Name,"fz")==0)
		return _connects[iNodeLocal]->DOFVar(Name).NaturalVal;

	// Stress, strains, internal values, etc.
	else
	{
		// Vectors for extrapolation
		LinAlg::Vector<double>    ip_values (_n_int_pts);
		LinAlg::Vector<double> nodal_values (_n_nodes);

		// Get integration point values
		for (size_t i=0; i<_n_int_pts; i++)
			ip_values(i) = _a_model[i]->Val(Name);

		// Extrapolation
		Extrapolate (ip_values, nodal_values);

		// Output single value
		return nodal_values (iNodeLocal);
	}
}

inline double EquilibElem::Val(char const * Name) const
{
	// Get integration point values
	double sum = 0.0;
	for (size_t i=0; i<_n_int_pts; i++)
		sum += _a_model[i]->Val(Name);

	// Output single value at CG
	return sum/_n_int_pts;
}

inline void EquilibElem::Deactivate()
{
	throw new Fatal("EquilibElem::Deactivate: Feature not implemented yet");
}

// Derived methods to assemble DAS matrices

inline void EquilibElem::Order1MatMap(size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const
{
	// Size of Ke
	int n_rows = _ndim*_n_nodes; // == n_cols

	// Mounting a map of positions from Ke to Global
	int idx_Ke = 0;                // position (idx) inside Ke matrix
	RowsMap       .Resize(n_rows); // size=Ke.Rows()=Ke.Cols()
	RowsEssenPresc.Resize(n_rows); // size=Ke.Rows()=Ke.Cols()

	// Fill map of Ke position to K position of DOFs components
	for (size_t i=0; i<_n_nodes; ++i)
	{
		RowsMap        [idx_Ke] = _connects[i]->DOFVar("ux").EqID; 
		RowsEssenPresc [idx_Ke] = _connects[i]->DOFVar("ux").IsEssenPresc; 
		idx_Ke++;
		RowsMap        [idx_Ke] = _connects[i]->DOFVar("uy").EqID; 
		RowsEssenPresc [idx_Ke] = _connects[i]->DOFVar("uy").IsEssenPresc; 
		idx_Ke++;
		if (_ndim==3)
		{
			RowsMap        [idx_Ke] = _connects[i]->DOFVar("uz").EqID; 
			RowsEssenPresc [idx_Ke] = _connects[i]->DOFVar("uz").IsEssenPresc; 
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
	Ke.Resize(_ndim*_n_nodes, _ndim*_n_nodes); // sum(Bt*D*B*det(J)*w)
	Ke.SetValues(0.0);

	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs; // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;      // Jacobian matrix
	LinAlg::Matrix<double> B;      // strain-displacement matrix
	LinAlg::Matrix<double> D;      // Constitutive matrix

	// Loop along integration points
	for (size_t i=0; i<_n_int_pts; ++i)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i].r;
		double s = _a_int_pts[i].s;
		double t = _a_int_pts[i].t;
		double w = _a_int_pts[i].w;

		Derivs   (r,s,t, derivs); // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);     // Calculate J (Jacobian) matrix for i Integration Point
		B_Matrix (derivs,J, B);   // Calculate B matrix for i Integration Point

		// Constitutive tensor 
		_a_model[i]->TgStiffness(D); 

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

	// geometry type: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	switch (_geom())
	{
		case 1: // 1D
		{
			// Derivatives and determinand of Jacobian
			LinAlg::Matrix<double> nat_derivs(_ndim, _ndim*_n_nodes); 
			nat_derivs.SetValues(0.0);
			double det_J = det(J);
			for (size_t i=0; i<_n_nodes; i++)
			for (int    j=0; j<_ndim;    j++)
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
			B.Resize(n_scomps,_ndim*_n_nodes);
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
			B.Resize(n_scomps,_ndim*_n_nodes);
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
			B.Resize(n_scomps,_ndim*_n_nodes);
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

inline void EquilibElem::_calc_initial_internal_forces()
{
	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> F(_ndim*_n_nodes);
	F.SetValues(0.0);

	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs;  // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;       // Jacobian matrix
	LinAlg::Matrix<double> B;       // strain-displacement matrix
	LinAlg::Vector<double> sig;     // Stress vector in Mandel's notation 

	// Loop along integration points
	for (size_t i=0; i<_n_int_pts; ++i)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i].r;
		double s = _a_int_pts[i].s;
		double t = _a_int_pts[i].t; // only for 3D
		double w = _a_int_pts[i].w;

		Derivs   (r,s,t, derivs); // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);     // Calculate J (Jacobian) matrix for i Integration Point
		B_Matrix (derivs, J, B);  // Calculate B matrix for i Integration Point

		_a_model[i]->Sig(sig); 

		// Calculate internal force vector;
		F += trn(B)*sig*det(J)*w;
	}

	// Update nodal NaturVals
	for (size_t i=0; i<_n_nodes; ++i)
	{
		// Assemble (local/element) displacements vector.
		              _connects[i]->DOFVar("fx").NaturalVal += F(i*_ndim  ); // NaturalVal must be set to zero during AddDOF routine
		              _connects[i]->DOFVar("fy").NaturalVal += F(i*_ndim+1);
		if (_ndim==3) _connects[i]->DOFVar("fz").NaturalVal += F(i*_ndim+2);
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
