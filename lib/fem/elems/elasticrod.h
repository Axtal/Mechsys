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

#ifndef MECHSYS_FEM_ELASTICROD_H
#define MECHSYS_FEM_ELASTICROD_H

// MechSys
#include "fem/element.h"
#include "util/lineparser.h"

using std::cout;
using std::endl;

namespace FEM
{

class ElasticRod : public Element
{

public:
	// Constructor
	ElasticRod();

	// Destructor
	virtual ~ElasticRod() {}

	// Derived methods
	String    Name            () const { return "ElasticRod"; }
	bool      IsReady         ()const;
	bool      IsEssential     (char const * DOFName) const;
	void      SetModel        (char const * ModelName, char const * Prms, char const * Inis);
	Element * SetNode         (int iNodeLocal, int iNodeGlobal);
	void      UpdateState     (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
	void      BackupState     () { _Sa_bkp=_Sa;  _Ea_bkp=_Ea; }
	void      RestoreState    () { _Sa=_Sa_bkp;  _Ea=_Ea_bkp; }
	void      SetGeometryType (int Geom) { _geom=Geom;  _ndim_prob=(_geom==3?3:2); }
	void      SetProperties   (Array<double> const & EleProps) { _unit_weight=EleProps[0]; }
	void      GetLabels       (Array<String> & Labels) const;

	// Derived methods (GEOMETRIC)
	int  VTKCellType () const { return 4; } // VTK_
	void Shape       (double r, double s, double t, LinAlg::Vector<double> & Shape) const {}
	void Derivs      (double r, double s, double t, LinAlg::Matrix<double> & Derivs) const {}
	void FaceShape   (double r, double s, LinAlg::Vector<double> & Shape) const {}
	void FaceDerivs  (double r, double s, LinAlg::Matrix<double> & Derivs) const {}

	// Derived methods to assemble DAS matrices
	size_t nOrder1Matrices () const { return 1; };
	void   Order1MatMap    (size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const;
	void   Order1Matrix    (size_t Index, LinAlg::Matrix<double> & Ke) const; // Stiffness

	// Access methods
	double Val(int iNodeLocal, char const * Name) const;
	double Val(                char const * Name) const;

private:
	// Data
	double _E;           ///< Young modulus
	double _Sa;          ///< Axial stress (compressive is positive)
	double _Ea;          ///< Axial strain
	double _Sa_bkp;      ///< Backup of axial stress
	double _Ea_bkp;      ///< Backup axial strain
	double _A;           ///< Cross-sectional area
	double _unit_weight; ///< Unit weight
	int    _ndim_prob;   ///< Problem dimension

	// Private methods
	void _calc_initial_internal_forces ();

}; // class ElasticRod


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */
	
inline ElasticRod::ElasticRod()
{
	// Setup nodes number
	_n_nodes        = 2;
	_n_int_pts      = 0;
	_n_face_nodes   = 0;
	_n_face_int_pts = 0;

	// Allocate nodes (connectivity)
	_connects.Resize(_n_nodes);

	// Setup pointer to the array of Integration Points
	_a_int_pts = NULL;

	// Default values
	_E           = 10000.0;
	_Sa          = 0.0;
	_Ea          = 0.0;
	_A           = 1.0;
	_unit_weight = 0.0;
}

// Derived methods

inline bool ElasticRod::IsReady() const
{
	if (_connects.Size()==static_cast<size_t>(_n_nodes) && _geom>0) return true;
	else return false;
}

inline bool ElasticRod::IsEssential(char const * DOFName) const
{
	if (strcmp(DOFName,"ux")==0 || strcmp(DOFName,"uy")==0) return true;
	if (_ndim_prob==3           && strcmp(DOFName,"uz")==0) return true;
	else return false;
}

inline void ElasticRod::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	/* Prms   "E=20000.0" */
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Check parameters
	int count = 0;
	if (names.Size()==1 && values.Size()==1)
		if (names[0]=="E") { _E = values[0];  count++; }
	if (count!=1) throw new Fatal("ElasticRod::SetModel: Parameters definition is incorrect. The syntax must be as in:\n\t E=10000.0\n");

	/* Inis   "Sa=0.0  A=1.0" */
	lp.Reset(Inis);
	lp.BreakExpressions(names,values);

	// Check initial values
	count = 0;
	if (names.Size()==2 && values.Size()==2)
	{
		for (size_t i=0; i<names.Size(); ++i)
		{
				 if (names[i]=="Sa") { _Sa = values[i];  count++; } // axial stress
			else if (names[i]=="A")  { _A  = values[i];  count++; } // area
		}
	}
	if (count!=2) throw new Fatal("LinElastic::SetInis: Initial values definition is incorrect. The syntax must be as in:\n\t Sa=0.0  A=1.0\n");

	// If the initial state was not set, compute initial internal forces
	_calc_initial_internal_forces();
}

inline Element * ElasticRod::SetNode(int iNodeLocal, int iNodeGlobal)
{
	// Connects
	_connects[iNodeLocal] = Nodes[iNodeGlobal];

	// Add Degree of Freedom to a node (Essential, Natural)
	                   Nodes[iNodeGlobal]->AddDOF("ux", "fx");
	                   Nodes[iNodeGlobal]->AddDOF("uy", "fy");
	if (_ndim_prob==3) Nodes[iNodeGlobal]->AddDOF("uz", "fz");

	// Shared
	Nodes[iNodeGlobal]->SetSharedBy(_my_id);

	return this;
}

inline void ElasticRod::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
{
	// Allocate (local/element) displacements vector
	LinAlg::Vector<double> dU(_ndim_prob*_n_nodes); // Delta disp. of this element
	
	// Assemble (local/element) displacements vector
	for (int i=0; i<_n_nodes; ++i)
	{
		                   dU(i*_ndim_prob  ) = dUglobal(_connects[i]->DOFVar("ux").EqID);
		                   dU(i*_ndim_prob+1) = dUglobal(_connects[i]->DOFVar("uy").EqID);
		if (_ndim_prob==3) dU(i*_ndim_prob+2) = dUglobal(_connects[i]->DOFVar("uz").EqID);
	}
	
	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> dF(_ndim_prob*_n_nodes); // Delta internal force of this element
	dF.SetValues(0.0);

	/*
	// Compute dF
	LinAlg::Matrix<double> Ke;
	Order1Matrix(0,Ke);
	LinAlg::Gemv(1.0,Ke,dU, 0.0,dF); // dF <- Ke*dU
	_Sa += dF(1)-dF(0);
	*/


	/*
	// Conversion matrix
	double L;
	LinAlg::Matrix<double> Te(_ndim_prob*_n_nodes, _ndim_prob*_n_nodes);
	if (_ndim_prob==2)
	{
		double x21 = _connects[1]->X()-_connects[0]->X();
		double y21 = _connects[1]->Y()-_connects[0]->Y();
		       L   = sqrt(x21*x21+y21*y21);
		       Te  =  x21/L, y21/L,   0.0,   0.0,
		             -y21/L, x21/L,   0.0,   0.0,
		                0.0,   0.0, x21/L, y21/L,
		                0.0,   0.0,-y21/L, x21/L;
	}
	else
	{
		double x21 = _connects[1]->X()-_connects[0]->X();
		double y21 = _connects[1]->Y()-_connects[0]->Y();
		double z21 = _connects[1]->Z()-_connects[0]->Z();
		       L   = sqrt(x21*x21+y21*y21+z21*z21);
	}

	// Local coordinates
	LinAlg::Vector<double>    dU_loc(_ndim_prob*_n_nodes);
	LinAlg::Vector<double> dFint_loc(_ndim_prob*_n_nodes); dFint_loc.SetValues(0.0);
	dU_loc = Te*dU;
	double elongation = dU_loc(0);
	_Sa += _E*elongation;
	dFint_loc(0) = _A*_E*elongation;
	dF = inv(Te)*dFint_loc;
	*/


	// Return internal forces
	for (int i=0; i<_n_nodes; ++i)
	{
		// Sum up contribution to internal forces vector
		                   dFint(_connects[i]->DOFVar("fx").EqID) += dF(i*_ndim_prob  );
		                   dFint(_connects[i]->DOFVar("fy").EqID) += dF(i*_ndim_prob+1);
		if (_ndim_prob==3) dFint(_connects[i]->DOFVar("fz").EqID) += dF(i*_ndim_prob+2);
	}
}

inline void ElasticRod::GetLabels(Array<String> & Labels) const
{
	// Get labels of all values to output
	switch (_geom) // 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	{
		case 2: // 2D(plane-strain)
		case 5: // 2D(plane-stress)
		{
			Labels.Resize(6); // 6 values to output
			Labels[0]="ux"; Labels[1]="uy"; // Displacements
			Labels[2]="fx"; Labels[3]="fy"; // Forces
			Labels[4]="Ea";                 // Strains (axial strain)
			Labels[5]="Sa";                 // Stress (axial stress)
			return;
		}
		case 3: // 3D
		{
			Labels.Resize(8); // 8 values to output
			Labels[0]="ux"; Labels[1]="uy"; Labels[2]="uz"; // Displacements
			Labels[3]="fx"; Labels[4]="fy"; Labels[5]="fz"; // Forces
			Labels[6]="Ea";                                 // Strains (axial strain)
			Labels[7]="Sa";                                 // Stress (axial stress)
			return;
		}
		case 1: // 1D
		case 4: // 2D(axis-symmetric)
		default:
			throw new Fatal("EquilibElem::GetLabels: GeometryType==%d is not implemented yet",_geom);
	}
}

inline double ElasticRod::Val(int iNodeLocal, char const * Name) const
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
		     if (strcmp(Name,"Sa")==0) return _Sa;
		else if (strcmp(Name,"Ea")==0) return _Ea;
		else throw new Fatal("ElasticRod::Val: value named %s is not available for this element",Name);
	}
}

inline double ElasticRod::Val(char const * Name) const
{
		 if (strcmp(Name,"Sa")==0) return _Sa;
	else if (strcmp(Name,"Ea")==0) return _Ea;
	else throw new Fatal("ElasticRod::Val: value named %s is not available for this element",Name);
}

// Derived methods to assemble DAS matrices

inline void ElasticRod::Order1MatMap(size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const
{
	// Size of Ke
	int n_rows = _ndim_prob*_n_nodes; // == n_cols

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
		if (_ndim_prob==3)
		{
			RowsMap        [idx_Ke] = _connects[i_node]->DOFVar("uz").EqID; 
			RowsEssenPresc [idx_Ke] = _connects[i_node]->DOFVar("uz").IsEssenPresc; 
			idx_Ke++;
		}
	}
	ColsMap        = RowsMap;
	ColsEssenPresc = RowsEssenPresc;
}

inline void ElasticRod::Order1Matrix(size_t index, LinAlg::Matrix<double> & Ke) const
{
	// Resize Ke
	Ke.Resize(_ndim_prob*_n_nodes, _ndim_prob*_n_nodes); // sum(Bt*D*B*det(J)*w)
	Ke.SetValues(0.0);

	double L;
	if (_ndim_prob==2) // 2D
	{
		double x21 = _connects[1]->X()-_connects[0]->X();
		double y21 = _connects[1]->Y()-_connects[0]->Y();
		       L   = sqrt(x21*x21+y21*y21);
		       Ke  =  x21*x21, x21*y21,-x21*x21,-x21*y21,
		              y21*x21, y21*y21,-y21*x21,-y21*y21,
		             -x21*x21,-x21*y21, x21*x21, x21*y21,
		             -y21*x21,-y21*y21, y21*x21, y21*y21;
	}
	else // 3D
	{
		double x21 = _connects[1]->X()-_connects[0]->X();
		double y21 = _connects[1]->Y()-_connects[0]->Y();
		double z21 = _connects[1]->Z()-_connects[0]->Z();
		       L   = sqrt(x21*x21+y21*y21+z21*z21);
		       Ke  =  x21*x21, x21*y21, x21*z21,-x21*x21,-x21*y21,-x21*z21,
		              y21*x21, y21*y21, y21*z21,-y21*x21,-y21*y21,-y21*z21,
		              z21*x21, z21*y21, z21*z21,-z21*x21,-z21*y21,-z21*z21,
		             -x21*x21,-x21*y21,-x21*z21, x21*x21, x21*y21, x21*z21,
		             -y21*x21,-y21*y21,-y21*z21, y21*x21, y21*y21, y21*z21,
		             -z21*x21,-z21*y21,-z21*z21, z21*x21, z21*y21, z21*z21;
	}
	Ke = Ke * (_E*_A/pow(L,3.0));
}

/* private */

inline void ElasticRod::_calc_initial_internal_forces()
{
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new ElasticRod element
Element * ElasticRodMaker()
{
	return new ElasticRod();
}

// Register a ElasticRod element into ElementFactory array map
int ElasticRodRegister()
{
	ElementFactory["ElasticRod"] = ElasticRodMaker;
	return 0;
}

// Execute the autoregistration
int __ElasticRod_dummy_int  = ElasticRodRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_ELASTICROD_H
