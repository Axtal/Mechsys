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
	// Constants
	static int NDIM;
	
	// Constructor
	ElasticRod();

	// Destructor
	virtual ~ElasticRod() {}

	// Derived methods
	String    Name            () const { return "ElasticRod"; }
	bool      IsEssential     (char const * DOFName) const;
	void      SetModel        (char const * ModelName, char const * Prms, char const * Inis);
	Element * SetNode         (int iNodeLocal, int iNodeGlobal);
	void      UpdateState     (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
	void      BackupState     () { _N_bkp = _N; }
	void      RestoreState    () { _N = _N_bkp; }
	void      SetGeometryType (int Geom) { _geom = Geom; }
	void      SetProperties   (Array<double> const & EleProps) { _unit_weight=EleProps[0]; }
	String    OutCenter       (bool PrintCaptionOnly) const { String res; res.Printf("%f",_N); return res; }
	void      OutNodes        (LinAlg::Matrix<double> & Values, Array<String> & Labels) const;

	// Derived methods (GEOMETRIC)
	int  VTKCellType    () const { return 4; } // VTK_
	void Shape          (double r, double s, double t, LinAlg::Vector<double> & Shape) const {}
	void Derivs         (double r, double s, double t, LinAlg::Matrix<double> & Derivs) const {}
	void FaceShape      (double r, double s, LinAlg::Vector<double> & Shape) const {}
	void FaceDerivs     (double r, double s, LinAlg::Matrix<double> & Derivs) const {}
	void Dist2FaceNodes (Array<Node*> const & FaceConnects, double const FaceValue, LinAlg::Vector<double> & NodalValues) const {}

	// Derived methods to assemble DAS matrices
	size_t nOrder1Matrices () const { return 1; };
	void   Order1MatMap    (size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const;
	void   Order1Matrix    (size_t Index, LinAlg::Matrix<double> & Ke) const; // Stiffness

private:
	// Data
	double _E;                 ///< Young modulus
	double _N;                 ///< Normal stress (compressive is positive)
	double _N_bkp;             ///< Backup of normal stress
	double _A;                 ///< Cross-sectional area
	double _unit_weight;       ///< Unit weight
	bool   _initial_state_set; ///< Was the initial state set?

	// Private methods
	void _calc_initial_internal_forces ();

}; // class ElasticRod

int ElasticRod::NDIM = 3;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */
	
inline ElasticRod::ElasticRod()
	: _initial_state_set(false)
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
}

// Derived methods

inline bool ElasticRod::IsEssential(char const * DOFName) const
{
	if (strncmp(DOFName,"ux",2)==0 || strncmp(DOFName,"uy",2)==0 || strncmp(DOFName,"uz",2)==0) return true;
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
	if (count!=1) throw Fatal("ElasticRod::SetModel: Parameters definition is incorrect. The syntax must be as in:\n\t E=10000.0\n");

	/* Inis   "N=0.0  A=1.0" */
	lp.Reset(Inis);
	lp.BreakExpressions(names,values);

	// Check initial values
	count = 0;
	if (names.Size()==2 && values.Size()==2)
	{
		for (size_t i=0; i<names.Size(); ++i)
		{
				 if (names[i]=="N") { _N = values[i];  count++; }
			else if (names[i]=="A") { _A = values[i];  count++; }
		}
	}
	if (count!=2) throw Fatal("LinElastic::SetInis: Initial values definition is incorrect. The syntax must be as in:\n\t N=0.0  A=1.0\n");

	// If the initial state was not set, compute initial internal forces
	if (!_initial_state_set) _calc_initial_internal_forces();
}

inline Element * ElasticRod::SetNode(int iNodeLocal, int iNodeGlobal)
{
	// Connects
	_connects[iNodeLocal] = Nodes[iNodeGlobal];

	// Add Degree of Freedom to a node (Essential, Natural)
	Nodes[iNodeGlobal]->AddDOF("ux", "fx");
	Nodes[iNodeGlobal]->AddDOF("uy", "fy");
	Nodes[iNodeGlobal]->AddDOF("uz", "fz");

	// Shared
	Nodes[iNodeGlobal]->SetSharedBy(_my_id);

	return this;
}

inline void ElasticRod::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
{
	// Allocate (local/element) displacements vector
	LinAlg::Vector<double> dU(NDIM*_n_nodes); // Delta disp. of this element
	
	// Assemble (local/element) displacements vector
	for (int i=0; i<_n_nodes; ++i)
	{
		dU(i*NDIM  ) = dUglobal(_connects[i]->DOFVar("ux").EqID);
		dU(i*NDIM+1) = dUglobal(_connects[i]->DOFVar("uy").EqID);
		dU(i*NDIM+2) = dUglobal(_connects[i]->DOFVar("uz").EqID);
	}
	
	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> dF(NDIM*_n_nodes); // Delta internal force of this element
	dF.SetValues(0.0);

	// Compute dF
	LinAlg::Matrix<double> Ke;
	Order1Matrix(0,Ke);
	LinAlg::Gemv(1.0,Ke,dU, 0.0,dF); // dF <- Ke*dU
	_N += dF(1)-dF(0);

	// Return internal forces
	for (int i=0; i<_n_nodes; ++i)
	{
		// Sum up contribution to internal forces vector
		dFint(_connects[i]->DOFVar("fx").EqID) += dF(i*NDIM  );
		dFint(_connects[i]->DOFVar("fy").EqID) += dF(i*NDIM+1);
		dFint(_connects[i]->DOFVar("fz").EqID) += dF(i*NDIM+2);
	}
}

inline void ElasticRod::OutNodes(LinAlg::Matrix<double> & Values, Array<String> & Labels) const
{
	int const DATA_COMPS=6;
	Values.Resize(_n_nodes,DATA_COMPS);
	Labels.Resize(DATA_COMPS);
	Labels[ 0] = "ux"; Labels[ 1] = "uy"; Labels[ 2] = "uz"; 
	Labels[ 3] = "fx"; Labels[ 4] = "fy"; Labels[ 5] = "fz";
	for (int i_node=0; i_node<_n_nodes; i_node++)
	{
		Values(i_node,0) = _connects[i_node]->DOFVar("ux").EssentialVal;
		Values(i_node,1) = _connects[i_node]->DOFVar("uy").EssentialVal;
		Values(i_node,2) = _connects[i_node]->DOFVar("uz").EssentialVal;
		Values(i_node,3) = _connects[i_node]->DOFVar("fx").NaturalVal;
		Values(i_node,4) = _connects[i_node]->DOFVar("fy").NaturalVal;
		Values(i_node,5) = _connects[i_node]->DOFVar("fz").NaturalVal;
	}
}

// Derived methods to assemble DAS matrices

inline void ElasticRod::Order1MatMap(size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const
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
		RowsMap        [idx_Ke] = _connects[i_node]->DOFVar("ux").EqID; 
		RowsEssenPresc [idx_Ke] = _connects[i_node]->DOFVar("ux").IsEssenPresc; 
		idx_Ke++;
		RowsMap        [idx_Ke] = _connects[i_node]->DOFVar("uy").EqID; 
		RowsEssenPresc [idx_Ke] = _connects[i_node]->DOFVar("uy").IsEssenPresc; 
		idx_Ke++;
		RowsMap        [idx_Ke] = _connects[i_node]->DOFVar("uz").EqID; 
		RowsEssenPresc [idx_Ke] = _connects[i_node]->DOFVar("uz").IsEssenPresc; 
		idx_Ke++;
	}
	ColsMap        = RowsMap;
	ColsEssenPresc = RowsEssenPresc;
}

inline void ElasticRod::Order1Matrix(size_t index, LinAlg::Matrix<double> & Ke) const
{
	// Resize Ke
	Ke.Resize(NDIM*_n_nodes, NDIM*_n_nodes); // sum(Bt*D*B*det(J)*w)
	Ke.SetValues(0.0);

	double x21 = _connects[1]->X()-_connects[0]->X();
	double y21 = _connects[1]->Y()-_connects[0]->Y();
	double z21 = _connects[1]->Z()-_connects[0]->Z();
	double L   = sqrt(x21*x21+y21*y21+z21*z21);
	       Ke  =  x21*x21, x21*y21, x21*z21,-x21*x21,-x21*y21,-x21*z21,
                  y21*x21, y21*y21, y21*z21,-y21*x21,-y21*y21,-y21*z21,
                  z21*x21, z21*y21, z21*z21,-z21*x21,-z21*y21,-z21*z21,
                 -x21*x21,-x21*y21,-x21*z21, x21*x21, x21*y21, x21*z21,
                 -y21*x21,-y21*y21,-y21*z21, y21*x21, y21*y21, y21*z21,
                 -z21*x21,-z21*y21,-z21*z21, z21*x21, z21*y21, z21*z21;
	       Ke  = Ke * (_E*_A/pow(L,3.0));
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
