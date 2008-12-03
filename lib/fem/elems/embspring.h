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

#ifndef MECHSYS_FEM_EMBSPRING_H
#define MECHSYS_FEM_EMBSPRING_H

// MechSys
#include "fem/equilibelem.h"
#include "util/exception.h"
#include "fem/elems/vtkCellType.h"

namespace FEM
{

class EmbSpring : public EquilibElem
{
public:
	// Constructor
	EmbSpring(Element * Tress, Node * Ext, LinAlg::Vector<double> & Direction);

	// Derived methods
	char const * Name() const { return "EmbSpring"; }

	// Derived methods
	bool   CheckModel   () const;
	void   SetModel     (char const * ModelName, char const * Prms, char const * Inis);
	void   UpdateState  (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
	void   CalcDepVars  () const;
	double Val          (int iNodeLocal, char const * Name) const;
	double Val          (char const * Name) const;
	void   Order1Matrix (size_t Index, LinAlg::Matrix<double> & Ke) const; ///< Stiffness
	void   B_Matrix     (LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const;
	int    VTKCellType  () const { return VTK_POLY_VERTEX; }
	void   VTKConnect   (String & Nodes) const;
	void   OutInfo(std::ostream & os) const;
	

private:
	// Data
	double _K;   ///< Spring stiffness
	double _LA;  ///< Representative contact area

	// Private methods
	int  _geom     () const { return 2; }              ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void _initialize();

	void _mount_T_matrix(Vector<double> const & Shape, Vector<double> const & Direction, Matrix<double> & T) const;
	void _calc_initial_internal_state();
	void _get_another_vector(Vector<double> const & V1, Vector<double> & V2) const;
	void _get_normal_vector (Vector<double> const & V1, Vector<double> const & V2, Vector<double> & V3) const;

	Element *              _tress;
	Node    *              _ext; 
	LinAlg::Vector<double> _direct;
	

}; // class EmbSpring


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */
EmbSpring::EmbSpring(Element * Tress, Node * Ext, LinAlg::Vector<double> & Direction): _K(-1), _LA(-1)
{
	_tress  = Tress;
	_ext    = Ext;
	_direct = Direction;
	_n_nodes = _tress->NNodes() + 1;
	_connects.Resize(_n_nodes);
}

inline bool EmbSpring::CheckModel() const
{
	if (_K<0.0 || _LA<0.0) return false;
	return true;
}

inline void EmbSpring::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// Check _ndim
	if (_ndim<1) throw new Fatal("EmbSpring::SetModel: The space dimension (SetDim) must be set before calling this method");
	if (CheckConnect()==false) throw new Fatal("EmbSpring::SetModel: Connectivity is not correct. Connectivity MUST be set before calling this method");

	/* "E=1 K=1" */
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Set
	for (size_t i=0; i<names.Size(); ++i)
	{
		     if (names[i]=="K" )  _K = values[i];
		else if (names[i]=="LA") _LA = values[i];
		else throw new Fatal("EmbSpring::SetModel: Parameter name (%s) is invalid",names[i].CStr());
	}
}

inline void EmbSpring::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
{
	// Allocate (local/element) displacements vector
	LinAlg::Vector<double> du(_nd*_n_nodes); // Delta disp. of this element

	// Assemble (local/element) displacements vector
	for (size_t i=0; i<_n_nodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		du(i*_nd+j) = dUglobal(_connects[i]->DOFVar(UD[_d][j]).EqID);

	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> df; // Delta internal force of this element

	LinAlg::Matrix<double> Ke;
	Order1Matrix(0,Ke);
	df = Ke * du;

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_n_nodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		dFint(_connects[i]->DOFVar(UD[_d][j]).EqID) += df(i*_nd+j);
}

inline void EmbSpring::CalcDepVars() const
{
	// Element displacements vector
	//_uL.Resize(_nd*_n_nodes);
	//for (size_t i=0; i<_n_nodes; ++i)
	//for (int    j=0; j<_nd;      ++j)
	//	_uL(i*_nd+j) = _connects[i]->DOFVar(UD[_d][j]).EssentialVal;

	// Transform to rod-local coordinates
	//LinAlg::Matrix<double> T;
	//_transf_mat(T);
	//_uL = T * _uL;
}

inline double EmbSpring::Val(int iNodeLocal, char const * Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).NaturalVal;

	//if (_uL.Size()<1) throw new Fatal("EmbSpring::Val: Please, call CalcDepVars() before calling this method");
	//double l = (iNodeLocal==0 ? 0 : 1.0);
	//     if (strcmp(Name,"N" )==0) return N(l);
	//else if (strcmp(Name,"Ea")==0) return    (_uL(_nd)-_uL(0))/_L;
	//else if (strcmp(Name,"Sa")==0) return _E*(_uL(_nd)-_uL(0))/_L;
	//else throw new Fatal("EmbSpring::Val: This element does not have a Val named %s",Name);
	return 0;
}

inline double EmbSpring::Val(char const * Name) const
{
	throw new Fatal("EmbSpring::Val: Feature not available");
}

inline void EmbSpring::Order1Matrix(size_t Index, LinAlg::Matrix<double> & Ke) const
{
	//              T   T                       T   T                       T   T
	//       K = [T0]*[B]*k0*[B]*[T0]*Area + [T1]*[B]*k1*[B]*[T1]*Area + [T2]*[B]*k2*[B]*[T2]*Area
	//        

	// Perpendicular stiffness
	double kp = 1.0E3*_K;

	Vector<double> shape;
	double x = _ext->X(), y = _ext->Y(), z = _ext->Z();
	double r, s, t;      // Local coordinates in tresspased element
	_tress->InverseMap(x, y, z, r, s, t);
	_tress->Shape(r, s, t, shape);

	// Mounting B Matrix
	Matrix<double> B(1,2); B = 1, -1;

	// Getting perpendicular vectors
	Vector<double> L0; 
	if (_ndim==2) { L0.Resize(2); L0 = _direct(0), _direct(1); }
	else L0 = _direct;

	Vector<double> AUX, L1, L2;

	_get_another_vector(L0, AUX);       // Find a vector (AUX) vector different from L0
	if (_ndim==2)
	{
		L1.Resize(2); L1 = -L0(1), L0(0); // Find a vector (L1)  perpendicular to L0

		// Mounting T Matrices (for all springs)
		Matrix<double> T0, T1;
		_mount_T_matrix(shape, L0, T0);
		_mount_T_matrix(shape, L1, T1);

		// Mounting Stiffness Matrix
		Ke = trn(T0)*trn(B)*_K*B*T0*_LA + trn(T1)*trn(B)*kp*B*T1*_LA;
	}
	else if (_ndim==3)
	{
		_get_normal_vector (L0, AUX, L1);   // Find a vector (L1)  perpendicular to L0 and AUX
		_get_normal_vector (L0, L1, L2);    // Find a vector (L2)  perpendicular to L0 and L1

		// Mounting T Matrices (for all springs)
		Matrix<double> T0, T1, T2;
		_mount_T_matrix(shape, L0, T0);
		_mount_T_matrix(shape, L1, T1);
		_mount_T_matrix(shape, L2, T2);

		// Perpendicular stiffness
		double kp = 1.0E4*_K;

		// Mounting Stiffness Matrix
		Ke = trn(T0)*trn(B)*_K*B*T0*_LA + trn(T1)*trn(B)*kp*B*T1*_LA + trn(T2)*trn(B)*kp*B*T2*_LA;
	}
}

inline void EmbSpring::B_Matrix(LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const
{
	throw new Fatal("EmbSpring::B_Matrix: Feature not available");
}

/* private */

inline void EmbSpring::_initialize()
{
	if (_ndim<1) throw new Fatal("EmbSpring::_initialize: For this element, _ndim must be greater than or equal to 1 (%d is invalid)",_ndim);
	_d  = _ndim-1; // Not used
	_nd = _ndim;   // Not used
	_nl = 0;
}

inline void EmbSpring::_calc_initial_internal_state()
{
	throw new Fatal("EmbSpring::_calc_initial_internal_state: Feature not available");
}

inline void EmbSpring::_mount_T_matrix(Vector<double> const & Shape, Vector<double> const & Direction, Matrix<double> & T) const 
{
	// Mounting T Matrix (2 x _n_nodes*_ndim)
	//
	//  T = [  N1*l N1*m N1*n N2*l N2*m ... | 0  0  0  ]
	//      [    0    0    0    0    0  ... | l  m  n  ]
	//      
	//  l, m, n = Direction(0), Direction(1), Direction(2)

	T.Resize(2, _ndim*_n_nodes);
	T.SetValues(0.0);
	for (size_t i=0; i<_n_nodes-1; i++)
	{
		              T(0,i*_ndim    ) = Shape(i)*Direction(0);
		              T(0,i*_ndim + 1) = Shape(i)*Direction(1);
		if (_ndim==3) T(0,i*_ndim + 2) = Shape(i)*Direction(2);
	}

	              T(1, (_n_nodes-1)*_ndim    ) = Direction(0);
	              T(1, (_n_nodes-1)*_ndim + 1) = Direction(1);
	if (_ndim==3) T(1, (_n_nodes-1)*_ndim + 2) = Direction(2);
}

inline void EmbSpring::_get_another_vector(Vector<double> const & V1, Vector<double> & V2) const 
{
	double sqrt2 = pow(2.0, 0.5); 
	double sqrt3 = pow(3.0, 0.5); 
	assert(V1.Size()==_ndim);
	V2.Resize(_ndim);

	if (_ndim==2)
	{
		if   (V1(0)!=0.0 && V1(1)!=0.0) { V2(0) = 1.0; V2(1) = 0.0; } // Taking a x direction vector 
		else { V2(0) = V2(1) = 1.0/sqrt2; }
	}
	else if (_ndim==3)
	{
		if   (V1(0)!=0.0 && V1(1)!=0.0 && V1(2)!=0.0) { V2(0) = 1.0; V2(1) = 0.0; V2(2) = 0.0; } // Taking a x direction vector 
		else { V2(0) = V2(1) = V2(2) = 1.0/sqrt3; }
	}
} 

inline void EmbSpring::_get_normal_vector(Vector<double> const & V1, Vector<double> const & V2, Vector<double> & V3) const 
{
	assert(V1.Size()==V2.Size() && _ndim==3);
	V3.Resize(_ndim);
	// Performing cross product
	V3(0) = V1(1)*V2(2) - V1(2)*V2(1);
	V3(1) = V1(2)*V2(0) - V1(0)*V2(2);
	V3(2) = V1(0)*V2(1) - V1(1)*V2(0);
	double norm = sqrt(pow(V3(0),2) + pow(V3(1),2) + pow(V3(2),2)); // Calculating the norm
	assert(norm != 0);
	V3 = V3/norm; // Normalizing the obtained normal vector
}

inline void EmbSpring::OutInfo(std::ostream & os) const
{
}

inline void EmbSpring::VTKConnect(String & Nodes) const 
{ 
	std::ostringstream os;
	for(size_t i=0; i<_connects.Size(); i++)
		os << _connects[i]->GetID() << " ";
	Nodes.Printf(os.str().c_str());
}

}; // namespace FEM

#endif // MECHSYS_FEM_EMBSPRING_H
