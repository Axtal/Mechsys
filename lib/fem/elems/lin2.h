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

#ifndef MECHSYS_FEM_LIN2_H
#define MECHSYS_FEM_LIN2_H

// MechSys
#include "fem/node.h"
#include "fem/element.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "vtkCellType.h"

namespace FEM
{

class Lin2: public virtual Element
{
public:
	// Constructor
	Lin2();

	// Destructor
	virtual ~Lin2() {}

	// Derived methods
	void SetIntPoints (int NumGaussPoints1D);
	int  VTKCellType  () const { return VTK_LINE; }
	void VTKConnect   (String & Nodes) const;
	void GetFaceNodes (int FaceID, Array<Node*> & FaceConnects) const { throw new Fatal("Lin2::GetFaceNodes: Method is not available for this element"); }
	void Shape        (double r, double s, double t, LinAlg::Vector<double> & Shape)  const;
	void Derivs       (double r, double s, double t, LinAlg::Matrix<double> & Derivs) const;
	void FaceShape    (double r, double s, LinAlg::Vector<double> & FaceShape)  const;
	void FaceDerivs   (double r, double s, LinAlg::Matrix<double> & FaceDerivs) const;

}; // class Lin2


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Lin2::Lin2()
{
	// Setup nodes number
	_n_nodes        = 2;
	_n_face_nodes   = 1;

	// Allocate nodes (connectivity)
	_connects.Resize    (_n_nodes);
	_connects.SetValues (NULL);

	// Integration points and Extrapolation Matrix
	SetIntPoints (/*NumGaussPoints1D*/2);
}

inline void Lin2::SetIntPoints(int NumGaussPoints1D)
{
	// Set IPs
	FEM::SetGaussIP (/*NDim*/1, NumGaussPoints1D, _a_int_pts);

	// Evaluation matrix [E]
	// [nodalvalues]T = [E] * [c0,c1]T
	//  nodalvalue_i  = c0 + c1*r_i
	// (i:node number)
	LinAlg::Matrix<double> eval_mat;
	eval_mat.Resize(_n_nodes, 2);
	eval_mat = 1.0, -1.0,
	           1.0,  1.0;

	// Set extrapolation matrix
	FEM::SetExtrapMatrix (/*NDim*/1, _a_int_pts, eval_mat, _extrap_mat);
}

inline void Lin2::VTKConnect(String & Nodes) const
{
	Nodes.Printf("%d %d",_connects[0]->GetID(),_connects[1]->GetID());
}

inline void Lin2::Shape(double r, double s, double t, LinAlg::Vector<double> & Shape) const
{
	/*  
	 *       0           |           1
	 *       @-----------+-----------@-> r
	 *      -1           |          +1
	 */
	Shape.Resize(/*NNodes*/2);
	Shape(0) = 0.5*(1.0-r);
	Shape(1) = 0.5*(1.0+r);
}

inline void Lin2::Derivs(double r, double s, double t, LinAlg::Matrix<double> & Derivs) const
{
	/*           _     _ T
	 *          |  dNi  |
	 * Derivs = |  ---  |   , where cj = r, s
	 *          |_ dcj _|
	 *
	 * Derivs(j,i), j=>local coordinate and i=>shape function
	 */
	Derivs.Resize(1,/*NNodes*/2);
	Derivs(0,0) = -0.5;
	Derivs(0,1) =  0.5;
}

inline void Lin2::FaceShape(double r, double s, LinAlg::Vector<double> & FaceShape) const
{
	throw new Fatal("Lin2::FaceShape: This method is not available for this element (Lin2)");
}

inline void Lin2::FaceDerivs(double r, double s, LinAlg::Matrix<double> & FaceDerivs) const
{
	throw new Fatal("Lin2::FaceDerivs: This method is not available for this element (Lin2)");
}


}; // namespace FEM

#endif // MECHSYS_FEM_LIN2_H
