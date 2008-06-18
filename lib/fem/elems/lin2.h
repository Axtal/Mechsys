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

namespace FEM
{

// Lin2 Constants
const int LIN2_NNODES      = 2;
const int LIN2_NINTPTS     = 2;
const int LIN2_NFACENODES  = 1;
const int LIN2_NFACEINTPTS = 1;
const Element::IntegPoint LIN2_INTPTS[]=
{{ -sqrt(3.0)/3.0, 0.0, 0.0, 1.0 },
 {  sqrt(3.0)/3.0, 0.0, 0.0, 1.0 }};
const Element::IntegPoint LIN2_FACEINTPTS[]=
{{ 1.0, 0.0, 0.0, 1.0 }};

class Lin2: public virtual Element
{
public:
	// Constructor
	Lin2();

	// Destructor
	virtual ~Lin2() {}

	// Derived methods
	int  VTKCellType () const { return 3; } // VTK_LINE
	void Shape       (double r, double s, double t, LinAlg::Vector<double> & Shape)  const;
	void Derivs      (double r, double s, double t, LinAlg::Matrix<double> & Derivs) const;
	void FaceShape   (double r, double s, LinAlg::Vector<double> & FaceShape)  const;
	void FaceDerivs  (double r, double s, LinAlg::Matrix<double> & FaceDerivs) const;

}; // class Lin2


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Lin2::Lin2()
{
	// Setup nodes number
	_n_dim          = 2;
	_n_nodes        = LIN2_NNODES;
	_n_int_pts      = LIN2_NINTPTS;
	_n_face_nodes   = LIN2_NFACENODES;
	_n_face_int_pts = LIN2_NFACEINTPTS;

	// Allocate nodes (connectivity)
	_connects.Resize(_n_nodes);

	// Setup pointer to the array of Integration Points
	_a_int_pts      = LIN2_INTPTS;
	_a_face_int_pts = LIN2_FACEINTPTS;
}

inline void Lin2::Shape(double r, double s, double t, LinAlg::Vector<double> & Shape) const
{
	/*  
	 *       0           |           1
	 *       @-----------+-----------@-> r
	 *      -1           |          +1
	 */
	Shape.Resize(LIN2_NNODES);
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
	Derivs.Resize(1,LIN2_NNODES);
	Derivs(0,0) = -0.5;
	Derivs(0,1) =  0.5;
}

inline void Lin2::FaceShape(double r, double s, LinAlg::Vector<double> & FaceShape) const
{
	/* Just a single point */
	FaceShape.Resize(LIN2_NFACENODES);
	FaceShape(0) = 1.0;
}

inline void Lin2::FaceDerivs(double r, double s, LinAlg::Matrix<double> & FaceDerivs) const
{
	/* Just a single point */
	FaceDerivs.Resize(1,LIN2_NFACENODES);
	FaceDerivs(0,0) = 0.0;
}


}; // namespace FEM

#endif // MECHSYS_FEM_LIN2_H
