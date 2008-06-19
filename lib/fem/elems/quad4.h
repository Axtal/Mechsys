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

#ifndef MECHSYS_FEM_QUAD4_H
#define MECHSYS_FEM_QUAD4_H

// MechSys
#include "fem/node.h"
#include "fem/element.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"

namespace FEM
{

// Quad4 Constants
const int QUAD4_NNODES      = 4;
const int QUAD4_NINTPTS     = 4;
const int QUAD4_NFACENODES  = 2;
const int QUAD4_NFACEINTPTS = 2;
const Element::IntegPoint QUAD4_INTPTS[]=
{{ -sqrt(3.0)/3.0, -sqrt(3.0)/3.0, 0.0, 1.0 },
 {  sqrt(3.0)/3.0, -sqrt(3.0)/3.0, 0.0, 1.0 },
 { -sqrt(3.0)/3.0,  sqrt(3.0)/3.0, 0.0, 1.0 },
 {  sqrt(3.0)/3.0,  sqrt(3.0)/3.0, 0.0, 1.0 }};
const Element::IntegPoint QUAD4_FACEINTPTS[]=
{{ -sqrt(3.0)/3.0, 0.0, 0.0, 1.0 },
 {  sqrt(3.0)/3.0, 0.0, 0.0, 1.0 }};
/*
const Element::IntegPoint QUAD4_INTPTS[]=
{{ -sqrt(3.0/5.0), -sqrt(3.0/5.0), 0.0, 25.0/81.0 },
 {           0.0 , -sqrt(3.0/5.0), 0.0, 40.0/81.0 },
 {  sqrt(3.0/5.0), -sqrt(3.0/5.0), 0.0, 25.0/81.0 },
 { -sqrt(3.0/5.0),           0.0 , 0.0, 40.0/81.0 },
 {           0.0 ,           0.0 , 0.0, 64.0/81.0 },
 {  sqrt(3.0/5.0),           0.0 , 0.0, 40.0/81.0 },
 { -sqrt(3.0/5.0),  sqrt(3.0/5.0), 0.0, 25.0/81.0 },
 {           0.0 ,  sqrt(3.0/5.0), 0.0, 40.0/81.0 },
 {  sqrt(3.0/5.0),  sqrt(3.0/5.0), 0.0, 25.0/81.0 }};
const Element::IntegPoint QUAD4_FACEINTPTS[]=
{{ -sqrt(3.0/5.0), 0.0, 0.0, 5.0/9.0 },
 {           0.0 , 0.0, 0.0, 8.0/9.0 },
 {  sqrt(3.0/5.0), 0.0, 0.0, 5.0/9.0 }};
*/

class Quad4: public virtual Element
{
public:
	// Constructor
	Quad4();

	// Destructor
	virtual ~Quad4() {}

	// Derived methods
	int  VTKCellType () const { return 9; } // VTK_QUAD
	void Shape       (double r, double s, double t, LinAlg::Vector<double> & Shape)  const;
	void Derivs      (double r, double s, double t, LinAlg::Matrix<double> & Derivs) const;
	void FaceShape   (double r, double s, LinAlg::Vector<double> & FaceShape)  const;
	void FaceDerivs  (double r, double s, LinAlg::Matrix<double> & FaceDerivs) const;

}; // class Quad4


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Quad4::Quad4()
{
	// Setup nodes number
	_n_nodes        = QUAD4_NNODES;
	_n_int_pts      = QUAD4_NINTPTS;
	_n_face_nodes   = QUAD4_NFACENODES;
	_n_face_int_pts = QUAD4_NFACEINTPTS;

	// Allocate nodes (connectivity)
	_connects.Resize(_n_nodes);

	// Setup pointer to the array of Integration Points
	_a_int_pts      = QUAD4_INTPTS;
	_a_face_int_pts = QUAD4_FACEINTPTS;
}

inline void Quad4::Shape(double r, double s, double t, LinAlg::Vector<double> & Shape) const
{
	/*      3                        2
	 *        @--------------------@
	 *        |               (1,1)|
	 *        |       s ^          |
	 *        |         |          |
	 *        |         |          |
	 *        |         +----> r   |
	 *        |       (0,0)        |
	 *        |                    |
	 *        |                    |
	 *        |(-1,-1)             |
	 *        @--------------------@
	 *      0                        1
	 */
	Shape.Resize (QUAD4_NNODES);

	double rp=1.0+r; double rm=1.0-r;
	double sp=1.0+s; double sm=1.0-s;

	Shape(0) = 0.25*rm*sm;
	Shape(1) = 0.25*rp*sm;
	Shape(2) = 0.25*rp*sp;
	Shape(3) = 0.25*rm*sp;
}

inline void Quad4::Derivs(double r, double s, double t, LinAlg::Matrix<double> & Derivs) const
{
	/*           _     _ T
	 *          |  dNi  |
	 * Derivs = |  ---  |   , where cj = r, s
	 *          |_ dcj _|
	 *
	 * Derivs(j,i), j=>local coordinate and i=>shape function
	 */
	Derivs.Resize (2, QUAD4_NNODES);

	double rp=1.0+r; double rm=1.0-r;
	double sp=1.0+s; double sm=1.0-s;

	Derivs(0,0) = -0.25*sm;   Derivs(1,0) = -0.25*rm;
	Derivs(0,1) =  0.25*sm;   Derivs(1,1) = -0.25*rp;
	Derivs(0,2) =  0.25*sp;   Derivs(1,2) =  0.25*rp;
	Derivs(0,3) = -0.25*sp;   Derivs(1,3) =  0.25*rm;
}

inline void Quad4::FaceShape(double r, double s, LinAlg::Vector<double> & FaceShape) const
{
	/*  
	 *       0           |           1
	 *       @-----------+-----------@-> r
	 *      -1           |          +1
	 */
	FaceShape.Resize(QUAD4_NFACENODES);
	FaceShape(0) = 0.5*(1.0-r);
	FaceShape(1) = 0.5*(1.0+r);
}

inline void Quad4::FaceDerivs(double r, double s, LinAlg::Matrix<double> & FaceDerivs) const
{
	/*           _     _ T
	 *          |  dNi  |
	 * Derivs = |  ---  |   , where cj = r, s
	 *          |_ dcj _|
	 *
	 * Derivs(j,i), j=>local coordinate and i=>shape function
	 */
	FaceDerivs.Resize(1,QUAD4_NFACENODES);
	FaceDerivs(0,0) = -0.5;
	FaceDerivs(0,1) =  0.5;
}


}; // namespace FEM

#endif // MECHSYS_FEM_QUAD4_H
