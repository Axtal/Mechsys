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

#ifndef MECHSYS_FEM_HEX8_H
#define MECHSYS_FEM_HEX8_H

// MechSys
#include "fem/element.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"

namespace FEM
{

// Hex8 Constants
const int HEX8_NNODES      = 8;
const int HEX8_NINTPTS     = 8;
const int HEX8_NFACENODES  = 4;
const int HEX8_NFACEINTPTS = 4;
const Element::IntegPoint HEX8_INTPTS[]=
{{ -0.577350,  -0.577350,  -0.577350,  1. }, 
{   0.577350,  -0.577350,  -0.577350,  1. },
{   0.577350,   0.577350,  -0.577350,  1. },
{  -0.577350,   0.577350,  -0.577350,  1. },
{  -0.577350,  -0.577350,   0.577350,  1. },
{   0.577350,  -0.577350,   0.577350,  1. },
{   0.577350,   0.577350,   0.577350,  1. },
{  -0.577350,   0.577350,   0.577350,  1. }};
const Element::IntegPoint HEX8_FACEINTPTS[]=
{{ -0.577350,  -0.577350,  0.0,  1. }, 
{   0.577350,  -0.577350,  0.0,  1. },
{   0.577350,   0.577350,  0.0,  1. },
{  -0.577350,   0.577350,  0.0,  1. }};

class Hex8 : public virtual Element
{
public:
	// Constructor
	Hex8();

	// Destructor
	virtual ~Hex8() {}

	// Derived methods
	int    VTKCellType   () const { return 12; } // VTK_HEXAHEDRON
	void   Shape         (double r, double s, double t, LinAlg::Vector<double> & Shape)  const;
	void   Derivs        (double r, double s, double t, LinAlg::Matrix<double> & Derivs) const;
	void   FaceShape     (double r, double s, LinAlg::Vector<double> & FaceShape)  const;
	void   FaceDerivs    (double r, double s, LinAlg::Matrix<double> & FaceDerivs) const;
	double BoundDistance (double r, double s, double t) const;

}; // class Hex8


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Hex8::Hex8()
{
	// Setup nodes number
	_n_nodes        = HEX8_NNODES;
	_n_int_pts      = HEX8_NINTPTS;
	_n_face_nodes   = HEX8_NFACENODES;
	_n_face_int_pts = HEX8_NFACEINTPTS;

	// Allocate nodes (connectivity)
	_connects.Resize(_n_nodes);

	// Setup pointer to the array of Integration Points
	_a_int_pts      = HEX8_INTPTS;
	_a_face_int_pts = HEX8_FACEINTPTS;
}

inline void Hex8::Shape(double r, double s, double t, LinAlg::Vector<double> & Shape) const
{
	/*                    t
	 *                    ^
	 *                    |     
	 *                   4                7
	 *                    @________________@
	 *                  ,'|              ,'|
	 *                ,'  |            ,'  |
	 *              ,'    |          ,'    |
	 *        5   ,'      |     6  ,'      |
	 *          @'===============@'        |
	 *          |         |      |         |
	 *          |         |      |         |            
	 *          |       0 @_____ | ________@ 3 --> s
	 *          |       ,'       |       ,' 
	 *          |     ,'         |     ,' 
	 *          |   ,'           |   ,' 
	 *          | ,'             | ,' 
	 *          @________________@'
	 *         1                  2 
	 *      ,'
	 *    |_
	 *   r
	 */
	Shape.Resize(8);
	Shape(0) = 0.125*(1.0-r-s+r*s-t+s*t+r*t-r*s*t);
	Shape(1) = 0.125*(1.0+r-s-r*s-t+s*t-r*t+r*s*t);
	Shape(2) = 0.125*(1.0+r+s+r*s-t-s*t-r*t-r*s*t);
	Shape(3) = 0.125*(1.0-r+s-r*s-t-s*t+r*t+r*s*t);
	Shape(4) = 0.125*(1.0-r-s+r*s+t-s*t-r*t+r*s*t);
	Shape(5) = 0.125*(1.0+r-s-r*s+t-s*t+r*t-r*s*t);
	Shape(6) = 0.125*(1.0+r+s+r*s+t+s*t+r*t+r*s*t);
	Shape(7) = 0.125*(1.0-r+s-r*s+t+s*t-r*t-r*s*t);
}

inline void Hex8::Derivs(double r, double s, double t, LinAlg::Matrix<double> & Derivs) const
{
	/*           _     _ T
	 *          |  dNi  |
	 * Derivs = |  ---  |   , where cj = r, s
	 *          |_ dcj _|
	 *
	 * Derivs(j,i), j=>local coordinate and i=>shape function
	 */
	Derivs.Resize(3,8);
	Derivs(0,0) = 0.125*(-1.0+s+t-s*t);   Derivs(1,0)=0.125*(-1.0+r+t-r*t);   Derivs(2,0)=0.125*(-1.0+r+s-r*s);
	Derivs(0,1) = 0.125*(+1.0-s-t+s*t);   Derivs(1,1)=0.125*(-1.0-r+t+r*t);   Derivs(2,1)=0.125*(-1.0-r+s+r*s);
	Derivs(0,2) = 0.125*(+1.0+s-t-s*t);   Derivs(1,2)=0.125*(+1.0+r-t-r*t);   Derivs(2,2)=0.125*(-1.0-r-s-r*s);
	Derivs(0,3) = 0.125*(-1.0-s+t+s*t);   Derivs(1,3)=0.125*(+1.0-r-t+r*t);   Derivs(2,3)=0.125*(-1.0+r-s+r*s);
	Derivs(0,4) = 0.125*(-1.0+s-t+s*t);   Derivs(1,4)=0.125*(-1.0+r-t+r*t);   Derivs(2,4)=0.125*(+1.0-r-s+r*s);
	Derivs(0,5) = 0.125*(+1.0-s+t-s*t);   Derivs(1,5)=0.125*(-1.0-r-t-r*t);   Derivs(2,5)=0.125*(+1.0+r-s-r*s);
	Derivs(0,6) = 0.125*(+1.0+s+t+s*t);   Derivs(1,6)=0.125*(+1.0+r+t+r*t);   Derivs(2,6)=0.125*(+1.0+r+s+r*s);
	Derivs(0,7) = 0.125*(-1.0-s-t-s*t);   Derivs(1,7)=0.125*(+1.0-r+t-r*t);   Derivs(2,7)=0.125*(+1.0-r+s-r*s);
}

inline void Hex8::FaceShape(double r, double s, LinAlg::Vector<double> & FaceShape) const
{
	/*           s
	 *           ^
	 *           |             
	 *         3 @-----------@ 2
	 *           |           |
	 *           |           |
	 *           |           |
	 *           |           |
	 *           |           |
	 *           @-----------@-> r
	 *           0           1
	 */

	FaceShape.Resize(4);
	FaceShape(0) = 0.25*(1.0-r-s+r*s);
	FaceShape(1) = 0.25*(1.0+r-s-r*s);
	FaceShape(2) = 0.25*(1.0+r+s+r*s);
	FaceShape(3) = 0.25*(1.0-r+s-r*s);
}

inline void Hex8::FaceDerivs(double r, double s, LinAlg::Matrix<double> & FaceDerivs) const
{
	/*           _     _ T
	 *          |  dNi  |
	 * Derivs = |  ---  |   , where cj = r, s
	 *          |_ dcj _|
	 *
	 * Derivs(j,i), j=>local coordinate and i=>shape function
	 */

	FaceDerivs.Resize(2,4);
	FaceDerivs(0,0) = 0.25*(-1.0+s);   FaceDerivs(1,0) = 0.25*(-1.0+r);
	FaceDerivs(0,1) = 0.25*(+1.0-s);   FaceDerivs(1,1) = 0.25*(-1.0-r);
	FaceDerivs(0,2) = 0.25*(+1.0+s);   FaceDerivs(1,2) = 0.25*(+1.0+r);
	FaceDerivs(0,3) = 0.25*(-1.0-s);   FaceDerivs(1,3) = 0.25*(+1.0-r);
}

inline double Hex8::BoundDistance(double r, double s, double t) const
{
	return std::min(std::min( 1-fabs(r) , 1-fabs(s) ), 1-fabs(t)) ;
}


}; // namespace FEM

#endif // MECHSYS_FEM_HEX8_H
