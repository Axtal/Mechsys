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

#ifndef MECHSYS_FEM_HEX20_H
#define MECHSYS_FEM_HEX20_H

// STL
#include <algorithm> 

// MechSys
#include "fem/element.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "vtkCellType.h"

namespace FEM
{

class Hex20 : public virtual Element
{
public:
	// Auxiliar structure to map local face IDs to local node IDs
	struct FaceMap
	{
		int n0; // node #0 local id
		int n1; // node #1 local id
		int n2; // node #2 local id
		int n3; // node #3 local id
		int n4; // node #4 local id
		int n5; // node #5 local id
		int n6; // node #6 local id
		int n7; // node #7 local id
	};
	static FaceMap Face2Node[];
	// Constructor
	Hex20();

	// Destructor
	virtual ~Hex20() {}

	// Derived methods
	void   SetIntPoints  (int NumGaussPoints1D);
	int    VTKCellType   () const { return VTK_QUADRATIC_HEXAHEDRON; }
	void   VTKConnect    (String & Nodes) const;
	void   GetFaceNodes  (int FaceID, Array<Node*> & FaceConnects) const;
	void   Shape         (double r, double s, double t, LinAlg::Vector<double> & Shape)  const;
	void   Derivs        (double r, double s, double t, LinAlg::Matrix<double> & Derivs) const;
	void   FaceShape     (double r, double s, LinAlg::Vector<double> & FaceShape)  const;
	void   FaceDerivs    (double r, double s, LinAlg::Matrix<double> & FaceDerivs) const;
	double BoundDistance (double r, double s, double t) const;
	void   LocalCoords   (LinAlg::Matrix<double> & coords) const;

}; // class Hex20

/* Local IDs
                  Vertices                            Edges                              Faces
    t
    |           4        15        7                   4      16
   ,+--s         @-------@--------@                 +-------+--------+                 +----------------+ 
 r'            ,'|              ,'|              6,'|            7 ,'|               ,'|              ,'| 
          12 @'  |         14 ,'  |             +'  |            ,'  |23           ,'  |  ___       ,'  | 
           ,'    |16        ,@    |19      18 ,'  20|          ,+    |           ,'    |,'5,'  [0],'    | 
     5   ,'      @      6 ,'      @         ,'      +  17    ,19     +         ,'      |~~~     ,'      | 
       @'=======@=======@'        |       +'=======+=======+'        |       +'===============+'  ,'|   | 
       |      13 |      |         |       |    5    |      |         |11     |   ,'|   |      |   |3|   | 
       |         |      |  11     |     21|        8|  0   |22  12   |       |   |2|   |      |   |,'   | 
    17 |       0 @- - - | @- - - -@       |         +- - - | +- - - -+       |   |,'   +- - - | +- - - -+ 
       @       ,'       @       ,' 3      +      2,'       +       ,'        |       ,'       |       ,'  
       |   8 @'      18 |     ,'          |     +'         |     ,'3         |     ,' [1]  ___|     ,'    
       |   ,'           |   ,@ 10        9|   14         10|   ,+            |   ,'      ,'4,'|   ,'      
       | ,'             | ,'              | ,'             | ,15             | ,'        ~~~  | ,'        
       @-------@--------@'                +-------+--------+'                +----------------+'          
     1         9         2                    1       13
*/
Hex20::FaceMap Hex20::Face2Node[]= {{ 0, 3, 7, 4, 11, 19, 15, 16 },
                                    { 1, 2, 6, 5,  9, 18, 13, 17 },
                                    { 0, 1, 5, 4,  8, 17, 12, 16 },
                                    { 2, 3, 7, 6, 10, 19, 14, 18 }, 
                                    { 0, 1, 2, 3,  8,  9, 10, 11 }, 
                                    { 4, 5, 6, 7, 12, 13, 14, 15 }};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Hex20::Hex20()
{
	// Setup nodes number
	_n_nodes        = 20;
	_n_face_nodes   = 8;

	// Allocate nodes (connectivity)
	_connects.Resize(_n_nodes);
	_connects.SetValues(NULL);

	// Integration Points and Extrapolation Matrix
	SetIntPoints (/*NumGaussPoints1D*/2);
}

inline void Hex20::SetIntPoints(int NumGaussPoints1D)
{
	// Setup pointer to the array of Integration Points
	if      (NumGaussPoints1D==2) _a_int_pts = HEX_IP2;
	else if (NumGaussPoints1D==3) _a_int_pts = HEX_IP3;
	else if (NumGaussPoints1D==4) _a_int_pts = HEX_IP4;
	else if (NumGaussPoints1D==5) _a_int_pts = HEX_IP5;
	else throw new Fatal("Hex20::SetIntPoints: Error in number of integration points.");

	_n_int_pts      = pow(NumGaussPoints1D, 3);
	_a_face_int_pts = QUAD_IP2;
	_n_face_int_pts = 4;
}

inline void Hex20::LocalCoords(LinAlg::Matrix<double> & coords) const 
{
	coords.Resize(20,4);
	coords = -1.0, -1.0, -1.0, 1.0, // nodes 0 to 7
	          1.0, -1.0, -1.0, 1.0,
	          1.0,  1.0, -1.0, 1.0,
	         -1.0,  1.0, -1.0, 1.0,
	         -1.0, -1.0,  1.0, 1.0,
	          1.0, -1.0,  1.0, 1.0,
	          1.0,  1.0,  1.0, 1.0,
	         -1.0,  1.0,  1.0, 1.0,

              0.0, -1.0, -1.0, 1.0, // nodes 8 to 11
              1.0,  0.0, -1.0, 1.0,
              0.0,  1.0, -1.0, 1.0,
             -1.0,  0.0, -1.0, 1.0,

              0.0, -1.0,  1.0, 1.0, // nodes 12 to 15
              1.0,  0.0,  1.0, 1.0,
              0.0,  1.0,  1.0, 1.0,
             -1.0,  0.0,  1.0, 1.0,

	         -1.0, -1.0,  0.0, 1.0, // nodes 16 to 19
	          1.0, -1.0,  0.0, 1.0,
	          1.0,  1.0,  0.0, 1.0,
	         -1.0,  1.0,  0.0, 1.0;
}

inline void Hex20::VTKConnect(String & Nodes) const
{
	Nodes.Printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
	             _connects[ 0]->GetID(),
	             _connects[ 1]->GetID(),
	             _connects[ 2]->GetID(),
	             _connects[ 3]->GetID(),
	             _connects[ 4]->GetID(),
	             _connects[ 5]->GetID(),
	             _connects[ 6]->GetID(),
	             _connects[ 7]->GetID(),
	             _connects[ 8]->GetID(),
	             _connects[ 9]->GetID(),
	             _connects[10]->GetID(),
	             _connects[11]->GetID(),
	             _connects[12]->GetID(),
	             _connects[13]->GetID(),
	             _connects[14]->GetID(),
	             _connects[15]->GetID(),
	             _connects[16]->GetID(),
	             _connects[17]->GetID(),
	             _connects[18]->GetID(),
	             _connects[19]->GetID());
}

inline void Hex20::GetFaceNodes(int FaceID, Array<Node*> & FaceConnects) const
{
	FaceConnects.Resize(/*NumFaceNodes*/8);
	FaceConnects[0] = _connects[Face2Node[FaceID].n0];
	FaceConnects[1] = _connects[Face2Node[FaceID].n1];
	FaceConnects[2] = _connects[Face2Node[FaceID].n2];
	FaceConnects[3] = _connects[Face2Node[FaceID].n3];
	FaceConnects[4] = _connects[Face2Node[FaceID].n4];
	FaceConnects[5] = _connects[Face2Node[FaceID].n5];
	FaceConnects[6] = _connects[Face2Node[FaceID].n6];
	FaceConnects[7] = _connects[Face2Node[FaceID].n7];
}

inline void Hex20::Shape(double r, double s, double t, LinAlg::Vector<double> & Shape) const
{
	Shape.Resize(20);

	double rp1=1.0+r; double rm1=1.0-r;
	double sp1=1.0+s; double sm1=1.0-s;
	double tp1=1.0+t; double tm1=1.0-t;

	Shape( 0) = 1/8.*rm1*sm1*tm1*(-r-s-t-2);
	Shape( 1) = 1/8.*rp1*sm1*tm1*( r-s-t-2);
	Shape( 2) = 1/8.*rp1*sp1*tm1*( r+s-t-2);
	Shape( 3) = 1/8.*rm1*sp1*tm1*(-r+s-t-2);
	Shape( 4) = 1/8.*rm1*sm1*tp1*(-r-s+t-2);
	Shape( 5) = 1/8.*rp1*sm1*tp1*( r-s+t-2);
	Shape( 6) = 1/8.*rp1*sp1*tp1*( r+s+t-2);
	Shape( 7) = 1/8.*rm1*sp1*tp1*(-r+s+t-2);
	Shape( 8) = 1/4.*(1-r*r)*sm1*tm1;
	Shape( 9) = 1/4.*rp1*(1-s*s)*tm1;
	Shape(10) = 1/4.*(1-r*r)*sp1*tm1;
	Shape(11) = 1/4.*rm1*(1-s*s)*tm1;
	Shape(12) = 1/4.*(1-r*r)*sm1*tp1;
	Shape(13) = 1/4.*rp1*(1-s*s)*tp1;
	Shape(14) = 1/4.*(1-r*r)*sp1*tp1;
	Shape(15) = 1/4.*rm1*(1-s*s)*tp1;
	Shape(16) = 1/4.*rm1*sm1*(1-t*t);
	Shape(17) = 1/4.*rp1*sm1*(1-t*t);
	Shape(18) = 1/4.*rp1*sp1*(1-t*t);
	Shape(19) = 1/4.*rm1*sp1*(1-t*t);
}

inline void Hex20::Derivs(double r, double s, double t, LinAlg::Matrix<double> & Derivs) const
{
	/*           _     _ T
	 *          |  dNi  |
	 * Derivs = |  ---  |   , where cj = r, s
	 *          |_ dcj _|
	 *
	 * Derivs(j,i), j=>local coordinate and i=>shape function
	 */

	Derivs.Resize(3,20);

	double rp1=1.0+r; double rm1=1.0-r;
	double sp1=1.0+s; double sm1=1.0-s;
	double tp1=1.0+t; double tm1=1.0-t;
	
	//Derivatives with respect to r
    Derivs(0, 0)= -.125*sm1*tm1*(-r-s-t-2)-0.125*rm1*sm1*tm1;
	Derivs(0, 1)=  .125*sm1*tm1*( r-s-t-2)+0.125*rp1*sm1*tm1;
	Derivs(0, 2)=  .125*sp1*tm1*( r+s-t-2)+0.125*rp1*sp1*tm1;
	Derivs(0, 3)= -.125*sp1*tm1*(-r+s-t-2)-0.125*rm1*sp1*tm1;
	Derivs(0, 4)= -.125*sm1*tp1*(-r-s+t-2)-0.125*rm1*sm1*tp1;
	Derivs(0, 5)=  .125*sm1*tp1*( r-s+t-2)+0.125*rp1*sm1*tp1;
	Derivs(0, 6)=  .125*sp1*tp1*( r+s+t-2)+0.125*rp1*sp1*tp1;
	Derivs(0, 7)= -.125*sp1*tp1*(-r+s+t-2)-0.125*rm1*sp1*tp1;
	Derivs(0, 8)= -.5*r*sm1*tm1;
	Derivs(0, 9)=  .25*(1-s*s)*tm1;
	Derivs(0,10)= -.5*r*sp1*tm1;
	Derivs(0,11)= -.25*(1-s*s)*tm1;
	Derivs(0,12)= -.5*r*sm1*tp1;
	Derivs(0,13)=  .25*(1-s*s)*tp1;
	Derivs(0,14)= -.5*r*sp1  *tp1;
	Derivs(0,15)= -.25*(1-s*s)*tp1;
	Derivs(0,16)= -.25*sm1*(1-t*t);
	Derivs(0,17)=  .25*sm1*(1-t*t);
	Derivs(0,18)=  .25*sp1*(1-t*t);
	Derivs(0,19)= -.25*sp1*(1-t*t);


	//Derivatives with respect to s
	Derivs(1, 0)= -.125*rm1*tm1*(-r-s-t-2)-0.125*rm1*sm1*tm1;
	Derivs(1, 1)= -.125*rp1*tm1*( r-s-t-2)-0.125*rp1*sm1*tm1;
	Derivs(1, 2)=  .125*rp1*tm1*( r+s-t-2)+0.125*rp1*sp1*tm1;
	Derivs(1, 3)=  .125*rm1*tm1*(-r+s-t-2)+0.125*rm1*sp1*tm1;
	Derivs(1, 4)= -.125*rm1*tp1*(-r-s+t-2)-0.125*rm1*sm1*tp1;
	Derivs(1, 5)= -.125*rp1*tp1*( r-s+t-2)-0.125*rp1*sm1*tp1;
	Derivs(1, 6)=  .125*rp1*tp1*( r+s+t-2)+0.125*rp1*sp1*tp1;
	Derivs(1, 7)=  .125*rm1*tp1*(-r+s+t-2)+0.125*rm1*sp1*tp1;
	Derivs(1, 8)= -.25*(1-r*r)*tm1;
	Derivs(1, 9)= -.5*s*rp1*tm1;
	Derivs(1,10)=  .25*(1-r*r)*tm1;
	Derivs(1,11)= -.5*s*rm1*tm1;
	Derivs(1,12)= -.25*(1-r*r)*tp1;
	Derivs(1,13)= -.5*s*rp1*tp1;
	Derivs(1,14)=  .25*(1-r*r)*tp1;
	Derivs(1,15)= -.5*s*rm1*tp1;
	Derivs(1,16)= -.25*rm1*(1-t*t);
	Derivs(1,17)= -.25*rp1*(1-t*t);
	Derivs(1,18)=  .25*rp1*(1-t*t);
	Derivs(1,19)=  .25*rm1*(1-t*t);

	//Derivatives with respect to t
	Derivs(2, 0)= -.125*rm1*sm1*(-r-s-t-2)-0.125*rm1*sm1*tm1;
	Derivs(2, 1)= -.125*rp1*sm1*( r-s-t-2)-0.125*rp1*sm1*tm1;
	Derivs(2, 2)= -.125*rp1*sp1*( r+s-t-2)-0.125*rp1*sp1*tm1;
	Derivs(2, 3)= -.125*rm1*sp1*(-r+s-t-2)-0.125*rm1*sp1*tm1;
	Derivs(2, 4)=  .125*rm1*sm1*(-r-s+t-2)+0.125*rm1*sm1*tp1;
	Derivs(2, 5)=  .125*rp1*sm1*( r-s+t-2)+0.125*rp1*sm1*tp1;
	Derivs(2, 6)=  .125*rp1*sp1*( r+s+t-2)+0.125*rp1*sp1*tp1;
	Derivs(2, 7)=  .125*rm1*sp1*(-r+s+t-2)+0.125*rm1*sp1*tp1;
	Derivs(2, 8)= -.25*(1-r*r)*sm1;
	Derivs(2, 9)= -.25*rp1*(1-s*s);
	Derivs(2,10)= -.25*(1-r*r)*sp1;
	Derivs(2,11)= -.25*rm1*(1-s*s);
	Derivs(2,12)=  .25*(1-r*r)*sm1;
	Derivs(2,13)=  .25*rp1*(1-s*s);
	Derivs(2,14)=  .25*(1-r*r)*sp1;
	Derivs(2,15)=  .25*rm1*(1-s*s);
	Derivs(2,16)= -.5*t*rm1*sm1;
	Derivs(2,17)= -.5*t*rp1*sm1;
	Derivs(2,18)= -.5*t*rp1*sp1;
	Derivs(2,19)= -.5*t*rm1*sp1;
}
	
inline void Hex20::FaceShape(double r, double s, LinAlg::Vector<double> & FaceShape) const
{
	/*                 ^ s
	 *                 |
	 *                 6       
	 *         3 @-----------@ 2
	 *           |           |
	 *           |           |
	 *         7 |           | 5  -> r
	 *           |           |
	 *           |           |
	 *           @-----------@
	 *           0     4     1
	 */
	FaceShape.Resize(8);

	double rp1=1.0+r; double rm1=1.0-r;
	double sp1=1.0+s; double sm1=1.0-s;

	FaceShape(0) = 0.25*rm1*sm1*(rm1+sm1-3.0);
	FaceShape(1) = 0.25*rp1*sm1*(rp1+sm1-3.0);
	FaceShape(2) = 0.25*rp1*sp1*(rp1+sp1-3.0);
	FaceShape(3) = 0.25*rm1*sp1*(rm1+sp1-3.0);
	FaceShape(4) = 0.50*sm1*(1.0-r*r);
	FaceShape(5) = 0.50*rp1*(1.0-s*s);
	FaceShape(6) = 0.50*sp1*(1.0-r*r);
	FaceShape(7) = 0.50*rm1*(1.0-s*s);
}

inline void Hex20::FaceDerivs(double r, double s, LinAlg::Matrix<double> & FaceDerivs) const
{
	/*           _     _ T
	 *          |  dNi  |
	 * Derivs = |  ---  |   , where cj = r, s
	 *          |_ dcj _|
	 *
	 * Derivs(j,i), j=>local coordinate and i=>shape function
	 */
	FaceDerivs.Resize(2,8);
	double rp1=1.0+r; double rm1=1.0-r;
	double sp1=1.0+s; double sm1=1.0-s;

	FaceDerivs(0,0) = - 0.25 * sm1 * (rm1 + rm1 + sm1 - 3.0);
	FaceDerivs(0,1) =   0.25 * sm1 * (rp1 + rp1 + sm1 - 3.0);
	FaceDerivs(0,2) =   0.25 * sp1 * (rp1 + rp1 + sp1 - 3.0);
	FaceDerivs(0,3) = - 0.25 * sp1 * (rm1 + rm1 + sp1 - 3.0);
	FaceDerivs(0,4) = - r * sm1;
	FaceDerivs(0,5) =   0.50 * (1.0 - s * s);
	FaceDerivs(0,6) = - r * sp1;
	FaceDerivs(0,7) = - 0.5 * (1.0 - s * s);

	FaceDerivs(1,0) = - 0.25 * rm1 * (sm1 + rm1 + sm1 - 3.0);
	FaceDerivs(1,1) = - 0.25 * rp1 * (sm1 + rp1 + sm1 - 3.0);
	FaceDerivs(1,2) =   0.25 * rp1 * (sp1 + rp1 + sp1 - 3.0);
	FaceDerivs(1,3) =   0.25 * rm1 * (sp1 + rm1 + sp1 - 3.0);
	FaceDerivs(1,4) = - 0.50 * (1.0 - r * r);
	FaceDerivs(1,5) = - s * rp1;
	FaceDerivs(1,6) =   0.50 * (1.0 - r * r);
	FaceDerivs(1,7) = - s * rm1;
}

inline double Hex20::BoundDistance(double r, double s, double t) const
{
	return std::min(std::min( 1-fabs(r) , 1-fabs(s) ), 1-fabs(t)) ;
}


}; // namespace FEM

#endif // MECHSYS_FEM_HEX20_H
