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
#include "fem/node.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "vtkCellType.h"

namespace FEM
{

class Hex20 : public GeomElem
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
	void   SetIPs     (int NIPs1D);
	int    VTKType    () const { return VTK_QUADRATIC_HEXAHEDRON; }
	void   VTKConn    (String & Nodes) const;
	void   GetFNodes  (int FaceID, Array<Node*> & FaceConnects) const;
	double BoundDist  (double r, double s, double t) const { return std::min(std::min(1-fabs(r),1-fabs(s)),1-fabs(t)); }
	void   Shape      (double r, double s, double t, Vec_t & N)  const;
	void   Derivs     (double r, double s, double t, Mat_t & dN) const;
	void   FaceShape  (double r, double s, Vec_t & FN)  const;
	void   FaceDerivs (double r, double s, Mat_t & FdN) const;

private:
	void _local_coords (Mat_t & coords) const;

}; // class Hex20


/* Local IDs
                  Vertices                             Edges                              Faces
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
Hex20::FaceMap Hex20::Face2Node[]= {{ 0, 4, 7, 3, 16, 15, 19, 11 },
                                    { 1, 2, 6, 5,  9, 18, 13, 17 },
                                    { 0, 1, 5, 4,  8, 17, 12, 16 },
                                    { 2, 3, 7, 6, 10, 19, 14, 18 }, 
                                    { 0, 3, 2, 1, 11, 10,  9,  8 }, 
                                    { 4, 5, 6, 7, 12, 13, 14, 15 }};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Hex20::Hex20()
{
	// Setup nodes number
	NNodes  = 20;
	NFNodes = 8;

	// Allocate nodes (connectivity)
	Conn.Resize    (NNodes);
	Conn.SetValues (NULL);

	// Integration points and Extrapolation Matrix
	SetIPs (/*NIPs1D*/2);
}

inline void Hex20::SetIPs(int NumGaussPoints1D)
{
	// Setup pointer to the array of Integration Points
	     if (NIPs1D==2) IPs = HEX_IP2;
	else if (NIPs1D==3) IPs = HEX_IP3;
	else if (NIPs1D==4) IPs = HEX_IP4;
	else if (NIPs1D==5) IPs = HEX_IP5;
	else throw new Fatal("Hex20::SetIPs: Number of integration points < %d > is invalid",NIPs1D);

	NIPs  = pow(NIPs1D, 3);
	FIPs  = QUAD_IP2;
	NFIPs = 4; 
}


inline void Hex20::VTKConnect(String & Nodes) const
{
	Nodes.Printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
	             Conn[ 0]->GetID(),
	             Conn[ 1]->GetID(),
	             Conn[ 2]->GetID(),
	             Conn[ 3]->GetID(),
	             Conn[ 4]->GetID(),
	             Conn[ 5]->GetID(),
	             Conn[ 6]->GetID(),
	             Conn[ 7]->GetID(),
	             Conn[ 8]->GetID(),
	             Conn[ 9]->GetID(),
	             Conn[10]->GetID(),
	             Conn[11]->GetID(),
	             Conn[12]->GetID(),
	             Conn[13]->GetID(),
	             Conn[14]->GetID(),
	             Conn[15]->GetID(),
	             Conn[16]->GetID(),
	             Conn[17]->GetID(),
	             Conn[18]->GetID(),
	             Conn[19]->GetID());
}

inline void Hex20::GetFaceNodes(int FaceID, Array<Node*> & FaceConnects) const
{
	FaceConnects.Resize(/*NumFaceNodes*/8);
	FaceConnects[0] = Conn[Face2Node[FaceID].n0];
	FaceConnects[1] = Conn[Face2Node[FaceID].n1];
	FaceConnects[2] = Conn[Face2Node[FaceID].n2];
	FaceConnects[3] = Conn[Face2Node[FaceID].n3];
	FaceConnects[4] = Conn[Face2Node[FaceID].n4];
	FaceConnects[5] = Conn[Face2Node[FaceID].n5];
	FaceConnects[6] = Conn[Face2Node[FaceID].n6];
	FaceConnects[7] = Conn[Face2Node[FaceID].n7];
}

inline void Hex20::Shape(double r, double s, double t, Vec_t & N) const
{
	N.Resize(20);

	double rp1=1.0+r; double rm1=1.0-r;
	double sp1=1.0+s; double sm1=1.0-s;
	double tp1=1.0+t; double tm1=1.0-t;

	N( 0) = 1/8.*rm1*sm1*tm1*(-r-s-t-2);
	N( 1) = 1/8.*rp1*sm1*tm1*( r-s-t-2);
	N( 2) = 1/8.*rp1*sp1*tm1*( r+s-t-2);
	N( 3) = 1/8.*rm1*sp1*tm1*(-r+s-t-2);
	N( 4) = 1/8.*rm1*sm1*tp1*(-r-s+t-2);
	N( 5) = 1/8.*rp1*sm1*tp1*( r-s+t-2);
	N( 6) = 1/8.*rp1*sp1*tp1*( r+s+t-2);
	N( 7) = 1/8.*rm1*sp1*tp1*(-r+s+t-2);
	N( 8) = 1/4.*(1-r*r)*sm1*tm1;
	N( 9) = 1/4.*rp1*(1-s*s)*tm1;
	N(10) = 1/4.*(1-r*r)*sp1*tm1;
	N(11) = 1/4.*rm1*(1-s*s)*tm1;
	N(12) = 1/4.*(1-r*r)*sm1*tp1;
	N(13) = 1/4.*rp1*(1-s*s)*tp1;
	N(14) = 1/4.*(1-r*r)*sp1*tp1;
	N(15) = 1/4.*rm1*(1-s*s)*tp1;
	N(16) = 1/4.*rm1*sm1*(1-t*t);
	N(17) = 1/4.*rp1*sm1*(1-t*t);
	N(18) = 1/4.*rp1*sp1*(1-t*t);
	N(19) = 1/4.*rm1*sp1*(1-t*t);
}

inline void Hex20::Derivs(double r, double s, double t, Mat_t & dN) const
{
	/*       _     _ T
	 *      |  dNi  |
	 * dN = |  ---  |   , where cj = r, s
	 *      |_ dcj _|
	 *
	 * dN(j,i), j=>local coordinate and i=>shape function
	 */

	dN.Resize(3,20);

	double rp1=1.0+r; double rm1=1.0-r;
	double sp1=1.0+s; double sm1=1.0-s;
	double tp1=1.0+t; double tm1=1.0-t;
	
	//Derivatives with respect to r
	dN(0, 0)= -.125*sm1*tm1*(-r-s-t-2)-0.125*rm1*sm1*tm1;
	dN(0, 1)=  .125*sm1*tm1*( r-s-t-2)+0.125*rp1*sm1*tm1;
	dN(0, 2)=  .125*sp1*tm1*( r+s-t-2)+0.125*rp1*sp1*tm1;
	dN(0, 3)= -.125*sp1*tm1*(-r+s-t-2)-0.125*rm1*sp1*tm1;
	dN(0, 4)= -.125*sm1*tp1*(-r-s+t-2)-0.125*rm1*sm1*tp1;
	dN(0, 5)=  .125*sm1*tp1*( r-s+t-2)+0.125*rp1*sm1*tp1;
	dN(0, 6)=  .125*sp1*tp1*( r+s+t-2)+0.125*rp1*sp1*tp1;
	dN(0, 7)= -.125*sp1*tp1*(-r+s+t-2)-0.125*rm1*sp1*tp1;
	dN(0, 8)= -.5*r*sm1*tm1;
	dN(0, 9)=  .25*(1-s*s)*tm1;
	dN(0,10)= -.5*r*sp1*tm1;
	dN(0,11)= -.25*(1-s*s)*tm1;
	dN(0,12)= -.5*r*sm1*tp1;
	dN(0,13)=  .25*(1-s*s)*tp1;
	dN(0,14)= -.5*r*sp1  *tp1;
	dN(0,15)= -.25*(1-s*s)*tp1;
	dN(0,16)= -.25*sm1*(1-t*t);
	dN(0,17)=  .25*sm1*(1-t*t);
	dN(0,18)=  .25*sp1*(1-t*t);
	dN(0,19)= -.25*sp1*(1-t*t);

	//Derivatives with respect to s
	dN(1, 0)= -.125*rm1*tm1*(-r-s-t-2)-0.125*rm1*sm1*tm1;
	dN(1, 1)= -.125*rp1*tm1*( r-s-t-2)-0.125*rp1*sm1*tm1;
	dN(1, 2)=  .125*rp1*tm1*( r+s-t-2)+0.125*rp1*sp1*tm1;
	dN(1, 3)=  .125*rm1*tm1*(-r+s-t-2)+0.125*rm1*sp1*tm1;
	dN(1, 4)= -.125*rm1*tp1*(-r-s+t-2)-0.125*rm1*sm1*tp1;
	dN(1, 5)= -.125*rp1*tp1*( r-s+t-2)-0.125*rp1*sm1*tp1;
	dN(1, 6)=  .125*rp1*tp1*( r+s+t-2)+0.125*rp1*sp1*tp1;
	dN(1, 7)=  .125*rm1*tp1*(-r+s+t-2)+0.125*rm1*sp1*tp1;
	dN(1, 8)= -.25*(1-r*r)*tm1;
	dN(1, 9)= -.5*s*rp1*tm1;
	dN(1,10)=  .25*(1-r*r)*tm1;
	dN(1,11)= -.5*s*rm1*tm1;
	dN(1,12)= -.25*(1-r*r)*tp1;
	dN(1,13)= -.5*s*rp1*tp1;
	dN(1,14)=  .25*(1-r*r)*tp1;
	dN(1,15)= -.5*s*rm1*tp1;
	dN(1,16)= -.25*rm1*(1-t*t);
	dN(1,17)= -.25*rp1*(1-t*t);
	dN(1,18)=  .25*rp1*(1-t*t);
	dN(1,19)=  .25*rm1*(1-t*t);

	//Derivatives with respect to t
	dN(2, 0)= -.125*rm1*sm1*(-r-s-t-2)-0.125*rm1*sm1*tm1;
	dN(2, 1)= -.125*rp1*sm1*( r-s-t-2)-0.125*rp1*sm1*tm1;
	dN(2, 2)= -.125*rp1*sp1*( r+s-t-2)-0.125*rp1*sp1*tm1;
	dN(2, 3)= -.125*rm1*sp1*(-r+s-t-2)-0.125*rm1*sp1*tm1;
	dN(2, 4)=  .125*rm1*sm1*(-r-s+t-2)+0.125*rm1*sm1*tp1;
	dN(2, 5)=  .125*rp1*sm1*( r-s+t-2)+0.125*rp1*sm1*tp1;
	dN(2, 6)=  .125*rp1*sp1*( r+s+t-2)+0.125*rp1*sp1*tp1;
	dN(2, 7)=  .125*rm1*sp1*(-r+s+t-2)+0.125*rm1*sp1*tp1;
	dN(2, 8)= -.25*(1-r*r)*sm1;
	dN(2, 9)= -.25*rp1*(1-s*s);
	dN(2,10)= -.25*(1-r*r)*sp1;
	dN(2,11)= -.25*rm1*(1-s*s);
	dN(2,12)=  .25*(1-r*r)*sm1;
	dN(2,13)=  .25*rp1*(1-s*s);
	dN(2,14)=  .25*(1-r*r)*sp1;
	dN(2,15)=  .25*rm1*(1-s*s);
	dN(2,16)= -.5*t*rm1*sm1;
	dN(2,17)= -.5*t*rp1*sm1;
	dN(2,18)= -.5*t*rp1*sp1;
	dN(2,19)= -.5*t*rm1*sp1;
}

inline void Hex20::FaceShape(double r, double s, Vec_t & FaceShape) const
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

	FN(0) = 0.25*rm1*sm1*(rm1+sm1-3.0);
	FN(1) = 0.25*rp1*sm1*(rp1+sm1-3.0);
	FN(2) = 0.25*rp1*sp1*(rp1+sp1-3.0);
	FN(3) = 0.25*rm1*sp1*(rm1+sp1-3.0);
	FN(4) = 0.50*sm1*(1.0-r*r);
	FN(5) = 0.50*rp1*(1.0-s*s);
	FN(6) = 0.50*sp1*(1.0-r*r);
	FN(7) = 0.50*rm1*(1.0-s*s);
}

inline void Hex20::FaceDerivs(double r, double s, Mat_t & FdN) const
{
	/*          _     _ T
	 *         |  dNi  |
	 *   FdN = |  ---  |   , where cj = r, s
	 *         |_ dcj _|
	 *
	 *   FdN(j,i), j=>local coordinate and i=>shape function
	 */
	FdN.Resize(2,8);
	double rp1=1.0+r; double rm1=1.0-r;
	double sp1=1.0+s; double sm1=1.0-s;

	FdN(0,0) = - 0.25 * sm1 * (rm1 + rm1 + sm1 - 3.0);
	FdN(0,1) =   0.25 * sm1 * (rp1 + rp1 + sm1 - 3.0);
	FdN(0,2) =   0.25 * sp1 * (rp1 + rp1 + sp1 - 3.0);
	FdN(0,3) = - 0.25 * sp1 * (rm1 + rm1 + sp1 - 3.0);
	FdN(0,4) = - r * sm1;
	FdN(0,5) =   0.50 * (1.0 - s * s);
	FdN(0,6) = - r * sp1;
	FdN(0,7) = - 0.5 * (1.0 - s * s);

	FdN(1,0) = - 0.25 * rm1 * (sm1 + rm1 + sm1 - 3.0);
	FdN(1,1) = - 0.25 * rp1 * (sm1 + rp1 + sm1 - 3.0);
	FdN(1,2) =   0.25 * rp1 * (sp1 + rp1 + sp1 - 3.0);
	FdN(1,3) =   0.25 * rm1 * (sp1 + rm1 + sp1 - 3.0);
	FdN(1,4) = - 0.50 * (1.0 - r * r);
	FdN(1,5) = - s * rp1;
	FdN(1,6) =   0.50 * (1.0 - r * r);
	FdN(1,7) = - s * rm1;
}

inline void Hex20::_local_coords(Mat_t & coords) const 
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


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
GeomElem * Hex20Maker() { return new Hex20(); }

// Register element
int Hex20Register() { GeomElemFactory["Hex20"]=Hex20Maker; return 0; }

// Call register
int __Hex20_dummy_int  = Hex20Register();

}; // namespace FEM

#endif // MECHSYS_FEM_HEX20_H
