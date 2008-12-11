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
#include "fem/node.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "vtkCellType.h"

namespace FEM
{

class Hex8 : public GeomElem
{
public:
	// Auxiliar structure to map local face IDs to local node IDs
	struct FaceMap
	{
		int n0; // node #0 local id
		int n1; // node #1 local id
		int n2; // node #2 local id
		int n3; // node #3 local id
	};
	static FaceMap Face2Node[];

	// Constructor
	Hex8();

	// Derived methods
	void   SetIPs     (int NIPs1D);
	int    VTKType    () const { return VTK_HEXAHEDRON; }
	void   VTKConn    (String & Nodes) const;
	void   GetFNodes  (int FaceID, Array<Node*> & FaceConnects) const;
	double BoundDist  (double r, double s, double t) const { return std::min(std::min(1-fabs(r),1-fabs(s)),1-fabs(t)); }
	void   Shape      (double r, double s, double t, Vec_t & N)  const;
	void   Derivs     (double r, double s, double t, Mat_t & dN) const;
	void   FaceShape  (double r, double s, Vec_t & FN)  const;
	void   FaceDerivs (double r, double s, Mat_t & FdN) const;
private:
	void _local_coords (Mat_t & coords) const;

}; // class Hex8

/* Local IDs
                  Vertices                             Edges                              Faces
    z
    |           4                  7
   ,+--y         @________________@                 +_______(4)______+                 +________________+ 
 x'            ,'|              ,'|               ,'|              ,'|               ,'|              ,'| 
             ,'  |            ,'  |             ,'  |            ,'  |             ,'  |  ___       ,'  | 
           ,'    |          ,'    |           (6)  (8)         (7)   |           ,'    |,'5,'  [0],'    | 
     5   ,'      |      6 ,'      |         ,'      |        ,'     (11)       ,'      |~~~     ,'      | 
       @'===============@'        |       +'==========(5)==+'        |       +'===============+'  ,'|   | 
       |         |      |         |       |         |      |         |       |   ,'|   |      |   |3|   | 
       |         |      |         |       |         |      |         |       |   |2|   |      |   |,'   | 
       |       0 @______|_________@       |         +______|_(0)_____+       |   |,'   +______|_________+ 
       |       ,'       |       ,' 3     (9)      ,'       |       ,'        |       ,'       |       ,'  
       |     ,'         |     ,'          |    (2)        (10)   ,'          |     ,' [1]  ___|     ,'    
       |   ,'           |   ,'            |   ,'           |   (3)           |   ,'      ,'4,'|   ,'      
       | ,'             | ,'              | ,'             | ,'              | ,'        ~~~  | ,'        
       @________________@'                +______(1)_______+'                +________________+'          
     1                   2
*/
Hex8::FaceMap Hex8::Face2Node[]= {{ 0, 3, 7, 4 },
                                  { 1, 2, 6, 5 },
                                  { 0, 1, 5, 4 },
                                  { 2, 3, 7, 6 }, 
                                  { 0, 1, 2, 3 }, 
                                  { 4, 5, 6, 7 }};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Hex8::Hex8()
{
	// Setup nodes number
	NNodes  = 8;
	NFNodes = 4;

	// Allocate nodes (connectivity)
	Conn.Resize    (NNodes);
	Conn.SetValues (NULL);

	// Integration points and Extrapolation Matrix
	SetIPs (/*NIPs1D*/2);
}

inline void Hex8::SetIPs(int NIPs1D)
{
	// Setup pointer to the array of Integration Points
	     if (NIPs1D==2) IPs = HEX_IP2;
	else if (NIPs1D==3) IPs = HEX_IP3;
	else if (NIPs1D==4) IPs = HEX_IP4;
	else if (NIPs1D==5) IPs = HEX_IP5;
	else throw new Fatal("Hex8::SetIPs: Number of integration points < %d > is invalid",NIPs1D);

	NIPs  = pow(NIPs1D, 3);
	FIPs  = QUAD_IP2;
	NFIPs = 4; 
}


inline void Hex8::VTKConnect(String & Nodes) const
{
	Nodes.Printf("%d %d %d %d %d %d %d %d",Conn[1]->GetID(),
	                                       Conn[2]->GetID(),
	                                       Conn[3]->GetID(),
	                                       Conn[0]->GetID(),
	                                       Conn[5]->GetID(),
	                                       Conn[6]->GetID(),
	                                       Conn[7]->GetID(),
	                                       Conn[4]->GetID());
}

inline void Hex8::GetFaceNodes(int FaceID, Array<Node*> & FaceConnects) const
{
	FaceConnects.Resize(/*NumFaceNodes*/4);
	FaceConnects[0] = Conn[Face2Node[FaceID].n0];
	FaceConnects[1] = Conn[Face2Node[FaceID].n1];
	FaceConnects[2] = Conn[Face2Node[FaceID].n2];
	FaceConnects[3] = Conn[Face2Node[FaceID].n3];
}

inline void Hex8::Shape(double r, double s, double t, Vec_t & N) const
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
	N.Resize(/*NumNodes*/8);
	N(0) = 0.125*(1.0-r-s+r*s-t+s*t+r*t-r*s*t);
	N(1) = 0.125*(1.0+r-s-r*s-t+s*t-r*t+r*s*t);
	N(2) = 0.125*(1.0+r+s+r*s-t-s*t-r*t-r*s*t);
	N(3) = 0.125*(1.0-r+s-r*s-t-s*t+r*t+r*s*t);
	N(4) = 0.125*(1.0-r-s+r*s+t-s*t-r*t+r*s*t);
	N(5) = 0.125*(1.0+r-s-r*s+t-s*t+r*t-r*s*t);
	N(6) = 0.125*(1.0+r+s+r*s+t+s*t+r*t+r*s*t);
	N(7) = 0.125*(1.0-r+s-r*s+t+s*t-r*t-r*s*t);
}

inline void Hex8::Derivs(double r, double s, double t, Mat_t & dN) const
{
	/*       _     _ T
	 *      |  dNi  |
	 * dN = |  ---  |   , where cj = r, s
	 *      |_ dcj _|
	 *
	 * dN(j,i), j=>local coordinate and i=>shape function
	 */
	dN.Resize(3,/*NumNodes*/8);
	dN(0,0) = 0.125*(-1.0+s+t-s*t);   dN(1,0)=0.125*(-1.0+r+t-r*t);   dN(2,0)=0.125*(-1.0+r+s-r*s);
	dN(0,1) = 0.125*(+1.0-s-t+s*t);   dN(1,1)=0.125*(-1.0-r+t+r*t);   dN(2,1)=0.125*(-1.0-r+s+r*s);
	dN(0,2) = 0.125*(+1.0+s-t-s*t);   dN(1,2)=0.125*(+1.0+r-t-r*t);   dN(2,2)=0.125*(-1.0-r-s-r*s);
	dN(0,3) = 0.125*(-1.0-s+t+s*t);   dN(1,3)=0.125*(+1.0-r-t+r*t);   dN(2,3)=0.125*(-1.0+r-s+r*s);
	dN(0,4) = 0.125*(-1.0+s-t+s*t);   dN(1,4)=0.125*(-1.0+r-t+r*t);   dN(2,4)=0.125*(+1.0-r-s+r*s);
	dN(0,5) = 0.125*(+1.0-s+t-s*t);   dN(1,5)=0.125*(-1.0-r-t-r*t);   dN(2,5)=0.125*(+1.0+r-s-r*s);
	dN(0,6) = 0.125*(+1.0+s+t+s*t);   dN(1,6)=0.125*(+1.0+r+t+r*t);   dN(2,6)=0.125*(+1.0+r+s+r*s);
	dN(0,7) = 0.125*(-1.0-s-t-s*t);   dN(1,7)=0.125*(+1.0-r+t-r*t);   dN(2,7)=0.125*(+1.0-r+s-r*s);
}

inline void Hex8::FaceShape(double r, double s, Vec_t & FN) const
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

	FN.Resize(/*NumFaceNodes*/4);
	FN(0) = 0.25*(1.0-r-s+r*s);
	FN(1) = 0.25*(1.0+r-s-r*s);
	FN(2) = 0.25*(1.0+r+s+r*s);
	FN(3) = 0.25*(1.0-r+s-r*s);
}

inline void Hex8::FaceDerivs(double r, double s, Mat_t & FdN) const
{
	/*          _     _ T
	 *         |  dNi  |
	 *   FdN = |  ---  |   , where cj = r, s
	 *         |_ dcj _|
	 *
	 *   FdN(j,i), j=>local coordinate and i=>shape function
	 */

	FdN.Resize(2,/*NumFaceNodes*/4);
	FdN(0,0) = 0.25*(-1.0+s);   FdN(1,0) = 0.25*(-1.0+r);
	FdN(0,1) = 0.25*(+1.0-s);   FdN(1,1) = 0.25*(-1.0-r);
	FdN(0,2) = 0.25*(+1.0+s);   FdN(1,2) = 0.25*(+1.0+r);
	FdN(0,3) = 0.25*(-1.0-s);   FdN(1,3) = 0.25*(+1.0-r);
}

inline void Hex8::_local_coords(Mat_t & coords) const 
{
	coords.Resize(8,4);
	coords = -1.0, -1.0, -1.0, 1.0,
	         +1.0, -1.0, -1.0, 1.0,
	         +1.0, +1.0, -1.0, 1.0,
	         -1.0, +1.0, -1.0, 1.0,
	         -1.0, -1.0, +1.0, 1.0,
	         +1.0, -1.0, +1.0, 1.0,
	         +1.0, +1.0, +1.0, 1.0,
	         -1.0, +1.0, +1.0, 1.0;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
GeomElem * Hex8Maker() { return new Hex8(); }

// Register element
int Hex8Register() { GeomElemFactory["Hex8"]=Hex8Maker; return 0; }

// Call register
int __Hex8_dummy_int  = Hex8Register();

}; // namespace FEM

#endif // MECHSYS_FEM_HEX8_H
