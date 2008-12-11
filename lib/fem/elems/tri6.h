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

#ifndef MECHSYS_FEM_TRI6_H
#define MECHSYS_FEM_TRI6_H

// MechSys
#include "fem/node.h"
#include "fem/geomelem.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "vtkCellType.h"

namespace FEM
{

class Tri6: public GeomElem
{
public:
	// Auxiliar structure to map local face IDs to local node IDs
	struct FaceMap
	{
		int L; // Left node local id
		int R; // Right node local id
		int M; // Mid node local id
	};
	static FaceMap Face2Node[];

	// Constructor
	Tri6();

	// Derived methods
	void   SetIPs     (int NIPs1D);
	int    VTKType    () const { return VTK_QUADRATIC_TRIANGLE; }
	void   VTKConn    (String & Nodes) const;
	void   GetFNodes  (int FaceID, Array<Node*> & FaceConnects) const;
	double BoundDist  (double r, double s, double t) const { return std::min(std::min(r,s), 1-r-s); }
	void   Shape      (double r, double s, double t, Vec_t & N)  const;
	void   Derivs     (double r, double s, double t, Mat_t & dN) const;
	void   FaceShape  (double r, double s, Vec_t & FN)  const;
	void   FaceDerivs (double r, double s, Mat_t & FdN) const;

private:
	void _local_coords (Mat_t & coords) const;

}; // class Tri6


/* Local IDs
             Nodes                 Faces

   y           2
   |           @                     @
   +--x       / \                   / \
           5 /   \ 4               /   \
            @     @             2 /     \ 1
           /       \             /       \
          /         \           /         \
         @-----@-----@         @-----------@
        0      3      1              0
*/
Tri6::FaceMap Tri6::Face2Node[]= {{ 0, 1, 3 },
                                  { 1, 2, 4 },
                                  { 2, 0, 5 }}; // order of nodes is important


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Tri6::Tri6()
{
	// Setup nodes number
	NNodes  = 6;
	NFNodes = 3;

	// Allocate nodes (connectivity)
	Conn.Resize    (NNodes);
	Conn.SetValues (NULL);

	// Integration points and Extrapolation Matrix
	SetIPs (/*NIPsTotal*/3);
}

inline void Tri6::SetIPs(int NIPsTotal)
{
	// Setup pointer to the array of Integration Points
	if      (NIPsTotal==3)  IPs = TRI_IP3;
	else if (NIPsTotal==4)  IPs = TRI_IP4;
	else if (NIPsTotal==6)  IPs = TRI_IP6;
	else if (NIPsTotal==7)  IPs = TRI_IP7;
	else if (NIPsTotal==13) IPs = TRI_IP13;
	else throw new Fatal("tri6::SetIntPoints: Error in number of integration points.");

	NIPs  = NIPsTotal ;
	FIPs  = LIN_IP2;
	NFIPs = 2;
}

inline void Tri6::VTKConnect(String & Nodes) const
{
	Nodes.Printf("%d %d %d %d %d %d",Conn[0]->GetID(),
	                                 Conn[1]->GetID(),
	                                 Conn[2]->GetID(),
	                                 Conn[3]->GetID(),
	                                 Conn[4]->GetID(),
	                                 Conn[5]->GetID());
}

inline void Tri6::GetFaceNodes(int FaceID, Array<Node*> & FaceConnects) const
{
	FaceConnects.Resize(3);
	FaceConnects[0] = Conn[Face2Node[FaceID].L];
	FaceConnects[1] = Conn[Face2Node[FaceID].R];
	FaceConnects[2] = Conn[Face2Node[FaceID].M];
}

inline void Tri6::Shape(double r, double s, double t, Vec_t & N) const
{

	/*    s
	 *    ^
	 *    |
	 *  2
	 *    @,(0,1)
	 *    | ',
	 *    |   ',
	 *    |     ',
	 *    |       ',   4
	 *  5 @          @
	 *    |           ',
	 *    |             ',
	 *    |               ',
	 *    |(0,0)            ', (1,0)
	 *    @---------@---------@  --> r
	 *  0           3          1
	 */
	N.Resize(/*NumNodes*/6);
	N(0) = 1.0-(r+s)*(3.0-2.0*(r+s));
	N(1) = r*(2.0*r-1.0);
	N(2) = s*(2.0*s-1.0);
	N(3) = 4.0*r*(1.0-(r+s));
	N(4) = 4.0*r*s;
	N(5) = 4.0*s*(1.0-(r+s));
}

inline void Tri6::Derivs(double r, double s, double t, Mat_t & dN) const
{
	/*         _     _ T
	 *        |  dNi  |
	 *   dN = |  ---  |   , where cj = r, s
	 *        |_ dcj _|
	 *  
	 *   dN(j,i), j=>local coordinate and i=>shape function
	 */
	dN.Resize(2, /*NumNodes*/6);

	dN(0,0) = -3.0 + 4.0 * (r + s);       dN(1,0) = -3.0 + 4.0*(r + s);
	dN(0,1) =  4.0 * r - 1.;              dN(1,1) =  0.0 ;
	dN(0,2) =  0.0;                       dN(1,2) =  4.0 * s - 1.0;
	dN(0,3) =  4.0 - 8.0 * r - 4.0 * s;   dN(1,3) = -4.0 * r;
	dN(0,4) =  4.0 * s;                   dN(1,4) =  4.0 * r;
	dN(0,5) = -4.0 * s;                   dN(1,5) =  4.0 - 4.0 * r - 8.0*s;
}

inline void Tri6::FaceShape(double r, double s, Vec_t & FN) const
{
	/*
	 *       @-----------@-----------@-> r
	 *       0           2           1
	 *       |           |           |
	 *      r=-1         r=0        r=+1
	 */
	FN.Resize(/*NumFaceNodes*/3);
	FN(0) = 0.5 * (r*r-r);
	FN(1) = 0.5 * (r*r+r);
	FN(2) = 1.0 -  r*r;
}

inline void Tri6::FaceDerivs(double r, double s, Mat_t & FdN) const
{
	FdN.Resize(1,/*NumFaceNodes*/3);
	FdN(0,0) =  r  - 0.5;
	FdN(0,1) =  r  + 0.5;
	FdN(0,2) = -2.0* r;
}

inline void Tri6::_local_coords(Mat_t & C) const
{
	C.Resize(6,3);
	C =  0.0,  0.0, 1.0,
	     1.0,  0.0, 1.0,
	     0.0,  1.0, 1.0,
	     0.5,  0.0, 1.0,
	     0.5,  0.5, 1.0,
	     0.0,  0.5, 1.0;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
GeomElem * Tri6Maker() { return new Tri6(); }

// Register element
int Tri6Register() { GeomElemFactory["Tri6"]=Tri6Maker;  return 0; }

// Call register
int __Tri6_dummy_int  = Tri6Register();

}; // namespace FEM

#endif // MECHSYS_FEM_TRI6_H
