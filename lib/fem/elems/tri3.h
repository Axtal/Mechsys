/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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

#ifndef MECHSYS_FEM_TRI3_H
#define MECHSYS_FEM_TRI3_H

// MechSys
#include "fem/node.h"
#include "fem/geomelem.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "vtkCellType.h"

namespace FEM
{

class Tri3: public GeomElem
{
public:
	// Auxiliar structure to map local face IDs to local node IDs
	struct FaceMap
	{
		int L; // Left node local id
		int R; // Right node local id
	};
	static FaceMap Face2Node[];

	// Constructor
	Tri3();

	// Derived methods
	void   SetIPs     (int NIPs1D);
	int    VTKType    () const { return VTK_TRIANGLE; }
	void   VTKConn    (String & Nodes) const;
	void   GetFNodes  (int FaceID, Array<Node*> & FaceConnects) const;
	double BoundDist  (double r, double s, double t) const { return std::min(std::min(r,s), 1-r-s); }
	void   Shape      (double r, double s, double t, Vec_t & N)  const;
	void   Derivs     (double r, double s, double t, Mat_t & dN) const;
	void   FaceShape  (double r, double s, Vec_t & FN)  const;
	void   FaceDerivs (double r, double s, Mat_t & FdN) const;

private:
	void _local_coords (Mat_t & coords) const;

}; // class Tri3


/* Local IDs
             Nodes                 Faces 
                           
   y           2                                                        
   |           @                     @                                  
   +--x       / \                   / \                                  
             /   \                 /   \                                 
            /     \             2 /     \ 1                              
           /       \             /       \                               
          /         \           /         \                              
         @-----------@         @-----------@                             
        0             1              0               
*/
Tri3::FaceMap Tri3::Face2Node[]= {{ 0, 1 },
                                  { 1, 2 },
                                  { 2, 0 }}; // order of nodes is important


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Tri3::Tri3()
{
	// Setup nodes number
	NNodes  = 3;
	NFNodes = 2;

	// Allocate nodes (connectivity)
	Conn.Resize    (NNodes);
	Conn.SetValues (NULL);

	// Integration points and Extrapolation Matrix
	SetIPs (/*NIPsTotal*/3);
}

inline void Tri3::SetIPs(int NIPsTotal)
{
	// Setup pointer to the array of Integration Points
	if      (NIPsTotal==3)  IPs = TRI_IP3;
	else if (NIPsTotal==4)  IPs = TRI_IP4;
	else if (NIPsTotal==6)  IPs = TRI_IP6;
	else if (NIPsTotal==7)  IPs = TRI_IP7;
	else if (NIPsTotal==13) IPs = TRI_IP13;
	else throw new Fatal("tri3::SetIntPoints: Error in number of integration points.");

	NIPs  = NIPsTotal ;
	FIPs  = LIN_IP2;
	NFIPs = 2;
}

inline void Tri3::VTKConn(String & Nodes) const
{
	Nodes.Printf("%d %d %d",Conn[0]->ID(),
	                        Conn[1]->ID(),
	                        Conn[2]->ID());
}

inline void Tri3::GetFNodes(int FaceID, Array<Node*> & FaceConnects) const
{
	FaceConnects.Resize(2);
	FaceConnects[0] = Conn[Face2Node[FaceID].L];
	FaceConnects[1] = Conn[Face2Node[FaceID].R];
}

inline void Tri3::Shape(double r, double s, double t, Vec_t & N) const
{

	/*    s
	 *    ^
	 *    |
	 *  2 
	 *    @,(0,1)
	 *    | ',
	 *    |   ',
	 *    |     ',
	 *    |       ',    
	 *    |         ',
	 *    |           ', 
	 *    |             ', 
	 *    |               ', 
	 *    |(0,0)            ', (1,0)
	 *    @-------------------@  --> r
	 *  0                      1
	 */
	N.Resize(/*NumNodes*/3);
	N(0) = 1.0-r-s;
	N(1) = r;
	N(2) = s;
}

inline void Tri3::Derivs(double r, double s, double t, Mat_t & dN) const
{
	/*         _     _ T
	 *        |  dNi  |
	 *   dN = |  ---  |   , where cj = r, s
	 *        |_ dcj _|
	 *  
	 *   dN(j,i), j=>local coordinate and i=>shape function
	 */
	dN.Resize(2, /*NumNodes*/3);
	dN(0,0) = -1.0;    dN(1,0) = -1.0;
	dN(0,1) =  1.0;    dN(1,1) =  0.0;
	dN(0,2) =  0.0;    dN(1,2) =  1.0;
}

inline void Tri3::FaceShape(double r, double s, Vec_t & FN) const
{
	/*  
	 *       0           |           1
	 *       @-----------+-----------@-> r
	 *      -1           |          +1
	 */
	FN.Resize(/*NumFaceNodes*/2);
	FN(0) = 0.5*(1.0-r);
	FN(1) = 0.5*(1.0+r);
}

inline void Tri3::FaceDerivs(double r, double s, Mat_t & FdN) const
{
	/*          _     _ T
	 *         |  dNi  |
	 *   FdN = |  ---  |   , where cj = r, s
	 *         |_ dcj _|
	 *
	 *   FdN(j,i), j=>local coordinate and i=>shape function
	 */
	FdN.Resize(1,/*NumFNodes*/2);
	FdN(0,0) = -0.5;
	FdN(0,1) =  0.5;
}

inline void Tri3::_local_coords(Mat_t & C) const 
{
	C.Resize(3,3);
	C = 0.0, 0.0, 1.0,
	    1.0, 0.0, 1.0,
	    0.0, 1.0, 1.0;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
GeomElem * Tri3Maker() { return new Tri3(); }

// Register element
int Tri3Register() { GeomElemFactory["Tri3"]=Tri3Maker;  return 0; }

// Call register
int __Tri3_dummy_int  = Tri3Register();

}; // namespace FEM

#endif // MECHSYS_FEM_TRI3_H
