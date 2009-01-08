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

#ifndef MECHSYS_FEM_QUAD8_H
#define MECHSYS_FEM_QUAD8_H

// MechSys
#include "fem/node.h"
#include "fem/geomelem.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "vtkCellType.h"

namespace FEM
{

class Quad8: public GeomElem
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
	Quad8();

	// Derived methods
	void   SetIPs     (int NIPs1D);
	int    VTKType    () const { return VTK_QUADRATIC_QUAD; }
	void   VTKConn    (String & Nodes) const;
	void   GetFNodes  (int FaceID, Array<Node*> & FaceConnects) const;
	double BoundDist  (double r, double s, double t) const { return std::min(1-fabs(r),1-fabs(s)); }
	void   Shape      (double r, double s, double t, Vec_t & N)  const;
	void   Derivs     (double r, double s, double t, Mat_t & dN) const;
	void   FaceShape  (double r, double s, Vec_t & FN)  const;
	void   FaceDerivs (double r, double s, Mat_t & FdN) const;

private:
	void _local_coords (Mat_t & coords) const;


}; // class Quad8

/* Local IDs
                 Nodes                  Faces(edges)
   y
   |        3      6      2                  y+
   +--x      @-----@-----@             +----(3)----+
             |           |             |           |
             |           |             |           |
           7 @           @ 5       x- (0)         (1) x+
             |           |             |           |
             |           |             |           |
             @-----@-----@             +----(2)----+
            0      4      1                  y-
*/
Quad8::FaceMap Quad8::Face2Node[]= {{ 3, 0, 7 },
                                    { 1, 2, 5 },
                                    { 0, 1, 4 },
                                    { 2, 3, 6 }}; // order of nodes is important


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Quad8::Quad8()
{
	// Setup nodes number
	NNodes  = 8;
	NFNodes = 3;

	// Allocate nodes (connectivity)
	Conn.Resize    (NNodes);
	Conn.SetValues (NULL);

	// Integration points and Extrapolation Matrix
	SetIPs (/*NIPs1D*/2);
}

inline void Quad8::SetIPs(int NIPs1D)
{
	// Setup pointer to the array of Integration Points
	     if (NIPs1D==2) IPs = QUAD_IP2;
	else if (NIPs1D==3) IPs = QUAD_IP3;
	else if (NIPs1D==4) IPs = QUAD_IP4;
	else if (NIPs1D==5) IPs = QUAD_IP5;
	else throw new Fatal("Quad8::SetIPs: Number of integration points < %d > is invalid",NIPs1D);

	NIPs  = pow(NIPs1D, 2);
	FIPs  = LIN_IP2;
	NFIPs = 2;
}


inline void Quad8::VTKConn(String & Nodes) const
{
	Nodes.Printf("%d %d %d %d %d %d %d %d",Conn[0]->ID(),
	                                       Conn[1]->ID(),
	                                       Conn[2]->ID(),
	                                       Conn[3]->ID(),
	                                       Conn[4]->ID(),
	                                       Conn[5]->ID(),
	                                       Conn[6]->ID(),
	                                       Conn[7]->ID());
}

inline void Quad8::GetFNodes(int FaceID, Array<Node*> & FaceConnects) const
{
	FaceConnects.Resize(3);
	FaceConnects[0] = Conn[Face2Node[FaceID].L];
	FaceConnects[1] = Conn[Face2Node[FaceID].R];
	FaceConnects[2] = Conn[Face2Node[FaceID].M];
}

inline void Quad8::Shape(double r, double s, double t, Vec_t & N) const
{
	/*      3           6            2
	 *        @---------@----------@
	 *        |               (1,1)|
	 *        |       s ^          |
	 *        |         |          |
	 *        |         |          |
	 *      7 @         +----> r   @ 5
	 *        |       (0,0)        |
	 *        |                    |
	 *        |                    |
	 *        |(-1,-1)             |
	 *        @---------@----------@
	 *      0           4            1
	 */
	N.Resize (/*NumNodes*/8);

	double rp1=1.0+r; double rm1=1.0-r;
	double sp1=1.0+s; double sm1=1.0-s;

	N(0) = 0.25*rm1*sm1*(rm1+sm1-3.0);
	N(1) = 0.25*rp1*sm1*(rp1+sm1-3.0);
	N(2) = 0.25*rp1*sp1*(rp1+sp1-3.0);
	N(3) = 0.25*rm1*sp1*(rm1+sp1-3.0);
	N(4) = 0.50*sm1*(1.0-r*r);
	N(5) = 0.50*rp1*(1.0-s*s);
	N(6) = 0.50*sp1*(1.0-r*r);
	N(7) = 0.50*rm1*(1.0-s*s);
}

inline void Quad8::Derivs(double r, double s, double t, Mat_t & dN) const
{
	/*           _     _ T
	 *          |  dNi  |
	 * dN = |  ---  |   , where cj = r, s
	 *          |_ dcj _|
	 *
	 * dN(j,i), j=>local coordinate and i=>shape function
	 */
	dN.Resize (2, /*NumNodes*/8);

	double rp1=1.0+r; double rm1=1.0-r;
	double sp1=1.0+s; double sm1=1.0-s;

	dN(0,0) = - 0.25 * sm1 * (rm1 + rm1 + sm1 - 3.0);
	dN(0,1) =   0.25 * sm1 * (rp1 + rp1 + sm1 - 3.0);
	dN(0,2) =   0.25 * sp1 * (rp1 + rp1 + sp1 - 3.0);
	dN(0,3) = - 0.25 * sp1 * (rm1 + rm1 + sp1 - 3.0);
	dN(0,4) = - r * sm1;
	dN(0,5) =   0.50 * (1.0 - s * s);
	dN(0,6) = - r * sp1;
	dN(0,7) = - 0.5 * (1.0 - s * s);

	dN(1,0) = - 0.25 * rm1 * (sm1 + rm1 + sm1 - 3.0);
	dN(1,1) = - 0.25 * rp1 * (sm1 + rp1 + sm1 - 3.0);
	dN(1,2) =   0.25 * rp1 * (sp1 + rp1 + sp1 - 3.0);
	dN(1,3) =   0.25 * rm1 * (sp1 + rm1 + sp1 - 3.0);
	dN(1,4) = - 0.50 * (1.0 - r * r);
	dN(1,5) = - s * rp1;
	dN(1,6) =   0.50 * (1.0 - r * r);
	dN(1,7) = - s * rm1;
}

inline void Quad8::FaceShape(double r, double s, Vec_t & FN) const
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

inline void Quad8::FaceDerivs(double r, double s, Mat_t & FdN) const
{
	/*          _     _ T
	 *         |  dNi  |
	 *   FdN = |  ---  |   , where cj = r, s
	 *         |_ dcj _|
	 *
	 *   FdN(j,i), j=>local coordinate and i=>shape function
	 */
	FdN.Resize(1,/*NumFaceNodes*/3);
	FdN(0,0) =  r  - 0.5;
	FdN(0,1) =  r  + 0.5;
	FdN(0,2) = -2.0* r;
}

inline void Quad8::_local_coords(Mat_t & C) const 
{
	C.Resize(8,3);
	C = -1.0, -1.0, 1.0,   
	    +1.0, -1.0, 1.0,   
	    +1.0, +1.0, 1.0,   
	    -1.0, +1.0, 1.0,   
	     0.0, -1.0, 1.0,   
	    +1.0,  0.0, 1.0,   
	     0.0, +1.0, 1.0,   
	    -1.0, +0.0, 1.0;   
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
GeomElem * Quad8Maker() { return new Quad8(); }

// Register element
int Quad8Register() { GeomElemFactory["Quad8"]=Quad8Maker;  return 0; }

// Call register
int __Quad8_dummy_int  = Quad8Register();
}; // namespace FEM

#endif // MECHSYS_FEM_QUAD8_H
