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
#include "fem/geomelem.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "util/exception.h"
#include "fem/elems/vtkCellType.h"

namespace FEM
{

class Quad4: public GeomElem
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
	Quad4 ();

	// Derived methods
	void   SetIPs     (int NIPs1D);
	int    VTKType    () const { return VTK_QUAD; }
	void   VTKConn    (String & Nodes) const;
	void   GetFNodes  (int FaceID, Array<Node*> & FaceConnects) const;
	double BoundDist  (double r, double s, double t) const { return std::min(1-fabs(r),1-fabs(s)); }
	void   Shape      (double r, double s, double t, Vec_t & N)  const;
	void   Derivs     (double r, double s, double t, Mat_t & dN) const;
	void   FaceShape  (double r, double s, Vec_t & FN)  const;
	void   FaceDerivs (double r, double s, Mat_t & FdN) const;

private:
	void _local_coords (Mat_t & coords) const;

}; // class Quad4


/* Local IDs
                 Nodes                  Faces(edges)
   y
   |        3             2                  y+
   +--x      @-----------@             +----(3)----+
             |           |             |           |
             |           |             |           |
             |           |         x- (0)         (1) x+
             |           |             |           |
             |           |             |           |
             @-----------@             +----(2)----+
            0             1                  y-
*/
Quad4::FaceMap Quad4::Face2Node[]= {{ 3, 0 },
                                    { 1, 2 },
                                    { 0, 1 },
                                    { 2, 3 }}; // order of nodes is important


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Quad4::Quad4()
{
	// Setup nodes number
	NNodes  = 4;
	NFNodes = 2;

	// Allocate nodes (connectivity)
	Conn.Resize    (NNodes);
	Conn.SetValues (NULL);

	// Integration points
	SetIPs (/*NIPs1D*/2);
}

inline void Quad4::SetIPs(int NIPs1D)
{
	// Setup pointer to the array of Integration Points
	     if (NIPs1D==2) IPs = QUAD_IP2;
	else if (NIPs1D==3) IPs = QUAD_IP3;
	else if (NIPs1D==4) IPs = QUAD_IP4;
	else if (NIPs1D==5) IPs = QUAD_IP5;
	else throw new Fatal("Quad4::SetIPs: Number of integration points < %d > is invalid",NIPs1D);

	NIPs  = pow(NIPs1D, 2);
	FIPs  = LIN_IP2;
	NFIPs = 2;
}

inline void Quad4::VTKConn(String & Nodes) const
{
	Nodes.Printf("%d %d %d %d",Conn[0]->GetID(),
	                           Conn[1]->GetID(),
	                           Conn[2]->GetID(),
	                           Conn[3]->GetID());
}

inline void Quad4::GetFNodes(int FaceID, Array<Node*> & FaceConnects) const
{
	FaceConnects.Resize(2);
	FaceConnects[0] = Conn[Face2Node[FaceID].L];
	FaceConnects[1] = Conn[Face2Node[FaceID].R];
}

inline void Quad4::Shape(double r, double s, double t, Vec_t & N) const
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
	N.Resize (/*NNodes*/4);

	double rp=1.0+r; double rm=1.0-r;
	double sp=1.0+s; double sm=1.0-s;

	N(0) = 0.25*rm*sm;
	N(1) = 0.25*rp*sm;
	N(2) = 0.25*rp*sp;
	N(3) = 0.25*rm*sp;
}

inline void Quad4::Derivs(double r, double s, double t, Mat_t & dN) const
{
	/*         _     _ T
	 *        |  dNi  |
	 *   dN = |  ---  |   , where cj = r, s
	 *        |_ dcj _|
	 *  
	 *   dN(j,i), j=>local coordinate and i=>shape function
	 */
	dN.Resize (2, /*NNodes*/4);

	double rp=1.0+r; double rm=1.0-r;
	double sp=1.0+s; double sm=1.0-s;

	dN(0,0) = -0.25*sm;   dN(1,0) = -0.25*rm;
	dN(0,1) =  0.25*sm;   dN(1,1) = -0.25*rp;
	dN(0,2) =  0.25*sp;   dN(1,2) =  0.25*rp;
	dN(0,3) = -0.25*sp;   dN(1,3) =  0.25*rm;
}

inline void Quad4::FaceShape(double r, double s, Vec_t & FN) const
{
	/*  
	 *       0           |           1
	 *       @-----------+-----------@-> r
	 *      -1           |          +1
	 */
	FN.Resize(/*NumFNodes*/2);
	FN(0) = 0.5*(1.0-r);
	FN(1) = 0.5*(1.0+r);
}

inline void Quad4::FaceDerivs(double r, double s, Mat_t & FdN) const
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

inline void Quad4::_local_coords(Mat_t & C) const 
{
	C.Resize(4,3);
	C = -1.0, -1.0, 1.0,
	     1.0, -1.0, 1.0,
	     1.0,  1.0, 1.0,
	    -1.0,  1.0, 1.0;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
GeomElem * Quad4Maker() { return new Quad4(); }

// Register element
int Quad4Register() { GeomElemFactory["Quad4"]=Quad4Maker;  return 0; }

// Call register
int __Quad4_dummy_int  = Quad4Register();

}; // namespace FEM

#endif // MECHSYS_FEM_QUAD4_H
