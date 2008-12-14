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

#ifndef MECHSYS_FEM_LIN3_H
#define MECHSYS_FEM_LIN3_H

// MechSys
#include "fem/node.h"
#include "fem/geomelem.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "util/exception.h"
#include "fem/elems/vtkCellType.h"

namespace FEM
{

class Lin3 : public GeomElem
{
public:
	// Constructor
	Lin3 ();

	// Derived methods
	void   SetIPs     (int NIPs1D);
	int    VTKType    () const { return VTK_QUADRATIC_EDGE; }
	void   VTKConn    (String & Nodes) const { Nodes.Printf("%d %d %d",Conn[0]->GetID(),Conn[2]->GetID(),Conn[1]->GetID()); }
	void   Shape      (double r, double s, double t, Vec_t & N)  const;
	void   Derivs     (double r, double s, double t, Mat_t & dN) const;

private:
	void _local_coords (Mat_t & coords) const;

}; // class Lin3


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Lin3::Lin3()
{
	// Setup nodes number
	NNodes  = 3;
	NFNodes = 0;

	// Allocate nodes (connectivity)
	Conn.Resize    (NNodes);
	Conn.SetValues (NULL);

	// Integration points
	SetIPs (/*NIPs1D*/2);
}

inline void Lin3::SetIPs(int NIPs1D)
{
	// Setup pointer to the array of Integration Points
	     if (NIPs1D==2) IPs = LIN_IP2;
	else if (NIPs1D==3) IPs = LIN_IP3;
	else throw new Fatal("Lin3::SetIntPoints: Error in number of integration points.");

	NIPs  = NIPs1D;
	NFIPs = 0;
}

inline void Lin3::Shape(double r, double s, double t, Vec_t & N) const 
{
	/*
	 *       -----o=========o=========o----->  r
	 *            0         1         2
	 */

	N.Resize(3);
	N(0) = 0.5*(r*r-r);
	N(1) = 1.0 - r*r;
	N(2) = 0.5*(r*r+r);
}

inline void Lin3::Derivs(double r, double s, double t, Mat_t & dN) const 
{
	dN.Resize(1,3);
	dN(0,0) =  r-0.5;
	dN(0,1) = -2.0*r;
	dN(0,2) =  r+0.5;
}

inline void Lin3::LocalCoords(Mat_t & coords) const 
{
	if (_ndim==2)
	{
		coords.Resize(3,3);
		coords =  -1.0,  0.0, 1.0,
				   0.0,  0.0, 1.0,
				   1.0,  0.0, 1.0;
	}
	else
	{
		coords.Resize(3,4);
		coords =  -1.0,  0.0, 0.0, 1.0,
				   0.0,  0.0, 0.0, 1.0,
				   1.0,  0.0, 0.0, 1.0;
	}
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate new element
Element * Lin3Maker() { return new Lin3(); }

// Register new element
int Lin3Register() { ElementFactory["Lin3"]=Lin3Maker;  return 0; }

// Call register
int __Lin3_dummy_int  = Lin3Register();


}; // namespace FEM

#endif // MECHSYS_FEM_LIN3_H
