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

#ifndef MECHSYS_FEM_LIN2_H
#define MECHSYS_FEM_LIN2_H

// MechSys
#include "fem/node.h"
#include "fem/geomelem.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "vtkCellType.h"

namespace FEM
{

class Lin2: public GeomElem
{
public:
	// Constructor
	Lin2 ();

	// Derived methods
	void   SetIPs     (int NIPs1D);
	int    VTKType    () const { return VTK_QUAD; }
	void   VTKConn    (String & Nodes) const { Nodes.Printf("%d %d %d",Conn[0]->GetID(),Conn[2]->GetID(),Conn[1]->GetID()); }
	void   GetFNodes  (int FaceID, Array<Node*> & FaceConnects)  const;
	void   Shape      (double r, double s, double t, Vec_t & N)  const;
	void   Derivs     (double r, double s, double t, Mat_t & dN) const;
	void   FaceShape  (double r, double s, Vec_t & FN)           const;
	void   FaceDerivs (double r, double s, Mat_t & FdN)          const;

private:
	void _local_coords (Mat_t & coords) const;

}; // class Rod3


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Lin2::Lin2()
{
	// Setup nodes number
	NNodes  = 2;
	NFNodes = 0;

	// Allocate nodes (connectivity)
	Conn.Resize    (NNodes);
	Conn.SetValues (NULL);

	// Integration points 
	SetIPs (/*NIPs1D*/2);
	NIPs  = 0;
}

inline void Lin2::SetIPs(int NIPs1D)
{
	throw new Fatal("Lin2::SetIPs: This feature was not implemented yet.");
}

inline void Lin2::Shape(double r, double s, double t, Vec_t & N) const 
{
	throw new Fatal("Lin2::Shape: This feature was not implemented yet.");
}

inline void Lin2::Derivs(double r, double s, double t, Mat_t & dN) const 
{
	throw new Fatal("Lin2::Derivs: This feature was not implemented yet.");
}

inline void Lin2::FaceShape  (double r, double s, Vec_t & FN) const
{
	throw new Fatal("Lin2::FaceShape: This feature was not availible for this element.");
}

inline void Lin2::FaceDerivs (double r, double s, Mat_t & FdN) const
{
	throw new Fatal("Lin2::FaceDerivs: This feature was not availible for this element.");
}

inline void Lin2::GetFNodes (int FaceID, Array<Node*> & FaceConnects) const
{
	throw new Fatal("Lin2::GetFNodes: This feature was not availible for this element.");

}

inline void Lin2::_local_coords(Mat_t & coords) const 
{
	if (NDim==2)
	{
		coords.Resize(2,3);
		coords =  -1.0,  0.0, 1.0,
				   1.0,  0.0, 1.0;
	}
	else
	{
		coords.Resize(2,4);
		coords =  -1.0,  0.0, 0.0, 1.0,
				   1.0,  0.0, 0.0, 1.0;
	}
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
GeomElem * Lin2Maker() { return new Lin2(); }

// Register element
int Lin2Register() { GeomElemFactory["Lin2"] = Lin2Maker; return 0; }

// Call register
int __Lin2_dummy_int  = Lin2Register();

}; // namespace FEM

#endif // MECHSYS_FEM_LIN2_H
