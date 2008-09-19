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
#include "fem/element.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "vtkCellType.h"

namespace FEM
{

// Quad8 Constants
const int QUAD8_NNODES      = 8;
const int QUAD8_NINTPTS     = 4;
const int QUAD8_NFACENODES  = 3;
const int QUAD8_NFACEINTPTS = 2;
const Element::IntegPoint QUAD8_INTPTS[]=
{{ -sqrt(3.0)/3.0, -sqrt(3.0)/3.0, 0.0, 1.0 },
 {  sqrt(3.0)/3.0, -sqrt(3.0)/3.0, 0.0, 1.0 },
 { -sqrt(3.0)/3.0,  sqrt(3.0)/3.0, 0.0, 1.0 },
 {  sqrt(3.0)/3.0,  sqrt(3.0)/3.0, 0.0, 1.0 }};
const Element::IntegPoint QUAD8_FACEINTPTS[]=
{{ -sqrt(3.0)/3.0, 0.0, 0.0, 1.0 },
 {  sqrt(3.0)/3.0, 0.0, 0.0, 1.0 }};

class Quad8: public virtual Element
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

	// Destructor
	virtual ~Quad8() {}

	// Derived methods
	int  VTKCellType  () const { return VTK_QUADRATIC_QUAD; }
	void VTKConnect   (String & Nodes) const;
	void GetFaceNodes (int FaceID, Array<Node*> & FaceConnects) const;
	void Shape        (double r, double s, double t, LinAlg::Vector<double> & Shape)  const;
	void Derivs       (double r, double s, double t, LinAlg::Matrix<double> & Derivs) const;
	void FaceShape    (double r, double s, LinAlg::Vector<double> & FaceShape)  const;
	void FaceDerivs   (double r, double s, LinAlg::Matrix<double> & FaceDerivs) const;

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
Quad8::FaceMap Quad8::Face2Node[]= {{ 0, 3, 7 },
                                    { 1, 2, 5 },
                                    { 0, 1, 4 },
                                    { 2, 3, 6 }};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Quad8::Quad8()
{
	// Setup nodes number
	_n_nodes        = QUAD8_NNODES;
	_n_int_pts      = QUAD8_NINTPTS;
	_n_face_nodes   = QUAD8_NFACENODES;
	_n_face_int_pts = QUAD8_NFACEINTPTS;

	// Allocate nodes (connectivity)
	_connects.Resize(_n_nodes);
	_connects.SetValues(NULL);

	// Setup pointer to the array of Integration Points
	_a_int_pts      = QUAD8_INTPTS;
	_a_face_int_pts = QUAD8_FACEINTPTS;
}

inline void Quad8::VTKConnect(String & Nodes) const
{
	Nodes.Printf("%d %d %d %d %d %d %d %d",_connects[0]->GetID(),
	                                       _connects[1]->GetID(),
	                                       _connects[2]->GetID(),
	                                       _connects[3]->GetID(),
	                                       _connects[4]->GetID(),
	                                       _connects[5]->GetID(),
	                                       _connects[6]->GetID(),
	                                       _connects[7]->GetID());
}

inline void Quad8::GetFaceNodes(int FaceID, Array<Node*> & FaceConnects) const
{
	FaceConnects.Resize(3);
	FaceConnects[0] = _connects[Face2Node[FaceID].L];
	FaceConnects[1] = _connects[Face2Node[FaceID].R];
	FaceConnects[2] = _connects[Face2Node[FaceID].M];
}

inline void Quad8::Shape(double r, double s, double t, LinAlg::Vector<double> & Shape) const
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
	Shape.Resize (QUAD8_NNODES);

	double rp1=1.0+r; double rm1=1.0-r;
	double sp1=1.0+s; double sm1=1.0-s;

	Shape(0) = 0.25*rm1*sm1*(rm1+sm1-3.0);
	Shape(1) = 0.25*rp1*sm1*(rp1+sm1-3.0);
	Shape(2) = 0.25*rp1*sp1*(rp1+sp1-3.0);
	Shape(3) = 0.25*rm1*sp1*(rm1+sp1-3.0);
	Shape(4) = 0.50*sm1*(1.0-r*r);
	Shape(5) = 0.50*rp1*(1.0-s*s);
	Shape(6) = 0.50*sp1*(1.0-r*r);
	Shape(7) = 0.50*rm1*(1.0-s*s);
}

inline void Quad8::Derivs(double r, double s, double t, LinAlg::Matrix<double> & Derivs) const
{
	/*           _     _ T
	 *          |  dNi  |
	 * Derivs = |  ---  |   , where cj = r, s
	 *          |_ dcj _|
	 *
	 * Derivs(j,i), j=>local coordinate and i=>shape function
	 */
	Derivs.Resize (2, QUAD8_NNODES);

	double rp1=1.0+r; double rm1=1.0-r;
	double sp1=1.0+s; double sm1=1.0-s;

	Derivs(0,0) = - 0.25 * sm1 * (rm1 + rm1 + sm1 - 3.0);
	Derivs(0,1) =   0.25 * sm1 * (rp1 + rp1 + sm1 - 3.0);
	Derivs(0,2) =   0.25 * sp1 * (rp1 + rp1 + sp1 - 3.0);
	Derivs(0,3) = - 0.25 * sp1 * (rm1 + rm1 + sp1 - 3.0);
	Derivs(0,4) = - r * sm1;
	Derivs(0,5) =   0.50 * (1.0 - s * s);
	Derivs(0,6) = - r * sp1;
	Derivs(0,7) = - 0.5 * (1.0 - s * s);

	Derivs(1,0) = - 0.25 * rm1 * (sm1 + rm1 + sm1 - 3.0);
	Derivs(1,1) = - 0.25 * rp1 * (sm1 + rp1 + sm1 - 3.0);
	Derivs(1,2) =   0.25 * rp1 * (sp1 + rp1 + sp1 - 3.0);
	Derivs(1,3) =   0.25 * rm1 * (sp1 + rm1 + sp1 - 3.0);
	Derivs(1,4) = - 0.50 * (1.0 - r * r);
	Derivs(1,5) = - s * rp1;
	Derivs(1,6) =   0.50 * (1.0 - r * r);
	Derivs(1,7) = - s * rm1;
}

inline void Quad8::FaceShape(double r, double s, LinAlg::Vector<double> & FaceShape) const
{
	/*
	 *       @-----------@-----------@-> r
	 *       0           2           1
	 *       |           |           |
	 *      r=-1         r=0        r=+1
	 */
	FaceShape.Resize(QUAD8_NFACENODES);
	FaceShape(0) = 0.5 * (r*r-r);
	FaceShape(1) = 0.5 * (r*r+r);
	FaceShape(2) = 1.0 -  r*r;
}

inline void Quad8::FaceDerivs(double r, double s, LinAlg::Matrix<double> & FaceDerivs) const
{
	/*           _     _ T
	 *          |  dNi  |
	 * Derivs = |  ---  |   , where cj = r, s
	 *          |_ dcj _|
	 *
	 * Derivs(j,i), j=>local coordinate and i=>shape function
	 */
	FaceDerivs.Resize(1,QUAD8_NFACENODES);
	FaceDerivs(0,0) =  r  - 0.5;
	FaceDerivs(0,1) =  r  + 0.5;
	FaceDerivs(0,2) = -2.0* r;
}


}; // namespace FEM

#endif // MECHSYS_FEM_QUAD8_H
