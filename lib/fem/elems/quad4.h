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
#include "fem/element.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "vtkCellType.h"

namespace FEM
{

class Quad4: public virtual Element
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
	Quad4();

	// Destructor
	virtual ~Quad4() {}

	// Derived methods
	void SetIntPoints (int NumGaussPoints1D);
	int  VTKCellType  () const { return VTK_QUAD; }
	void VTKConnect   (String & Nodes) const;
	void GetFaceNodes (int FaceID, Array<Node*> & FaceConnects) const;
	void Shape        (double r, double s, double t, LinAlg::Vector<double> & Shape)  const;
	void Derivs       (double r, double s, double t, LinAlg::Matrix<double> & Derivs) const;
	void FaceShape    (double r, double s, LinAlg::Vector<double> & FaceShape)  const;
	void FaceDerivs   (double r, double s, LinAlg::Matrix<double> & FaceDerivs) const;

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
Quad4::FaceMap Quad4::Face2Node[]= {{ 0, 3 },
                                    { 1, 2 },
                                    { 0, 1 },
                                    { 2, 3 }};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Quad4::Quad4()
{
	// Setup nodes number
	_n_nodes      = 4;
	_n_face_nodes = 2;

	// Allocate nodes (connectivity)
	_connects.Resize    (_n_nodes);
	_connects.SetValues (NULL);

	// Integration points and Extrapolation Matrix
	SetIntPoints (/*NumGaussPoints1D*/2);
}

inline void Quad4::SetIntPoints(int NumGaussPoints1D)
{
	// Set IPs
	FEM::SetGaussIP (/*NDim*/2, NumGaussPoints1D, _a_int_pts);
	FEM::SetGaussIP (/*NDim*/1, NumGaussPoints1D, _a_face_int_pts);

	// Evaluation matrix [E] (see pag. 629 of Burnett, D. S. (1988). Finite Element Analysis: From Concepts to Applications. Addison-Wesley. 844p.)
	// [nodalvalues]T = [E] * [c0,c1,c2]T
	//  nodalvalue_i  = c0 + c1*r_i + c2*s_i
	// (i:node number)
	LinAlg::Matrix<double> eval_mat;
	eval_mat.Resize(_n_nodes, 3);
	eval_mat = 1.0, -1.0, -1.0,
	           1.0,  1.0, -1.0,
	           1.0,  1.0,  1.0,
	           1.0, -1.0,  1.0;

	// Set extrapolation matrix
	FEM::SetExtrapMatrix (/*NDim*/2, _a_int_pts, eval_mat, _extrap_mat);
}

inline void Quad4::VTKConnect(String & Nodes) const
{
	Nodes.Printf("%d %d %d %d",_connects[0]->GetID(),
	                           _connects[1]->GetID(),
	                           _connects[2]->GetID(),
	                           _connects[3]->GetID());
}

inline void Quad4::GetFaceNodes(int FaceID, Array<Node*> & FaceConnects) const
{
	FaceConnects.Resize(2);
	FaceConnects[0] = _connects[Face2Node[FaceID].L];
	FaceConnects[1] = _connects[Face2Node[FaceID].R];
}

inline void Quad4::Shape(double r, double s, double t, LinAlg::Vector<double> & Shape) const
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
	Shape.Resize (/*NNodes*/4);

	double rp=1.0+r; double rm=1.0-r;
	double sp=1.0+s; double sm=1.0-s;

	Shape(0) = 0.25*rm*sm;
	Shape(1) = 0.25*rp*sm;
	Shape(2) = 0.25*rp*sp;
	Shape(3) = 0.25*rm*sp;
}

inline void Quad4::Derivs(double r, double s, double t, LinAlg::Matrix<double> & Derivs) const
{
	/*           _     _ T
	 *          |  dNi  |
	 * Derivs = |  ---  |   , where cj = r, s
	 *          |_ dcj _|
	 *
	 * Derivs(j,i), j=>local coordinate and i=>shape function
	 */
	Derivs.Resize (2, /*NNodes*/4);

	double rp=1.0+r; double rm=1.0-r;
	double sp=1.0+s; double sm=1.0-s;

	Derivs(0,0) = -0.25*sm;   Derivs(1,0) = -0.25*rm;
	Derivs(0,1) =  0.25*sm;   Derivs(1,1) = -0.25*rp;
	Derivs(0,2) =  0.25*sp;   Derivs(1,2) =  0.25*rp;
	Derivs(0,3) = -0.25*sp;   Derivs(1,3) =  0.25*rm;
}

inline void Quad4::FaceShape(double r, double s, LinAlg::Vector<double> & FaceShape) const
{
	/*  
	 *       0           |           1
	 *       @-----------+-----------@-> r
	 *      -1           |          +1
	 */
	FaceShape.Resize(/*NumFaceNodes*/2);
	FaceShape(0) = 0.5*(1.0-r);
	FaceShape(1) = 0.5*(1.0+r);
}

inline void Quad4::FaceDerivs(double r, double s, LinAlg::Matrix<double> & FaceDerivs) const
{
	/*           _     _ T
	 *          |  dNi  |
	 * Derivs = |  ---  |   , where cj = r, s
	 *          |_ dcj _|
	 *
	 * Derivs(j,i), j=>local coordinate and i=>shape function
	 */
	FaceDerivs.Resize(1,/*NumFaceNodes*/2);
	FaceDerivs(0,0) = -0.5;
	FaceDerivs(0,1) =  0.5;
}

}; // namespace FEM

#endif // MECHSYS_FEM_QUAD4_H
