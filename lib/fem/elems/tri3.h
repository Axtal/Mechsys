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

#ifndef MECHSYS_FEM_TRI3_H
#define MECHSYS_FEM_TRI3_H

// MechSys
#include "fem/element.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "vtkCellType.h"

namespace FEM
{

class Tri3: public virtual Element
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

	// Destructor
	virtual ~Tri3() {}

	// Derived methods
	void SetIntPoints (int NumGaussPoints1D);
	int  VTKCellType  () const { return VTK_TRIANGLE; }
	void VTKConnect   (String & Nodes) const;
	void GetFaceNodes (int FaceID, Array<Node*> & FaceConnects) const;
	void Shape        (double r, double s, double t, LinAlg::Vector<double> & Shape)  const;
	void Derivs       (double r, double s, double t, LinAlg::Matrix<double> & Derivs) const;
	void FaceShape    (double r, double s, LinAlg::Vector<double> & FaceShape)  const;
	void FaceDerivs   (double r, double s, LinAlg::Matrix<double> & FaceDerivs) const;
	void LocalCoords  (LinAlg::Matrix<double> & coords) const;

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
	_n_nodes      = 3;
	_n_face_nodes = 2;

	// Allocate nodes (connectivity)
	_connects.Resize    (_n_nodes);
	_connects.SetValues (NULL);

	// Integration Points and Extrapolation Matrix
	SetIntPoints (/*NumGaussPointsTotal*/3);
}

inline void Tri3::SetIntPoints(int NumGaussPointsTotal)
{
	// Setup pointer to the array of Integration Points
	if      (NumGaussPointsTotal==3)  _a_int_pts = TRI_IP3;
	else if (NumGaussPointsTotal==4)  _a_int_pts = TRI_IP4;
	else if (NumGaussPointsTotal==6)  _a_int_pts = TRI_IP6;
	else if (NumGaussPointsTotal==7)  _a_int_pts = TRI_IP7;
	else if (NumGaussPointsTotal==13) _a_int_pts = TRI_IP13;
	else throw new Fatal("tri3::SetIntPoints: Error in number of integration points.");

	_n_int_pts      = NumGaussPointsTotal;
	_a_face_int_pts = LIN_IP2;
	_n_face_int_pts = 2;
}

inline void Tri3::LocalCoords(LinAlg::Matrix<double> & coords) const 
{
	coords.Resize(3,3);
	coords = 0.0, 0.0, 1.0,
	         1.0, 0.0, 1.0,
	         0.0, 1.0, 1.0;
}

inline void Tri3::VTKConnect(String & Nodes) const
{
	Nodes.Printf("%d %d %d",_connects[0]->GetID(),
	                        _connects[1]->GetID(),
	                        _connects[2]->GetID());
}

inline void Tri3::GetFaceNodes(int FaceID, Array<Node*> & FaceConnects) const
{
	FaceConnects.Resize(2);
	FaceConnects[0] = _connects[Face2Node[FaceID].L];
	FaceConnects[1] = _connects[Face2Node[FaceID].R];
}

inline void Tri3::Shape(double r, double s, double t, LinAlg::Vector<double> & Shape) const
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
	Shape.Resize(/*NumNodes*/3);
    Shape(0) = 1.0-r-s;
    Shape(1) = r;
    Shape(2) = s;
}

inline void Tri3::Derivs(double r, double s, double t, LinAlg::Matrix<double> & Derivs) const
{
	/*           _     _ T
	 *          |  dNi  |
	 * Derivs = |  ---  |   , where cj = r, s
	 *          |_ dcj _|
	 *
	 * Derivs(j,i), j=>local coordinate and i=>shape function
	 */
	Derivs.Resize(2, /*NumNodes*/3);
    Derivs(0,0) = -1.0;    Derivs(1,0) = -1.0;
    Derivs(0,1) =  1.0;    Derivs(1,1) =  0.0;
	Derivs(0,2) =  0.0;    Derivs(1,2) =  1.0;
}

inline void Tri3::FaceShape(double r, double s, LinAlg::Vector<double> & FaceShape) const
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

inline void Tri3::FaceDerivs(double r, double s, LinAlg::Matrix<double> & FaceDerivs) const
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

#endif // MECHSYS_FEM_TRI3_H
