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
#include "fem/element.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"

namespace FEM
{

// Tri6 Constants
const int TRI6_NNODES      = 6;
const int TRI6_NINTPTS     = 3;
const int TRI6_NFACENODES  = 3;
const int TRI6_NFACEINTPTS = 2;
const Element::IntegPoint TRI6_INTPTS[]=
{{  0.166666666666666666666667,  0.166666666666666666666667,  0.0,  0.166666666666666666666667 }, 
 {  0.666666666666666666666667,  0.166666666666666666666667,  0.0,  0.166666666666666666666667 },
 {  0.166666666666666666666667,  0.666666666666666666666667,  0.0,  0.166666666666666666666667 }};
const Element::IntegPoint TRI6_FACEINTPTS[]=
{{ -0.577350269189625764509149, 0.0, 0.0, 1.0 },
 {  0.577350269189625764509149, 0.0, 0.0, 1.0 }};

class Tri6: public virtual Element
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

	// Destructor
	virtual ~Tri6() {}

	// Derived methods
	int  VTKCellType  () const { return 22; } // VTK_QUADRATIC_TRIANGLE
	void VTKConnect   (String & Nodes) const;
	void GetFaceNodes (int FaceID, Array<Node*> & FaceConnects) const;
	void Shape        (double r, double s, double t, LinAlg::Vector<double> & Shape)  const;
	void Derivs       (double r, double s, double t, LinAlg::Matrix<double> & Derivs) const;
	void FaceShape    (double r, double s, LinAlg::Vector<double> & FaceShape)  const;
	void FaceDerivs   (double r, double s, LinAlg::Matrix<double> & FaceDerivs) const;

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
                                  { 0, 2, 5 }};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Tri6::Tri6()
{
	// Setup nodes number
	_n_nodes        = TRI6_NNODES;
	_n_int_pts      = TRI6_NINTPTS;
	_n_face_nodes   = TRI6_NFACENODES;
	_n_face_int_pts = TRI6_NFACEINTPTS;

	// Allocate nodes (connectivity)
	_connects.Resize(_n_nodes);

	// Setup pointer to the array of Integration Points
	_a_int_pts      = TRI6_INTPTS;
	_a_face_int_pts = TRI6_FACEINTPTS;
}

inline void Tri6::VTKConnect(String & Nodes) const
{
	Nodes.Printf("%d %d %d %d %d %d",_connects[0]->GetID(),
	                                 _connects[1]->GetID(),
	                                 _connects[2]->GetID(),
	                                 _connects[3]->GetID(),
	                                 _connects[4]->GetID(),
	                                 _connects[5]->GetID());
}

inline void Tri6::GetFaceNodes(int FaceID, Array<Node*> & FaceConnects) const
{
	FaceConnects.Resize(3);
	FaceConnects[0] = _connects[Face2Node[FaceID].L];
	FaceConnects[1] = _connects[Face2Node[FaceID].R];
	FaceConnects[2] = _connects[Face2Node[FaceID].M];
}

inline void Tri6::Shape(double r, double s, double t, LinAlg::Vector<double> & Shape) const
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
	Shape.Resize(TRI6_NNODES);
    Shape(0) = 1.0-(r+s)*(3.0-2.0*(r+s));
    Shape(1) = r*(2.0*r-1.0);
    Shape(2) = s*(2.0*s-1.0);
    Shape(3) = 4.0*r*(1.0-(r+s));
    Shape(4) = 4.0*r*s;
    Shape(5) = 4.0*s*(1.0-(r+s));
}

inline void Tri6::Derivs(double r, double s, double t, LinAlg::Matrix<double> & Derivs) const
{
	/*           _     _ T
	 *          |  dNi  |
	 * Derivs = |  ---  |   , where cj = r, s
	 *          |_ dcj _|
	 *
	 * Derivs(j,i), j=>local coordinate and i=>shape function
	 */
	Derivs.Resize(2, TRI6_NNODES);

    Derivs(0,0) = -3.0 + 4.0 * (r + s);       Derivs(1,0) = -3.0 + 4.0*(r + s);
    Derivs(0,1) =  4.0 * r - 1.;              Derivs(1,1) =  0.0 ; 
	Derivs(0,2) =  0.0;                       Derivs(1,2) =  4.0 * s - 1.0;   
    Derivs(0,3) =  4.0 - 8.0 * r - 4.0 * s;   Derivs(1,3) = -4.0 * r;
    Derivs(0,4) =  4.0 * s;                   Derivs(1,4) =  4.0 * r;
    Derivs(0,5) = -4.0 * s;                   Derivs(1,5) =  4.0 - 4.0 * r - 8.0*s;
}

inline void Tri6::FaceShape(double r, double s, LinAlg::Vector<double> & FaceShape) const
{
	/*  
	 *  
	 *       @-----------@-----------@-> r
	 *       0           1           2
	 */
	FaceShape.Resize(TRI6_NFACENODES);
	FaceShape(0) = 0.5 * (r*r-r);
	FaceShape(1) = 1.0 -  r*r;
	FaceShape(2) = 0.5 * (r*r+r);
}

inline void Tri6::FaceDerivs(double r, double s, LinAlg::Matrix<double> & FaceDerivs) const
{
	/*           _     _ T
	 *          |  dNi  |
	 * Derivs = |  ---  |   , where cj = r, s
	 *          |_ dcj _|
	 *
	 * Derivs(j,i), j=>local coordinate and i=>shape function
	 */
	FaceDerivs.Resize(1,TRI6_NFACENODES);
	FaceDerivs(0,0) =  r  - 0.5;
	FaceDerivs(0,1) = -2.0* r;
	FaceDerivs(0,2) =  r  + 0.5;
}


}; // namespace FEM

#endif // MECHSYS_FEM_TRI6_H
