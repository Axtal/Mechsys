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

// STL
#include <iostream>

// MechSys
#include "fem/data.h"
#include "fem/node.h"
#include "fem/elems/quad8pstrain.h"
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "util/exception.h"

using FEM::Nodes;
using FEM::Elems;
using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	double H  = 2.0;   // height
	double L  = 2.0;   // length
	double E  = 207.0; // Young
	double nu = 0.3;   // Poisson
	double q  = 1.0;   // Load

	/*        | | | | | | | | | | |  q=1
	          V V V V V V V V V V V 
	  ---   3 @---------@---------@ 2
	   |      |         6         |
	   |      |                   |
	   |      |                   |
	          |                   |
	   H    7 @        e[0]       @ 5
	          |                   |
	   |      |                   |
	   |      |                   |
	   |      |         4         |
	  ---   0 @---------@---------@ 1
	         /_\       /_\       /_\
	         o o       ///       o o
	
	          |-------- L --------|
	 */

	// Input
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);
	else cout << "[1;32mYou may call this program as in:\t " << argv[0] << " LinSol\n  where LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n [0m[1;34m Now using LA (LAPACK)\n[0m" << endl;

	// 0) Problem dimension
	FEM::Dim = 2; // 2D

	// 1) Nodes
	FEM::AddNode(  0.0 ,   0.0); // 0
	FEM::AddNode(    L ,   0.0); // 1
	FEM::AddNode(    L ,     H); // 2
	FEM::AddNode(  0.0 ,     H); // 3
	FEM::AddNode(L/2.0 ,   0.0); // 4
	FEM::AddNode(    L , H/2.0); // 5
	FEM::AddNode(L/2.0 ,     H); // 6
	FEM::AddNode(  0.0 , H/2.0); // 7

	// 2) Elements
	FEM::AddElem("Quad8PStrain", /*IsActive*/true);

	// 3) Set connectivity (list of nodes must be LOCAL)
	Elems[0]->SetNode(0,0)->SetNode(1,1)->SetNode(2,2)->SetNode(3,3)->SetNode(4,4)->SetNode(5,5)->SetNode(6,6)->SetNode(7,7);

	// 4) Boundary conditions (must be after connectivity)
	Nodes[0]->Bry("uy",0.0);
	Nodes[4]->Bry("uy",0.0)->Bry("ux",0.0);
	Nodes[1]->Bry("uy",0.0);
	Elems[0]->Bry("fy",-q, 3, 2,3,6); // Actually, fy is traction == ty (list of nodes must be LOCAL)

	// 5) Parameters and initial values
	String prms; prms.Printf("E=%f  nu=%f",E,nu);
	Elems[0]->SetModel("LinElastic", prms.GetSTL().c_str(), "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0");

	// Stiffness
	Array<size_t>          map;
	Array<bool>            pre;
	LinAlg::Matrix<double> Ke0;
	Elems[0]->Order1Matrix(0,Ke0);
	//cout << "Ke0=\n" << Ke0 << endl;

	// 6) Solve
	FEM::Solver * sol = FEM::AllocSolver("AutoME");
	sol -> SetLinSol(linsol.GetSTL().c_str()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();

	// Output
	cout << "Node 3: ux = " << Nodes[3]->Val("ux") << " : uy = " << Nodes[3]->Val("uy") << " : fy = "  << Nodes[3]->Val("fy")  << endl;
	cout << "Node 6: ux = " << Nodes[6]->Val("ux") << " : uy = " << Nodes[6]->Val("uy") << " : fy = "  << Nodes[6]->Val("fy")  << endl;
	cout << "Node 2: ux = " << Nodes[2]->Val("ux") << " : uy = " << Nodes[2]->Val("uy") << " : fy = "  << Nodes[2]->Val("fy")  << endl << endl;;
	cout << "Node 0: ux = " << Nodes[0]->Val("ux") << " : uy = " << Nodes[0]->Val("uy") << " : fy = "  << Nodes[0]->Val("fy")  << endl;
	cout << "Node 4: ux = " << Nodes[4]->Val("ux") << " : uy = " << Nodes[4]->Val("uy") << " : fy = "  << Nodes[4]->Val("fy")  << endl;
	cout << "Node 1: ux = " << Nodes[1]->Val("ux") << " : uy = " << Nodes[1]->Val("uy") << " : fy = "  << Nodes[1]->Val("fy")  << endl << endl;;
	cout << "Elem 0: Sx = " << Elems[0]->Val("Sx") << " : Sy = " << Elems[0]->Val("Sy") << " : Sxy = " << Elems[0]->Val("Sxy") << endl;
	cout << "Elem 0: Ex = " << Elems[0]->Val("Ex") << " : Ey = " << Elems[0]->Val("Ey") << " : Exy = " << Elems[0]->Val("Exy") << endl;

	// Check
    double errors = 0.0;

	double Sy = q;
	double Ex = -nu*(1.0+nu)*Sy/E;
	double Ey =  (1.0-nu*nu)*Sy/E;

	errors += fabs(Elems[0]->Val("Ex" ) - (Ex));
	errors += fabs(Elems[0]->Val("Ey" ) - (Ey));
	errors += fabs(Elems[0]->Val("Exy") - (0.0));
	errors += fabs(Elems[0]->Val("Sx" ) - (0.0));
	errors += fabs(Elems[0]->Val("Sy" ) - (Sy ));
	errors += fabs(Elems[0]->Val("Sxy") - (0.0));

	errors += fabs(Nodes[0]->Val("ux") - ( 0.5*L*Ex));
	errors += fabs(Nodes[1]->Val("ux") - (-0.5*L*Ex));
	errors += fabs(Nodes[2]->Val("ux") - (-0.5*L*Ex));
	errors += fabs(Nodes[3]->Val("ux") - ( 0.5*L*Ex));
	errors += fabs(Nodes[4]->Val("ux") - (      0.0));
	errors += fabs(Nodes[5]->Val("ux") - (-0.5*L*Ex));
	errors += fabs(Nodes[6]->Val("ux") - (      0.0));
	errors += fabs(Nodes[7]->Val("ux") - ( 0.5*L*Ex));

	errors += fabs(Nodes[0]->Val("uy") - (      0.0));
	errors += fabs(Nodes[1]->Val("uy") - (      0.0));
	errors += fabs(Nodes[2]->Val("uy") - (    -H*Ey));
	errors += fabs(Nodes[3]->Val("uy") - (    -H*Ey));
	errors += fabs(Nodes[4]->Val("uy") - (      0.0));
	errors += fabs(Nodes[5]->Val("uy") - (-0.5*H*Ey));
	errors += fabs(Nodes[6]->Val("uy") - (    -H*Ey));
	errors += fabs(Nodes[7]->Val("uy") - (-0.5*H*Ey));

	errors += fabs(Nodes[3]->Val("fy") - (    -q*L/6.0));
	errors += fabs(Nodes[6]->Val("fy") - (-2.0*q*L/3.0));
	errors += fabs(Nodes[2]->Val("fy") - (    -q*L/6.0));

	errors += fabs(Nodes[3]->Val("fy")+Nodes[6]->Val("fy")+Nodes[2]->Val("fy")-(-q*L));

	if (fabs(errors)>1.0e-13) cout << "[1;31m\nErrors(" << linsol << ") = " << errors << "[0m\n" << endl;
	else                      cout << "[1;32m\nErrors(" << linsol << ") = " << errors << "[0m\n" << endl;

	return 0;
}
catch (Exception * e) 
{
	e->Cout();
	if (e->IsFatal()) {delete e; exit(1);}
	delete e;
}
catch (char const * m)
{
	std::cout << "Fatal: " << m << std::endl;
	exit (1);
}
catch (...)
{
	std::cout << "Some exception (...) ocurred\n";
} 
