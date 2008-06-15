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
#include "fem/elems/quad4equilib.h"
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "util/exception.h"
#include "util/numstreams.h"

using FEM::Nodes;
using FEM::Elems;
using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	/*       F1                 F3    F1=F2=F3 = 1.0
	          ^                   ^
	          |                   |
	
	        3 @-------------------@ 2
	          |                   |
	          |                   |
	          |                   |
	          |                   |
	          |        e[0]       |
	          |                   |
	          |                   |
	          |                   |
	          |                   |
	        0 @-------------------@ 1
	         /_\                 /_\
	         ///                 o o
	 */

	// Input
	if (argc!=2) throw new Message(_("Please, call this program as in:\n\t\t %s LinSol\n  where:\n   LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n"),argv[0]);

	// 0) Geometry type
	FEM::GeometryType = 5; // 2D(plane-stress)

	// 1) Nodes
	FEM::AddNode(0.0, 0.0); // 0
	FEM::AddNode(1.0, 0.0); // 1
	FEM::AddNode(1.0, 1.0); // 2
	FEM::AddNode(0.0, 1.0); // 3

	// 2) Elements
	FEM::AddElem("Quad4Equilib", /*IsActive*/true);
	FEM::AddElem("Quad4Equilib", /*IsActive*/true);

	// 3) Set connectivity
	Elems[0]->SetNode(0,0)->SetNode(1,1)->SetNode(2,2)->SetNode(3,3);

	// 4) Boundary conditions (must be after connectivity)
	Nodes[0]->Bry("ux",0.0)->Bry("uy",0.0);
	Nodes[1]->Bry("uy",0.0);

	// 5) Parameters and initial values
	Elems[0]->SetModel("LinElastic", "E=96.0 nu=0.333333333333333333333333", "Sx=0.0 Sy=0.0 Sxy=0.0");

	// Stiffness
	Array<size_t>          map;
	Array<bool>            pre;
	LinAlg::Matrix<double> Ke0;  Ke0.SetNS(Util::_6_3);
	Elems[0]->Order1Matrix(0,Ke0);
	cout << "Ke0=\n" << Ke0 << endl;

	/*
	// 6) Solve
	//FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	FEM::Solver * sol = FEM::AllocSolver("AutoME");
	sol -> SetLinSol(argv[1]) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();

	// Check
    double errors = 1.0;
	*/
	/*
	errors += fabs(Elems[0]->Val(0, "Sx") - ( 1.56432140e-01));
	errors += fabs(Elems[0]->Val(1, "Sx") - (-3.00686928e-01));
	errors += fabs(Elems[0]->Val(2, "Sx") - ( 1.44254788e-01));
	errors += fabs(Elems[0]->Val(3, "Sx") - (-3.19109076e-01));
	errors += fabs(Elems[0]->Val(4, "Sx") - (-3.31286428e-01));
	errors += fabs(Elems[0]->Val(5, "Sx") - ( 1.25832639e-01));
	*/
	/*
	if (fabs(errors)>1.0e-7) cout << "[1;31mErrors(" << argv[1] << ") = " << errors << "[0m\n" << endl;
	else                     cout << "[1;32mErrors(" << argv[1] << ") = " << errors << "[0m\n" << endl;
	*/

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
