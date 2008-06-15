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
#include "fem/elems/quad8equilib.h"
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
	/*        | | | | | | | | | | |  q=1
	          V V V V V V V V V V V 
	        3 @---------@---------@ 2
	          |         6         |
	          |                   |
	          |                   |
	          |                   |
	        7 @        e[0]       @ 5
	          |                   |
	          |                   |
	          |                   |
	          |         4         |
	        0 @---------@---------@ 1
	         /_\       /_\       /_\
	         o o       ///       o o
	 */

	// Input
	if (argc!=2) throw new Message(_("Please, call this program as in:\n\t\t %s LinSol\n  where:\n   LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n"),argv[0]);

	// 0) Geometry type
	FEM::GeometryType = 2; // 2D(plane-strain)

	// 1) Nodes
	FEM::AddNode(0.0, 0.0); // 0
	FEM::AddNode(2.0, 0.0); // 1
	FEM::AddNode(2.0, 2.0); // 2
	FEM::AddNode(0.0, 2.0); // 3
	FEM::AddNode(1.0, 0.0); // 4
	FEM::AddNode(2.0, 1.0); // 5
	FEM::AddNode(1.0, 2.0); // 6
	FEM::AddNode(0.0, 1.0); // 7

	// 2) Elements
	FEM::AddElem("Quad8Equilib", /*IsActive*/true);

	// 3) Set connectivity (list of nodes must be LOCAL)
	Elems[0]->SetNode(0,0)->SetNode(1,1)->SetNode(2,2)->SetNode(3,3)->SetNode(4,4)->SetNode(5,5)->SetNode(6,6)->SetNode(7,7);

	// 4) Boundary conditions (must be after connectivity)
	Nodes[0]->Bry("uy",0.0);
	Nodes[4]->Bry("uy",0.0)->Bry("ux",0.0);
	Nodes[1]->Bry("uy",0.0);
	Elems[0]->Bry("fy",-1.0, 3, 3,6,2); // Actually, fy is traction == ty (list of nodes must be LOCAL)

	// 5) Parameters and initial values
	Elems[0]->SetModel("LinElastic", "E=207.0 nu=0.3", "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0");

	// Stiffness
	Array<size_t>          map;
	Array<bool>            pre;
	LinAlg::Matrix<double> Ke0;
	Elems[0]->Order1Matrix(0,Ke0);
	cout << "Ke0=\n" << Ke0 << endl;

	// 6) Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol -> SetLinSol(argv[1]) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();

	// Output
	LinAlg::Matrix<double> values;
	Array<String>          labels;
	Elems[0]->OutNodes(values,labels);
	std::cout << labels << std::endl;
	std::cout << values << std::endl;

	// Check
    double errors = 1.0;
	/*
	errors += fabs(Elems[0]->Val(0, "Sx") - ( 1.56432140e-01));
	errors += fabs(Elems[0]->Val(1, "Sx") - (-3.00686928e-01));
	errors += fabs(Elems[0]->Val(2, "Sx") - ( 1.44254788e-01));
	errors += fabs(Elems[0]->Val(3, "Sx") - (-3.19109076e-01));
	errors += fabs(Elems[0]->Val(4, "Sx") - (-3.31286428e-01));
	errors += fabs(Elems[0]->Val(5, "Sx") - ( 1.25832639e-01));
	*/
	if (fabs(errors)>1.0e-7) cout << "[1;31mErrors(" << argv[1] << ") = " << errors << "[0m\n" << endl;
	else                     cout << "[1;32mErrors(" << argv[1] << ") = " << errors << "[0m\n" << endl;

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
