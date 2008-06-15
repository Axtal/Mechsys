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
#include "fem/elems/tri6equilib.h"
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
	/*       F1        F2       F3    F1=F2=F3 = 1.0
	          ^         ^         ^
	          |         |         |
	
	        2 @---------@---------@ 6
	          |',       7         |
	          |  ',        e[1]   |
	          |    ',             |
	          |      ',  4        |
	        5 @        '@         @ 8
	          |          ',       |
	          |   e[0]     ',     |
	          |              ',   |
	          |         3      ', |
	        0 @---------@---------@ 1
	         /_\       /_\       /_\
	         ///       o o       o o
	 */

	// Input
	if (argc!=2) throw new Message(_("Please, call this program as in:\n\t\t %s LinSol\n  where:\n   LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n"),argv[0]);

	// 0) Geometry type
	FEM::GeometryType = 2; // 2D(plane-strain)

	// 1) Nodes
	FEM::AddNode(0.0, 0.0); // 0
	FEM::AddNode(1.0, 0.0); // 1
	FEM::AddNode(0.0, 1.0); // 2
	FEM::AddNode(0.5, 0.0); // 3
	FEM::AddNode(0.5, 0.5); // 4
	FEM::AddNode(0.0, 0.5); // 5
	FEM::AddNode(1.0, 1.0); // 6
	FEM::AddNode(0.5, 1.0); // 7
	FEM::AddNode(1.0, 0.5); // 8

	// 2) Elements
	FEM::AddElem("Tri6Equilib", /*IsActive*/true);
	FEM::AddElem("Tri6Equilib", /*IsActive*/true);

	// 3) Set connectivity
	Elems[0]->SetNode(0,0)->SetNode(1,1)->SetNode(2,2)->SetNode(3,3)->SetNode(4,4)->SetNode(5,5);
	Elems[1]->SetNode(0,6)->SetNode(1,2)->SetNode(2,1)->SetNode(3,7)->SetNode(4,4)->SetNode(5,8);

	// 4) Boundary conditions (must be after connectivity)
	Nodes[0]->Bry("ux",0.0)->Bry("uy",0.0);
	Nodes[1]->Bry("uy",0.0);
	Nodes[3]->Bry("uy",0.0);
	Nodes[2]->Bry("fy",1.0);
	Nodes[7]->Bry("fy",1.0);
	Nodes[6]->Bry("fy",1.0);

	// 5) Parameters and initial values
	Elems[0]->SetModel("LinElastic", "E=10000.0 nu=0.25", "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0");
	Elems[1]->SetModel("LinElastic", "E=10000.0 nu=0.25", "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0");

	// Stiffness
	Array<size_t>          map;
	Array<bool>            pre;
	LinAlg::Matrix<double> Ke0;
	LinAlg::Matrix<double> Ke1;
	Elems[0]->Order1Matrix(0,Ke0);
	Elems[1]->Order1Matrix(0,Ke1);
	cout << "Ke0=\n" << Ke0 << endl;
	cout << "Ke1=\n" << Ke1 << endl;

	// 6) Solve
	//FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	FEM::Solver * sol = FEM::AllocSolver("AutoME");
	sol -> SetLinSol(argv[1]) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();

	// Output
	LinAlg::Matrix<double> values0;
	LinAlg::Matrix<double> values1;
	Array<String>          labels0;
	Array<String>          labels1;
	Elems[0]->OutNodes(values0, labels0);
	Elems[1]->OutNodes(values1, labels1);
	std::cout << labels0 << std::endl;
	std::cout << values0 << std::endl;
	std::cout << labels1 << std::endl;
	std::cout << values1 << std::endl;

	// Check
    double errors = 0.0;

	// Element 0 : using Element::Val
	errors += fabs(Elems[0]->Val(0, "Sx") - ( 1.56432140e-01));
	errors += fabs(Elems[0]->Val(1, "Sx") - (-3.00686928e-01));
	errors += fabs(Elems[0]->Val(2, "Sx") - ( 1.44254788e-01));
	errors += fabs(Elems[0]->Val(3, "Sx") - (-3.19109076e-01));
	errors += fabs(Elems[0]->Val(4, "Sx") - (-3.31286428e-01));
	errors += fabs(Elems[0]->Val(5, "Sx") - ( 1.25832639e-01));

	errors += fabs(Elems[0]->Val(0, "Sy") - (-2.05141549e-01));
	errors += fabs(Elems[0]->Val(1, "Sy") - ( 1.15872190e+00));
	errors += fabs(Elems[0]->Val(2, "Sy") - (-9.53580350e-01));
	errors += fabs(Elems[0]->Val(3, "Sy") - (-2.22127394e+00));
	errors += fabs(Elems[0]->Val(4, "Sy") - (-2.96971274e+00));
	errors += fabs(Elems[0]->Val(5, "Sy") - (-4.33357619e+00));

	errors += fabs(Elems[0]->Val(0, "Sxy") - (-1.56432140e-01));
	errors += fabs(Elems[0]->Val(1, "Sxy") - (-6.74437968e-02));
	errors += fabs(Elems[0]->Val(2, "Sxy") - ( 2.23875937e-01));
	errors += fabs(Elems[0]->Val(3, "Sxy") - (-4.90216486e-02));
	errors += fabs(Elems[0]->Val(4, "Sxy") - ( 3.31286428e-01));
	errors += fabs(Elems[0]->Val(5, "Sxy") - ( 2.42298085e-01));

	// Element 1 : using Element OutNodes
	errors += fabs(Elems[1]->Val(0, "Sx")  - ( 9.95732723e-01));
	errors += fabs(Elems[1]->Val(1, "Sx")  - ( 2.23875937e-01));
	errors += fabs(Elems[1]->Val(2, "Sx")  - (-1.21960866e+00));
	errors += fabs(Elems[1]->Val(3, "Sx")  - ( 1.39446295e+00));
	errors += fabs(Elems[1]->Val(4, "Sx")  - (-8.20878435e-01));
	errors += fabs(Elems[1]->Val(5, "Sx")  - (-4.90216486e-02));

	errors += fabs(Elems[1]->Val(0, "Sy")  - (-1.25426728e+00));
	errors += fabs(Elems[1]->Val(1, "Sy")  - ( 1.68328476e+00));
	errors += fabs(Elems[1]->Val(2, "Sy")  - (-4.29017485e-01));
	errors += fabs(Elems[1]->Val(3, "Sy")  - (-2.39612823e+00));
	errors += fabs(Elems[1]->Val(4, "Sy")  - (-1.57087843e+00));
	errors += fabs(Elems[1]->Val(5, "Sy")  - (-4.50843047e+00));

	errors += fabs(Elems[1]->Val(0, "Sxy") - (-9.95732723e-01));
	errors += fabs(Elems[1]->Val(1, "Sxy") - ( 1.29641965e+00));
	errors += fabs(Elems[1]->Val(2, "Sxy") - (-3.00686928e-01));
	errors += fabs(Elems[1]->Val(3, "Sxy") - ( 1.25832639e-01));
	errors += fabs(Elems[1]->Val(4, "Sxy") - ( 8.20878435e-01));
	errors += fabs(Elems[1]->Val(5, "Sxy") - (-1.47127394e+00));

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
