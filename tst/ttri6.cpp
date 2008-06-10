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
	/*                    
	          F1         F2        F3    F1=F2=F3 = 1.0
	          ^          ^         ^
	          |          |         | 

	          4          8
	          @----------@---------@ 7
	          |  \                 |
	          |    \               |
	          |      \             |
	          |        \           |
	        5 @          @ 3       @ 6
	          |            \       |   
	          |              \     |  
	          |                \   |
	          |                  \ |
	          @---------@----------@  
	          0         1          2  
	         /_\       /_\        /_\
			 ///       o o        o o
	 */

	// Input
	LinAlg::LinSol_T linsol;
	bool input_ok = false;
	if (argc==2) { if (atoi(argv[1])>=1 && atoi(argv[1])<=3) input_ok = true; }
	if (input_ok) linsol = static_cast<LinAlg::LinSol_T>(atoi(argv[1]));
	else throw new Message(_("Please, call this program as in:\n\t\t %s LinSol\n  where:\n   LinSol:\n \t1 => LAPACK_T  : DENSE\n \t2 => UMFPACK_T : SPARSE\n \t3 => SuperLU_T : SPARSE\n"),argv[0]);

	// 0) Geometry type
	FEM::GeometryType = 2; // 2D

	// 1) Nodes
	FEM::AddNode(0.0, 0.0); // 0
	FEM::AddNode(0.5, 0.0); // 1
	FEM::AddNode(1.0, 0.0); // 2
	FEM::AddNode(0.5, 0.5); // 3
	FEM::AddNode(0.0, 1.0); // 4
	FEM::AddNode(0.0, 0.5); // 5
	FEM::AddNode(1.0, 0.5); // 6
	FEM::AddNode(1.0, 1.0); // 7
	FEM::AddNode(0.5, 1.0); // 8

	// 2) Elements
	FEM::AddElem("Tri6Equilib", /*IsActive*/true);
	FEM::AddElem("Tri6Equilib", /*IsActive*/true);

	// 3) Set connectivity
	Elems[0]->SetNode(0, 0)->SetNode(1, 1)->SetNode(2, 2)->SetNode(3, 3)->SetNode(4, 4)->SetNode(5, 5);
	Elems[1]->SetNode(0, 7)->SetNode(1, 8)->SetNode(2, 4)->SetNode(3, 3)->SetNode(4, 2)->SetNode(5, 6);

	// 4) Boundary conditions (must be after connectivity)
	Nodes[0]->Bry("Dux", 0.0)->Bry("Duy", 0.0);
	Nodes[1]->Bry("Duy", 0.0);
	Nodes[2]->Bry("Duy", 0.0);
	Nodes[4]->Bry("Dfy", 1.0);
	Nodes[8]->Bry("Dfy", 1.0);
	Nodes[7]->Bry("Dfy", 1.0);

	// 5) Parameters and initial values
	Elems[0]->SetModel("LinElastic", "E=10000.0 nu=0.25", "Sx=0.0 Sy=0.0 Sxy=0.0");
	Elems[1]->SetModel("LinElastic", "E=10000.0 nu=0.25", "Sx=0.0 Sy=0.0 Sxy=0.0");

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
	sol->Solve(linsol, /*iStage*/0, /*NumDiv*/1, /*DeltaTime*/0);

	// Output
	String fn;
	     if (linsol==LinAlg::LAPACK_T  ) fn.append("LP");
	else if (linsol==LinAlg::UMFPACK_T ) fn.append("UM");
	else if (linsol==LinAlg::SuperLU_T ) fn.append("SLU");
	else if (linsol==LinAlg::SuperLUd_T) fn.append("SLUd");

	// Check
	LinAlg::Matrix<double> nodeinfo0;
	LinAlg::Matrix<double> nodeinfo1;
	Array<String> labels; 

	Elems[0]->OutNodes(nodeinfo0, labels);
	Elems[1]->OutNodes(nodeinfo1, labels);

	std::cout << labels    << std::endl;
	std::cout << nodeinfo0 << std::endl;
	std::cout << labels    << std::endl;
	std::cout << nodeinfo1 << std::endl;

	size_t sx_idx;
	size_t sy_idx;
	size_t sxy_idx;
	for (size_t i=0; i<labels.Size(); i++)
	{
		if (labels[i]=="Sx")  sx_idx  = i;
		if (labels[i]=="Sy")  sy_idx  = i;
		if (labels[i]=="Sxy") sxy_idx = i;
	} 

    double errors = 0.0;

	// Element 0
	errors += fabs(nodeinfo0(0, sx_idx) - ( 1.56432140e-01));
	errors += fabs(nodeinfo0(1, sx_idx) - (-3.19109076e-01));
	errors += fabs(nodeinfo0(2, sx_idx) - (-3.00686928e-01));
	errors += fabs(nodeinfo0(3, sx_idx) - (-3.31286428e-01));
	errors += fabs(nodeinfo0(4, sx_idx) - ( 1.44254788e-01));
	errors += fabs(nodeinfo0(5, sx_idx) - ( 1.25832639e-01));

	errors += fabs(nodeinfo0(0, sy_idx) - (-2.05141549e-01));
	errors += fabs(nodeinfo0(1, sy_idx) - (-2.22127394e+00));
	errors += fabs(nodeinfo0(2, sy_idx) - ( 1.15872190e+00));
	errors += fabs(nodeinfo0(3, sy_idx) - (-2.96971274e+00));
	errors += fabs(nodeinfo0(4, sy_idx) - (-9.53580350e-01));
	errors += fabs(nodeinfo0(5, sy_idx) - (-4.33357619e+00));

	errors += fabs(nodeinfo0(0, sxy_idx) - (-1.56432140e-01));
	errors += fabs(nodeinfo0(1, sxy_idx) - (-4.90216486e-02));
	errors += fabs(nodeinfo0(2, sxy_idx) - (-6.74437968e-02));
	errors += fabs(nodeinfo0(3, sxy_idx) - ( 3.31286428e-01));
	errors += fabs(nodeinfo0(4, sxy_idx) - ( 2.23875937e-01));
	errors += fabs(nodeinfo0(5, sxy_idx) - ( 2.42298085e-01));

	// Element 1
	errors += fabs(nodeinfo1(0, sx_idx)  - ( 9.95732723e-01));
	errors += fabs(nodeinfo1(1, sx_idx)  - ( 1.39446295e+00));
	errors += fabs(nodeinfo1(2, sx_idx)  - ( 2.23875937e-01));
	errors += fabs(nodeinfo1(3, sx_idx)  - (-8.20878435e-01));
	errors += fabs(nodeinfo1(4, sx_idx)  - (-1.21960866e+00));
	errors += fabs(nodeinfo1(5, sx_idx)  - (-4.90216486e-02));
                                      
	errors += fabs(nodeinfo1(0, sy_idx)  - (-1.25426728e+00));
	errors += fabs(nodeinfo1(1, sy_idx)  - (-2.39612823e+00));
	errors += fabs(nodeinfo1(2, sy_idx)  - ( 1.68328476e+00));
	errors += fabs(nodeinfo1(3, sy_idx)  - (-1.57087843e+00));
	errors += fabs(nodeinfo1(4, sy_idx)  - (-4.29017485e-01));
	errors += fabs(nodeinfo1(5, sy_idx)  - (-4.50843047e+00));
                                      
	errors += fabs(nodeinfo1(0, sxy_idx) - (-9.95732723e-01));
	errors += fabs(nodeinfo1(1, sxy_idx) - ( 1.25832639e-01));
	errors += fabs(nodeinfo1(2, sxy_idx) - ( 1.29641965e+00));
	errors += fabs(nodeinfo1(3, sxy_idx) - ( 8.20878435e-01));
	errors += fabs(nodeinfo1(4, sxy_idx) - (-3.00686928e-01));
	errors += fabs(nodeinfo1(5, sxy_idx) - (-1.47127394e+00));

	if (fabs(errors)>1.0e-7) cout << "[1;31mErrors(" << fn << ") = " << errors << "[0m\n" << endl;
	else                     cout << "[1;32mErrors(" << fn << ") = " << errors << "[0m\n" << endl;

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
