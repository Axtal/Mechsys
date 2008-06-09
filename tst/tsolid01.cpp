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
#include "fem/elems/hex8equilib.h"
#include "fem/solvers/autome.h"
#include "fem/solvers/forwardeuler.h"
#include "models/equilibs/linelastic.h"
#include "util/exception.h"

using FEM::Nodes;
using FEM::Elems;
using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	/*  Unit cube
	 
	      z
	      |__y      4________________7
	   x,'        ,'|              ,'|
	            ,'               ,'  |
	          ,'    |          ,'    |
	        ,'      .        ,'      |
	      5'_______________6'        |
	      |                |         |
	      |         |      |         |
	      |         0 -  - | -  -  - 3
	      |       ,        |       ,' 
	      |     ,          |     ,'   
	      |   ,            |   ,'     
	      | ,              | ,'       
	      1________________2'         
	*/

	// Input
	LinAlg::LinSol_T linsol;
	bool input_ok = false;
	if (argc==2) { if (atoi(argv[1])>=1 && atoi(argv[1])<=3) input_ok = true; }
	if (input_ok) linsol = static_cast<LinAlg::LinSol_T>(atoi(argv[1]));
	else throw new Message(_("Please, call this program as in:\n\t\t %s LinSol\n  where:\n   LinSol:\n \t1 => LAPACK_T  : DENSE\n \t2 => UMFPACK_T : SPARSE\n \t3 => SuperLU_T : SPARSE\n"),argv[0]);

	// 0) Geometry type
	FEM::GeometryType = 3; // 3D

	// 1) Nodes
	FEM::AddNode(0.0, 0.0, 0.0);
	FEM::AddNode(1.0, 0.0, 0.0);
	FEM::AddNode(1.0, 1.0, 0.0);
	FEM::AddNode(0.0, 1.0, 0.0);
	FEM::AddNode(0.0, 0.0, 1.0);
	FEM::AddNode(1.0, 0.0, 1.0);
	FEM::AddNode(1.0, 1.0, 1.0);
	FEM::AddNode(0.0, 1.0, 1.0);

	// 2) Elements
	FEM::AddElem("Hex8Equilib", /*IsActive*/true); // 0

	// 3) Set connectivity
	Elems[0]->SetNode(0, 0);
	Elems[0]->SetNode(1, 1);
	Elems[0]->SetNode(2, 2);
	Elems[0]->SetNode(3, 3);
	Elems[0]->SetNode(4, 4);
	Elems[0]->SetNode(5, 5);
	Elems[0]->SetNode(6, 6);
	Elems[0]->SetNode(7, 7);

	// 4) Boundary conditions (must be after connectivity)
	Nodes[0]->Bry("Dux" ,0.0)->Bry("Duy" ,0.0)->Bry("Duz" ,0.0);
	Nodes[1]->Bry("Duz" ,0.0)->Bry("Duy" ,0.0);
	Nodes[2]->Bry("Duz" ,0.0);
	Nodes[3]->Bry("Duz" ,0.0)->Bry("Dux" ,0.0);
	Nodes[4]->Bry("Dux" ,0.0)->Bry("Duy" ,0.0);
	Nodes[5]->Bry("Duy" ,0.0);
	Nodes[7]->Bry("Dux" ,0.0);
	Nodes[4]->Bry("Dfz", -10.0);
	Nodes[5]->Bry("Dfz", -15.0);
	Nodes[6]->Bry("Dfz", -20.0);
	Nodes[7]->Bry("Dfz", -25.0);

	// 5) Parameters and initial values
	Elems[0]->ReAllocateModel("LinElastic", "E=2000.0 nu=0.2", "Sx=10.0 Sy=10.0 Sz=10");

	// Stiffness
	Array<size_t>          map;
	Array<bool>            pre;
	LinAlg::Matrix<double> ke;
	Elems[0]->Order1Matrix(0,ke);
	cout << ke << endl;

	// 6) Solve
	//FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	FEM::Solver * sol = FEM::AllocSolver("AutoME");
	sol->Solve(linsol, /*iStage*/0, /*NumDiv*/1, /*DeltaTime*/0);

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
