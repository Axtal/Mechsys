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
	if (argc!=2) throw new Message(_("Please, call this program as in:\n\t\t %s LinSol\n  where:\n   LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n"),argv[0]);

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
	Nodes[0]->Bry("ux",   0.0)->Bry("uy" ,0.0)->Bry("uz" ,0.0);
	Nodes[1]->Bry("uz",   0.0)->Bry("uy" ,0.0);
	Nodes[2]->Bry("uz",   0.0);
	Nodes[3]->Bry("uz",   0.0)->Bry("ux" ,0.0);
	Nodes[4]->Bry("ux",   0.0)->Bry("uy" ,0.0);
	Nodes[5]->Bry("uy",   0.0);
	Nodes[7]->Bry("ux",   0.0);
	Nodes[4]->Bry("fz", -10.0);
	Nodes[5]->Bry("fz", -15.0);
	Nodes[6]->Bry("fz", -20.0);
	Nodes[7]->Bry("fz", -25.0);

	// 5) Parameters and initial values
	Elems[0]->SetModel("LinElastic", "E=2000.0 nu=0.2", "Sx=10.0 Sy=10.0 Sz=10");

	// Stiffness
	Array<size_t>          map;
	Array<bool>            pre;
	LinAlg::Matrix<double> ke;
	Elems[0]->Order1Matrix(0,ke);
	cout << ke << endl;

	// 6) Solve
	//FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	FEM::Solver * sol = FEM::AllocSolver("AutoME");
	sol -> SetLinSol(argv[1]) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();

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
