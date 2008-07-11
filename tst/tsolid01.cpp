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
#include "fem/geometry.h"
#include "fem/functions.h"
#include "fem/elems/hex8equilib.h"
#include "fem/solvers/autome.h"
#include "fem/solvers/forwardeuler.h"
#include "models/equilibs/linelastic.h"
#include "util/exception.h"

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
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);
	else cout << "[1;32mYou may call this program as in:\t " << argv[0] << " LinSol\n  where LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n [0m[1;34m Now using LA (LAPACK)\n[0m" << endl;

	// 0) Problem dimension
	FEM::Geom g(3); // 3D

	// 1) Nodes
	g.SetNNodes (8);
	g.SetNode   (0, 0.0, 0.0, 0.0);
	g.SetNode   (1, 1.0, 0.0, 0.0);
	g.SetNode   (2, 1.0, 1.0, 0.0);
	g.SetNode   (3, 0.0, 1.0, 0.0);
	g.SetNode   (4, 0.0, 0.0, 1.0);
	g.SetNode   (5, 1.0, 0.0, 1.0);
	g.SetNode   (6, 1.0, 1.0, 1.0);
	g.SetNode   (7, 0.0, 1.0, 1.0);

	// 2) Elements
	g.SetNElems (1);
	g.SetElem   (0, "Hex8Equilib", /*IsActive*/true);

	// 3) Set connectivity
	g.Ele(0)->Connect(0, g.Nod(0))
	        ->Connect(1, g.Nod(1))
	        ->Connect(2, g.Nod(2))
	        ->Connect(3, g.Nod(3))
	        ->Connect(4, g.Nod(4))
	        ->Connect(5, g.Nod(5))
	        ->Connect(6, g.Nod(6))
	        ->Connect(7, g.Nod(7));

	// 4) Boundary conditions (must be after connectivity)
	g.Nod(0)->Bry("ux",   0.0)->Bry("uy" ,0.0)->Bry("uz" ,0.0);
	g.Nod(1)->Bry("uz",   0.0)->Bry("uy" ,0.0);
	g.Nod(2)->Bry("uz",   0.0);
	g.Nod(3)->Bry("uz",   0.0)->Bry("ux" ,0.0);
	g.Nod(4)->Bry("ux",   0.0)->Bry("uy" ,0.0);
	g.Nod(5)->Bry("uy",   0.0);
	g.Nod(7)->Bry("ux",   0.0);
	g.Nod(4)->Bry("fz", -10.0);
	g.Nod(5)->Bry("fz", -15.0);
	g.Nod(6)->Bry("fz", -20.0);
	g.Nod(7)->Bry("fz", -25.0);

	// 5) Parameters and initial values
	g.Ele(0)->SetModel("LinElastic", "E=2000.0 nu=0.2", "Sx=10.0 Sy=10.0 Sz=10");

	// Stiffness
	Array<size_t>          map;
	Array<bool>            pre;
	LinAlg::Matrix<double> ke;
	g.Ele(0)->Order1Matrix(0,ke);
	cout << ke << endl;

	// 6) Solve
	//FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	FEM::Solver * sol = FEM::AllocSolver("AutoME");
	sol -> SetGeom(&g) -> SetLinSol(linsol.GetSTL().c_str()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
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
