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
#include "fem/elems/quad4pstress.h"
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "util/exception.h"
#include "util/numstreams.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	/*       F1                 F3    F1=F2=F3 = 1.0
	          ^                   ^
	          |                   |
	
	        3 @-------------------@ 2    Prof. Carlos Felippa
	          |                   |      IFEM.Ch23.pdf
	          |                   |      Pag. 23-7
	   L/2    |                   |
	          |                   |
	        0 @-------------------@ 1
	         /_\                 /_\
	         ///       L         o o
	 */

	// Input
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);
	else cout << "[1;32mYou may call this program as in:\t " << argv[0] << " LinSol\n  where LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n [0m[1;34m Now using LA (LAPACK)\n[0m" << endl;

	// 0) Geometry
	FEM::Geom g(2); // 2D

	// 1) Nodes
	g.SetNNodes (4);
	g.SetNode   (0, 0.0, 0.0);
	g.SetNode   (1, 1.0, 0.0);
	g.SetNode   (2, 1.0, 0.5);
	g.SetNode   (3, 0.0, 0.5);

	// 2) Elements
	g.SetNElems (1);
	g.SetElem   (0, "Quad4PStress", /*IsActive*/true);

	// 3) Set connectivity
	g.Ele(0)->SetNode(0, g.Nod(0))
	        ->SetNode(1, g.Nod(1))
	        ->SetNode(2, g.Nod(2))
	        ->SetNode(3, g.Nod(3));

	// 4) Boundary conditions (must be after connectivity)
	g.Nod(0)->Bry("ux",0.0)->Bry("uy",0.0);
	g.Nod(1)->Bry("uy",0.0);

	// 5) Parameters and initial values
	g.Ele(0)->SetModel("LinElastic", "E=96.0 nu=0.333333333333333333333333", "Sx=0.0 Sy=0.0 Sxy=0.0");

	// 6) Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol -> SetGeom(&g) -> SetLinSol(linsol.GetSTL().c_str()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();

	// Stiffness
	Array<size_t>          map;
	Array<bool>            pre;
	LinAlg::Matrix<double> Ke0;  Ke0.SetNS(Util::_6_3);
	g.Ele(0)->Order1Matrix(0,Ke0);
	cout << "Ke0=\n" << Ke0 << endl;

	// Check
	LinAlg::Matrix<double> Ke_correct; Ke_correct.Resize(8,8);
	Ke_correct =  42.0,  18.0,  -6.0,   0.0, -21.0, -18.0, -15.0,   0.0,
	              18.0,  78.0,   0.0,  30.0, -18.0, -39.0,   0.0, -69.0,
	              -6.0,   0.0,  42.0, -18.0, -15.0,   0.0, -21.0,  18.0,
	               0.0,  30.0, -18.0,  78.0,   0.0, -69.0,  18.0, -39.0,
	             -21.0, -18.0, -15.0,   0.0,  42.0,  18.0,  -6.0,   0.0,
	             -18.0, -39.0,   0.0, -69.0,  18.0,  78.0,   0.0,  30.0,
	             -15.0,   0.0, -21.0,  18.0,  -6.0,   0.0,  42.0, -18.0,
	               0.0, -69.0,  18.0, -39.0,   0.0,  30.0, -18.0,  78.0;

	// Check
    double errors = 0.0;
	for (int i=0; i<8; ++i)
	for (int j=0; j<8; ++j)
		errors += fabs(Ke0(i,j)-Ke_correct(i,j));

	if (fabs(errors)>1.0e-10) cout << "[1;31mErrors(" << linsol << ") = " << errors << "[0m\n" << endl;
	else                      cout << "[1;32mErrors(" << linsol << ") = " << errors << "[0m\n" << endl;

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
