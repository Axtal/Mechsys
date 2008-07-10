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
#include "fem/elems/tri6pstrain.h"
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "util/exception.h"

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
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);
	else cout << "[1;32mYou may call this program as in:\t " << argv[0] << " LinSol\n  where LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n [0m[1;34m Now using LA (LAPACK)\n[0m" << endl;

	// 0) Problem dimension
	FEM::Geom g(2); // 2D

	// 1) Nodes
	g.SetNNodes (9);
	g.SetNode   (0, 0.0, 0.0);
	g.SetNode   (1, 1.0, 0.0);
	g.SetNode   (2, 0.0, 1.0);
	g.SetNode   (3, 0.5, 0.0);
	g.SetNode   (4, 0.5, 0.5);
	g.SetNode   (5, 0.0, 0.5);
	g.SetNode   (6, 1.0, 1.0);
	g.SetNode   (7, 0.5, 1.0);
	g.SetNode   (8, 1.0, 0.5);

	// 2) Elements
	g.SetNElems (2);
	g.SetElem   (0, "Tri6PStrain", /*IsActive*/true);
	g.SetElem   (1, "Tri6PStrain", /*IsActive*/true);

	// 3) Set connectivity
	g.Ele(0)->SetNode(0, g.Nod(0))
	        ->SetNode(1, g.Nod(1))
			->SetNode(2, g.Nod(2))
			->SetNode(3, g.Nod(3))
			->SetNode(4, g.Nod(4))
			->SetNode(5, g.Nod(5));
	g.Ele(1)->SetNode(0, g.Nod(6))
	        ->SetNode(1, g.Nod(2))
	        ->SetNode(2, g.Nod(1))
	        ->SetNode(3, g.Nod(7))
	        ->SetNode(4, g.Nod(4))
	        ->SetNode(5, g.Nod(8));

	// 4) Boundary conditions (must be after connectivity)
	g.Nod(0)->Bry("ux",0.0)->Bry("uy",0.0);
	g.Nod(1)->Bry("uy",0.0);
	g.Nod(3)->Bry("uy",0.0);
	g.Nod(2)->Bry("fy",1.0);
	g.Nod(7)->Bry("fy",1.0);
	g.Nod(6)->Bry("fy",1.0);

	// 5) Parameters and initial values
	g.Ele(0)->SetModel("LinElastic", "E=10000.0 nu=0.25", "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0");
	g.Ele(1)->SetModel("LinElastic", "E=10000.0 nu=0.25", "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0");

	// Stiffness
	Array<size_t>          map;
	Array<bool>            pre;
	LinAlg::Matrix<double> Ke0;
	LinAlg::Matrix<double> Ke1;
	g.Ele(0)->Order1Matrix(0,Ke0);
	g.Ele(1)->Order1Matrix(0,Ke1);
	cout << "Ke0=\n" << Ke0 << endl;
	cout << "Ke1=\n" << Ke1 << endl;

	// 6) Solve
	//FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	FEM::Solver * sol = FEM::AllocSolver("AutoME");
	sol -> SetGeom(&g) -> SetLinSol(linsol.GetSTL().c_str()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();

	// Output
	LinAlg::Matrix<double> values0;
	LinAlg::Matrix<double> values1;
	Array<String>          labels0;
	Array<String>          labels1;
	g.Ele(0)->OutNodes (values0, labels0);
	g.Ele(1)->OutNodes (values1, labels1);
	std::cout << labels0;
	std::cout << values0 << std::endl;
	std::cout << labels1;
	std::cout << values1 << std::endl;

	// Check
    double errors = 0.0;

	// Element 0
	errors += fabs(g.Ele(0)->Val(0, "Sx") - ( 1.56432140e-01));
	errors += fabs(g.Ele(0)->Val(1, "Sx") - (-3.00686928e-01));
	errors += fabs(g.Ele(0)->Val(2, "Sx") - ( 1.44254788e-01));
	errors += fabs(g.Ele(0)->Val(3, "Sx") - (-3.19109076e-01));
	errors += fabs(g.Ele(0)->Val(4, "Sx") - (-3.31286428e-01));
	errors += fabs(g.Ele(0)->Val(5, "Sx") - ( 1.25832639e-01));

	errors += fabs(g.Ele(0)->Val(0, "Sy") - (-2.05141549e-01));
	errors += fabs(g.Ele(0)->Val(1, "Sy") - ( 1.15872190e+00));
	errors += fabs(g.Ele(0)->Val(2, "Sy") - (-9.53580350e-01));
	errors += fabs(g.Ele(0)->Val(3, "Sy") - (-2.22127394e+00));
	errors += fabs(g.Ele(0)->Val(4, "Sy") - (-2.96971274e+00));
	errors += fabs(g.Ele(0)->Val(5, "Sy") - (-4.33357619e+00));

	errors += fabs(g.Ele(0)->Val(0, "Sxy") - (-1.56432140e-01));
	errors += fabs(g.Ele(0)->Val(1, "Sxy") - (-6.74437968e-02));
	errors += fabs(g.Ele(0)->Val(2, "Sxy") - ( 2.23875937e-01));
	errors += fabs(g.Ele(0)->Val(3, "Sxy") - (-4.90216486e-02));
	errors += fabs(g.Ele(0)->Val(4, "Sxy") - ( 3.31286428e-01));
	errors += fabs(g.Ele(0)->Val(5, "Sxy") - ( 2.42298085e-01));

	// Element 1
	errors += fabs(g.Ele(1)->Val(0, "Sx")  - ( 9.95732723e-01));
	errors += fabs(g.Ele(1)->Val(1, "Sx")  - ( 2.23875937e-01));
	errors += fabs(g.Ele(1)->Val(2, "Sx")  - (-1.21960866e+00));
	errors += fabs(g.Ele(1)->Val(3, "Sx")  - ( 1.39446295e+00));
	errors += fabs(g.Ele(1)->Val(4, "Sx")  - (-8.20878435e-01));
	errors += fabs(g.Ele(1)->Val(5, "Sx")  - (-4.90216486e-02));

	errors += fabs(g.Ele(1)->Val(0, "Sy")  - (-1.25426728e+00));
	errors += fabs(g.Ele(1)->Val(1, "Sy")  - ( 1.68328476e+00));
	errors += fabs(g.Ele(1)->Val(2, "Sy")  - (-4.29017485e-01));
	errors += fabs(g.Ele(1)->Val(3, "Sy")  - (-2.39612823e+00));
	errors += fabs(g.Ele(1)->Val(4, "Sy")  - (-1.57087843e+00));
	errors += fabs(g.Ele(1)->Val(5, "Sy")  - (-4.50843047e+00));

	errors += fabs(g.Ele(1)->Val(0, "Sxy") - (-9.95732723e-01));
	errors += fabs(g.Ele(1)->Val(1, "Sxy") - ( 1.29641965e+00));
	errors += fabs(g.Ele(1)->Val(2, "Sxy") - (-3.00686928e-01));
	errors += fabs(g.Ele(1)->Val(3, "Sxy") - ( 1.25832639e-01));
	errors += fabs(g.Ele(1)->Val(4, "Sxy") - ( 8.20878435e-01));
	errors += fabs(g.Ele(1)->Val(5, "Sxy") - (-1.47127394e+00));

	if (fabs(errors)>1.0e-7) cout << "[1;31mErrors(" << linsol << ") = " << errors << "[0m\n" << endl;
	else                     cout << "[1;32mErrors(" << linsol << ") = " << errors << "[0m\n" << endl;

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
