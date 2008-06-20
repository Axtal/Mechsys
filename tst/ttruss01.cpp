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
#include "fem/elems/rod.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "models/equilibs/linelastic.h"
#include "util/exception.h"

using FEM::Nodes;
using FEM::Elems;
using std::cout;
using std::endl;
using LinAlg::Matrix;

int main(int argc, char **argv) try
{
	/*  Small truss
	                              1.0 ^
	                                  |
	                                  |2
	                                  o----> 2.0
	                                ,'|
	                              ,'  |
	                   E=200    ,'    |
	                   A=SQ2  ,'      | E=50
	                    [2] ,'        | A=1
	                      ,'          | [1]
	                    ,'            |
	                  ,'              |
	   y            ,'                |
	   |         0,'        [0]       |1
	   |         o--------------------o
	  (z)___x   /_\        E=100     /_\
	           ////        A=1       ___  
	*/

	// Input
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);
	else cout << "[1;32mYou may call this program as in:\t " << argv[0] << " LinSol\n  where LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n [0m[1;34m Now using LA (LAPACK)\n[0m" << endl;

	//////////////////////////////////////////////////////////////////////////////////////////////// 2D
	
	{
		// 0) Problem dimension
		FEM::Dim = 2; // 2D

		// 1) Nodes
		FEM::AddNode( 0.0,  0.0); // 0
		FEM::AddNode(10.0,  0.0); // 1
		FEM::AddNode(10.0, 10.0); // 2

		// 2) Elements
		FEM::AddElem("Rod", /*IsActive*/true); // 0
		FEM::AddElem("Rod", /*IsActive*/true); // 1
		FEM::AddElem("Rod", /*IsActive*/true); // 2

		// 3) Set connectivity
		Elems[0]->SetNode(0, 0)->SetNode(1, 1);
		Elems[1]->SetNode(0, 1)->SetNode(1, 2);
		Elems[2]->SetNode(0, 0)->SetNode(1, 2);

		// 4) Boundary conditions (must be after set connectivity)
		Nodes[0]->Bry("ux", 0.0)->Bry("uy", -0.5); // Essential
		Nodes[1]->                Bry("uy",  0.4); // Essential
		Nodes[2]->Bry("fx", 2.0)->Bry("fy",  1.0); // Natural

		// 5) Parameters and initial value
		Elems[0]->SetModel("LinElastic", "E=100.0 A=1.0"              , "Sx=0.0");
		Elems[1]->SetModel("LinElastic", "E= 50.0 A=1.0"              , "Sx=0.0");
		Elems[2]->SetModel("LinElastic", "E=200.0 A=1.414213562373095", "Sx=0.0");

		// Check
		double errors = 0;

		// Stiffness
		Array<size_t>  map;
		Array<bool>    pre;
		Matrix<double> Ke0;
		Matrix<double> Ke1;
		Matrix<double> Ke2;
		Elems[0]->Order1Matrix(0,Ke0);
		Elems[1]->Order1Matrix(0,Ke1);
		Elems[2]->Order1Matrix(0,Ke2);
		cout << "Ke0=\n" << Ke0 << endl;
		cout << "Ke1=\n" << Ke1 << endl;
		cout << "Ke2=\n" << Ke2 << endl;

		Matrix<double> Ke0_correct;  Ke0_correct.Resize(4,4);
		Ke0_correct =  10.0,   0.0, -10.0,   0.0,
		                0.0,   0.0,   0.0,   0.0,
		              -10.0,   0.0,  10.0,   0.0,
		                0.0,   0.0,   0.0,   0.0;
		Matrix<double> Ke1_correct;  Ke1_correct.Resize(4,4);
		Ke1_correct =   0.0,   0.0,   0.0,   0.0,
		                0.0,   5.0,   0.0,  -5.0,
		                0.0,   0.0,   0.0,   0.0,
		                0.0,  -5.0,   0.0,   5.0;
		Matrix<double> Ke2_correct;  Ke2_correct.Resize(4,4);
		Ke2_correct =  10.0,  10.0, -10.0, -10.0,
		               10.0,  10.0, -10.0, -10.0,
		              -10.0, -10.0,  10.0,  10.0,
		              -10.0, -10.0,  10.0,  10.0;

		for (int i=0; i<4; ++i)
		for (int j=0; j<4; ++j)
		{
			errors += fabs(Ke0(i,j)-Ke0_correct(i,j));
			errors += fabs(Ke1(i,j)-Ke1_correct(i,j));
			errors += fabs(Ke2(i,j)-Ke2_correct(i,j));
		}

		// 6) Solve
		FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
		//FEM::Solver * sol = FEM::AllocSolver("AutoME");
		sol -> SetLinSol(linsol.GetSTL().c_str()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
		sol -> Solve();
		cout << "GFE_Resid = " << FEM::GFE_Resid << endl;

		errors += fabs(Nodes[0]->Val("ux") - ( 0.0));
		errors += fabs(Nodes[0]->Val("uy") - (-0.5));
		errors += fabs(Nodes[1]->Val("ux") - ( 0.0));
		errors += fabs(Nodes[1]->Val("uy") - ( 0.4));
		errors += fabs(Nodes[2]->Val("ux") - (-0.5));
		errors += fabs(Nodes[2]->Val("uy") - ( 0.2));

		errors += fabs(Nodes[0]->Val("fx") - (-2.0));
		errors += fabs(Nodes[0]->Val("fy") - (-2.0));
		errors += fabs(Nodes[1]->Val("fx") - ( 0.0));
		errors += fabs(Nodes[1]->Val("fy") - ( 1.0));
		errors += fabs(Nodes[2]->Val("fx") - ( 2.0));
		errors += fabs(Nodes[2]->Val("fy") - ( 1.0));

		if (fabs(errors)>1.0e-13) cout << "[1;31m2D ==> Errors(" << linsol << ") = " << errors << "[0m\n" << endl;
		else                      cout << "[1;32m2D ==> Errors(" << linsol << ") = " << errors << "[0m\n" << endl;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////// 3D

	{
	}

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
