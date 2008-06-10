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
#include "fem/elems/elasticrod.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "util/exception.h"

using FEM::Nodes;
using FEM::Elems;
using std::cout;
using std::endl;

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
	LinAlg::LinSol_T linsol;
	bool input_ok = false;
	if (argc==2) { if (atoi(argv[1])>=1 && atoi(argv[1])<=3) input_ok = true; }
	if (input_ok) linsol = static_cast<LinAlg::LinSol_T>(atoi(argv[1]));
	else throw new Message(_("Please, call this program as in:\n\t\t %s LinSol\n  where:\n   LinSol:\n \t1 => LAPACK_T  : DENSE\n \t2 => UMFPACK_T : SPARSE\n \t3 => SuperLU_T : SPARSE\n"),argv[0]);

	// 0) Geometry type
	FEM::GeometryType = 2; // 2D

	// 1) Nodes
	FEM::AddNode( 0.0,  0.0); // 0
	FEM::AddNode(10.0,  0.0); // 1
	FEM::AddNode(10.0, 10.0); // 2

	// 2) Elements
	FEM::AddElem("ElasticRod", /*IsActive*/true); // 0
	FEM::AddElem("ElasticRod", /*IsActive*/true); // 1
	FEM::AddElem("ElasticRod", /*IsActive*/true); // 2

	// 3) Set connectivity
	Elems[0]->SetNode(0, 0)->SetNode(1, 1);
	Elems[1]->SetNode(0, 1)->SetNode(1, 2);
	Elems[2]->SetNode(0, 0)->SetNode(1, 2);

	// 4) Boundary conditions (must be after set connectivity)
	Nodes[0]->Bry("Dux", 0.0)->Bry("Duy", -0.5)->Bry("Duz", 0.0); // Essential
	Nodes[1]->                 Bry("Duy",  0.4)->Bry("Duz", 0.0); // Essential
	Nodes[2]                                   ->Bry("Duz", 0.0); // Essential
	Nodes[2]->Bry("Dfx", 2.0)->Bry("Dfy",  1.0);                  // Natural

	// 5) Parameters and initial values
	Elems[0]->SetModel("", "E=100.0", "N=0.0  A=1.0");
	Elems[1]->SetModel("", "E= 50.0", "N=0.0  A=1.0");
	Elems[2]->SetModel("", "E=200.0", "N=0.0  A=1.414213562373095");

	// Stiffness
	Array<size_t>          map;
	Array<bool>            pre;
	LinAlg::Matrix<double> Ke0;
	LinAlg::Matrix<double> Ke1;
	LinAlg::Matrix<double> Ke2;
	Elems[0]->Order1Matrix(0,Ke0);
	Elems[1]->Order1Matrix(0,Ke1);
	Elems[2]->Order1Matrix(0,Ke2);
	cout << "Ke0=\n" << Ke0 << endl;
	cout << "Ke1=\n" << Ke1 << endl;
	cout << "Ke2=\n" << Ke2 << endl;

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
	double errors = 0;

	errors += fabs(Nodes[0]->DOFVar("Dux").EssentialVal - ( 0.0));
	errors += fabs(Nodes[0]->DOFVar("Duy").EssentialVal - (-0.5));
	errors += fabs(Nodes[0]->DOFVar("Duz").EssentialVal - ( 0.0));
	errors += fabs(Nodes[1]->DOFVar("Dux").EssentialVal - ( 0.0));
	errors += fabs(Nodes[1]->DOFVar("Duy").EssentialVal - ( 0.4));
	errors += fabs(Nodes[1]->DOFVar("Duz").EssentialVal - ( 0.0));
	errors += fabs(Nodes[2]->DOFVar("Dux").EssentialVal - (-0.5));
	errors += fabs(Nodes[2]->DOFVar("Duy").EssentialVal - ( 0.2));
	errors += fabs(Nodes[2]->DOFVar("Duz").EssentialVal - ( 0.0));

	errors += fabs(Nodes[0]->DOFVar("Dfx").NaturalVal - (-2.0));
	errors += fabs(Nodes[0]->DOFVar("Dfy").NaturalVal - (-2.0));
	errors += fabs(Nodes[0]->DOFVar("Dfz").NaturalVal - ( 0.0));
	errors += fabs(Nodes[1]->DOFVar("Dfx").NaturalVal - ( 0.0));
	errors += fabs(Nodes[1]->DOFVar("Dfy").NaturalVal - ( 1.0));
	errors += fabs(Nodes[1]->DOFVar("Dfz").NaturalVal - ( 0.0));
	errors += fabs(Nodes[2]->DOFVar("Dfx").NaturalVal - ( 2.0));
	errors += fabs(Nodes[2]->DOFVar("Dfy").NaturalVal - ( 1.0));
	errors += fabs(Nodes[2]->DOFVar("Dfz").NaturalVal - ( 0.0));

	if (fabs(errors)>1.0e-14) cout << "[1;31mErrors(" << fn << ") = " << errors << "[0m\n" << endl;
	else                      cout << "[1;32mErrors(" << fn << ") = " << errors << "[0m\n" << endl;

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
