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
#include "fem/elems/quad4heat.h"
#include "models/heats/linheat.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "util/exception.h"
#include "util/numstreams.h"

using FEM::Nodes;
using FEM::Elems;
using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	/*       
	        
	
	          (T=1.0)  3 @-------------------@ 2  (T=1.0)
	                     |                   |      
	                     |                   |      
	              L/2    |                   |
	                     |                   |
	          (T=0.0)  0 @-------------------@ 1  (T=0.0)
	                              L            
	 */

	// Input
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);
	else cout << "[1;32mYou may call this program as in:\t " << argv[0] << " LinSol\n  where LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n [0m[1;34m Now using LA (LAPACK)\n[0m" << endl;

	// 0) Problem dimension
	FEM::Dim = 2; // 2D

	// 1) Nodes
	FEM::AddNode(0.0, 0.0); // 0
	FEM::AddNode(1.0, 0.0); // 1
	FEM::AddNode(1.0, 0.5); // 2
	FEM::AddNode(0.0, 0.5); // 3

	// 2) Elements
	FEM::AddElem("Quad4Heat", /*IsActive*/true);

	// 3) Set connectivity
	Elems[0]->SetNode(0,0)->SetNode(1,1)->SetNode(2,2)->SetNode(3,3);

	// 4) Boundary conditions (must be after connectivity)
	Nodes[0]->Bry("T",0.0);
	Nodes[1]->Bry("T",0.0);
	Nodes[2]->Bry("T",1.0);
	Nodes[3]->Bry("T",1.0);

	// 5) Parameters and initial values
	Elems[0]->SetModel("LinHeat", "k=1.0", "");

	// 6) Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol -> SetLinSol(linsol.GetSTL().c_str()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();

	// Stiffness
	Array<size_t>          map;
	Array<bool>            pre;
	LinAlg::Matrix<double> Ke0;  Ke0.SetNS(Util::_6_3);
	Elems[0]->Order1Matrix(0,Ke0);
	cout << "Ke0=\n" << Ke0 << endl;

	cout << "GFE_Resid = " << FEM::GFE_Resid << endl;

	FEM::WriteVTK("out.vtk");
	

	// Check

	// Check
    //double errors = 0.0;
	//for (int i=0; i<4; ++i)
	//for (int j=0; j<4; ++j)
	//	errors += fabs(Ke0(i,j)-Ke_correct(i,j));

	//if (fabs(errors)>1.0e-10) cout << "[1;31mErrors(" << linsol << ") = " << errors << "[0m\n" << endl;
	//else                      cout << "[1;32mErrors(" << linsol << ") = " << errors << "[0m\n" << endl;

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
