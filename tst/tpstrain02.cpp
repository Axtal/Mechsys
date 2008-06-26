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
#include "fem/elems/quad4pstrain.h"
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "util/exception.h"
#include "linalg/matrix.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using FEM::Nodes;
using FEM::Elems;
using LinAlg::Matrix;

int main(int argc, char **argv) try
{
	double H     = 2.0;   // height
	double L     = 2.0;   // length
	double E     = 207.0; // Young
	double nu    = 0.3;   // Poisson
	double q     = 1.0;   // Load
	size_t ndivx = 2;     // divisions along x (must be even)
	size_t ndivy = 2;     // divisions along y

	/*        | | | | | | | | | | | | | | | | |  q
	          V V V V V V V V V V V V V V V V V 
	  -+-     @-------@-------@-------@-------@
	   |      |20     |21     |22     |23     |24
	   |      |   e12 |   e13 |   e14 |   e15 |
	   |      |15     |16     |17     |18     |19
	   |      @-------@-------@-------@-------@
	   |      |       |       |       |       |
	   |      |   e8  |   e9  |   e10 |   e11 |
	          |10     |11     |12     |13     |14
	   H      @-------@-------@-------@-------@
	          |       |       |       |       |
	   |      |   e4  |   e5  |   e6  |   e7  |
	   |      |5      |6      |7      |8      |9
	   |      @-------@-------@-------@-------@
	   |      |       |       |       |       |
	   |      |   e0  |   e1  |   e2  |   e3  |
	   |      |0      |1      |2      |3      |4
	  -+-     @-------@-------@-------@-------@
	         /_\     /_\     /_\     /_\     /_\
	         o o     o o     ///     o o     o o
	
	          |-------------- L --------------|
	 */

	// Input
	String linsol("LA");
	if (argc==2)   linsol.Printf("%s",argv[1]);
	if (argc==3) { linsol.Printf("%s",argv[1]); ndivx = ndivy = atoi(argv[2]); }
	else cout << "[1;32mYou may call this program as in:\t " << argv[0] << " LinSol nDiv\n  where LinSol:\n \tLA   => LAPACK_T  : DENSE\n \tUM   => UMFPACK_T : SPARSE\n \tSLU  => SuperLU_T : SPARSE\n \tnDiv => Number of division along x and y (must be even)\n [0m[1;34m Now using LA (LAPACK)\n[0m" << endl;

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Block coordinates
	Matrix<double> c(2,8);
	c = 0.,  L, L, 0.,    L/2.,    L, L/2.,   0.,
		0., 0., H,  H,      0., H/2.,    H, H/2.;

	// Block weights
	Array<double> wx(ndivx);  wx = 1.0;
	Array<double> wy(ndivy);  wy = 1.0;

	// Block tags
	Vector<int> e_tags(4);
	e_tags = 0, 0, -10, -20;

	// Blocks
	Mesh::Block b;
	b.Set      (&c, &wx, &wy);
	b.SetETags (&e_tags);
	Array<Mesh::Block*> blocks;  blocks.Push(&b);

	// Generate
	cout << "\nMesh Generation: --------------------------------------------------------------" << endl;
	Mesh::Structured ms;
	clock_t start = std::clock(); // Initial time
	size_t  ne    = ms.Generate (blocks);
	clock_t total = std::clock() - start; // Time elapsed
	cout << "[1;33m"<<ne<<" elements[0m. Time elapsed = [1;31m" << static_cast<double>(total)/CLOCKS_PER_SEC << "[0m [1;32mseconds[0m" << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// 0) Problem dimension
	FEM::Dim = 2; // 2D

	// 1) Nodes
	for (size_t i=0; i<ms.Verts().Size(); ++i)
		FEM::AddNode (ms.Verts()[i]->C(0), ms.Verts()[i]->C(1));

	for (size_t i=0; i<ms.Elems().Size(); ++i)
	{
		// 2) Elements
		FEM::Element * e = FEM::AddElem("Quad4PStrain");

		// 3) Connectivity
		Mesh::Elem * me = ms.Elems()[i];
		e->SetNode(0, me->V[0]->MyID);
		e->SetNode(1, me->V[1]->MyID);
		e->SetNode(2, me->V[2]->MyID);
		e->SetNode(3, me->V[3]->MyID);

		// 5) Parameters and initial values
		String prms; prms.Printf("E=%f  nu=%f",E,nu);
		e->SetModel("LinElastic", prms.GetSTL().c_str(), "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0");
	}

	// 4) Boundary conditions
	// Nodes
	Array<double>       X, Y, Z;
	Array<char const *> node_vars;
	Array<double>       node_vals;
	X.Push(L/2.0);  Y.Push(0.0);  Z.Push(0.0);  node_vars.Push("ux");  node_vals.Push(0.0);
	FEM::SetNodeBrys (&ms, &X, &Y, &Z, &node_vars, &node_vals);
	// Faces
	Array<int>          face_tags;
	Array<char const *> face_vars;
	Array<double>       face_vals;
	face_tags.Push(-10);  face_vars.Push("uy");  face_vals.Push(0.0);
	face_tags.Push(-20);  face_vars.Push("fy");  face_vals.Push( -q);
	FEM::SetFaceBrys (&ms, &face_tags, &face_vars, &face_vals);

	// 6) Solve
	cout << "\nSolution: ---------------------------------------------------------------------" << endl;
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	//FEM::Solver * sol = FEM::AllocSolver("AutoME");
	start = std::clock(); // Initial time
	sol -> SetLinSol(linsol.GetSTL().c_str()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();
	total = std::clock() - start; // Time elapsed
	cout << "GFE_Resid = "<<FEM::GFE_Resid<<". Time elapsed = [1;31m"<<static_cast<double>(total)/CLOCKS_PER_SEC<<"[0m [1;32mseconds[0m"<<std::endl;

	// Check
    double errors = 0.0;

	double Sy = q;
	double Ex = -nu*(1.0+nu)*Sy/E;
	double Ey =  (1.0-nu*nu)*Sy/E;
	double Sz = (E/(1.0+nu))*(nu/(1.0-2.0*nu))*(Ex+Ey);

	// Stress and strains
	for (size_t i=0; i<Elems.Size(); ++i)
	{
		errors += fabs(Elems[i]->Val("Ex" ) - (Ex));
		errors += fabs(Elems[i]->Val("Ey" ) - (Ey));
		errors += fabs(Elems[i]->Val("Exy") - (0.0));
		errors += fabs(Elems[i]->Val("Sx" ) - (0.0));
		errors += fabs(Elems[i]->Val("Sy" ) - (Sy ));
		errors += fabs(Elems[i]->Val("Sz" ) - (Sz ));
		errors += fabs(Elems[i]->Val("Sxy") - (0.0));
	}

	// Displacements
	for (size_t i=0; i<Nodes.Size(); ++i)
	{
		if (fabs(Nodes[i]->Y())<1.0e-5) // bottom nodes
		{
			errors += fabs(Nodes[i]->Val("uy") - (0.0));
			if (fabs(Nodes[i]->X()-(L/2.0))<1.0e-5) // central node
				errors += fabs(Nodes[i]->Val("ux") - (0.0));
		}
		else if (fabs(Nodes[i]->Y()-(H/2.0))<1.0e-5) // mid nodes
			errors += fabs(Nodes[i]->Val("uy") - (-0.5*H*Ey));
		else if (fabs(Nodes[i]->Y()-(H))<1.0e-5) // top nodes
			errors += fabs(Nodes[i]->Val("uy") - (-H*Ey));
		if (fabs(Nodes[i]->X())<1.0e-5) // left nodes
			errors += fabs(Nodes[i]->Val("ux") - (0.5*L*Ex));
		if (fabs(Nodes[i]->X()-L)<1.0e-5) // right nodes
			errors += fabs(Nodes[i]->Val("ux") - (-0.5*L*Ex));
	}

	if (fabs(errors)>1.0e-10) cout << "[1;31m\nErrors(" << linsol << ") = " << errors << "[0m\n" << endl;
	else                      cout << "[1;32m\nErrors(" << linsol << ") = " << errors << "[0m\n" << endl;

	// Write file
	cout << "Write output file: ------------------------------------------------------------" << endl;
	start = std::clock(); // Initial time
	FEM::WriteVTUEquilib ("tpstrain02.vtu");
	total = std::clock() - start; // Time elapsed
	cout << "[1;34mFile <tpstrain02.vtu> saved.[0m" << endl;
	cout << "Time elapsed = [1;31m"<<static_cast<double>(total)/CLOCKS_PER_SEC<<"[0m [1;32mseconds[0m"<<std::endl;
	cout << endl;

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
