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
#include "fem/elems/quad4biot.h" // << plane strain
//#include "fem/elems/quad8pstrain.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "fem/output.h"
#include "util/exception.h"
#include "linalg/matrix.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_8s;
using boost::make_tuple;

void CallSolve (int iStage, FEM::Solver * Sol)
{
	double start = std::clock();
	Sol -> Solve();
	double total = std::clock() - start;
	double norm_resid = LinAlg::Norm(Sol->Resid());
	cout << "[1;31mStage # " << iStage << "[0m" << endl;
	cout << "\tTime elapsed (solution) = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";
	cout << "\t[1;35mNorm(Resid=DFext-DFint) = " << norm_resid << "[0m\n";
	cout << "\t[1;32mNumber of DOFs          = " << Sol->nDOF() << "[0m\n";
	// Output: VTU
	//Output o; o.VTU (&g, "tconsolid01.vtu");
	//cout << "[1;34mFile <tconsolid01.vtu> saved.[0m\n\n";
}

int main(int argc, char **argv) try
{
	// Constants
	double W     = 1.0;    // Width
	double H     = 10.0;   // Height
	double E     = 5000.0; // Young
	double nu    = 0.25;   // Poisson
	double gw    = 10.0;   // GammaW
	double k     = 1.0e-5; // Isotropic permeability
	int    ndivy = 10;     // number of divisions along x and y
	bool   is_o2 = false;  // use high order elements?

	// Input
	cout << "Input: " << argv[0] << "  is_o2  ndivy\n";
	if (argc>=2) is_o2 = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) ndivy =  atof(argv[2]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Blocks
	Mesh::Block b;
	b.SetTag    (-1); // tag to be replicated to all generated elements inside this block
	b.SetCoords (false, 4,                // Is3D, NNodes
	             0.0,   W,  W, 0.0,       // x coordinates
	             0.0, 0.0,  H,   H);      // y coordinates
	b.SetNx     (1);                      // x weights and num of divisions along x
	b.SetNy     (ndivy);                  // y weights and num of divisions along y
	b.SetETags  (4,  -10, -20, -30, -40); // edge tags
	Array<Mesh::Block*> blocks;
	blocks.Push (&b);

	// Generate
	Mesh::Structured mesh(/*Is3D*/false);
	if (is_o2) mesh.SetO2();                // Non-linear elements
	clock_t start = std::clock();           // Initial time
	size_t  ne    = mesh.Generate (blocks); // Discretize domain
	clock_t total = std::clock() - start;   // Time elapsed
	if (is_o2) cout << "\nNum of quadrangles (o2) = " << ne << endl;
	else       cout << "\nNumber of quadrangles   = " << ne << endl;
	cout << "Time elapsed (mesh)     = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Geom g(2); // 2D

	// Elements attributes
	String prms; prms.Printf("gw=%f E=%f nu=%f k=%f",gw,E,nu,k);
	FEM::EAtts_T eatts;
	if (is_o2) eatts.Push (make_tuple(-1, "Quad8Biot", "", prms.CStr(), "ZERO"));
	else       eatts.Push (make_tuple(-1, "Quad4Biot", "", prms.CStr(), "ZERO"));

	// Set geometry: nodes, elements, attributes, and boundaries
	FEM::SetNodesElems (&mesh, &eatts, &g);

	// Solver
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol->SetGeom(&g)->SetLinSol("UM");

	// Edges boundaries
	FEM::EBrys_T ebrys;

	// Stage # 1
	sol->SetNumDiv(4)->SetDeltaTime(10.0);
	ebrys.Resize (0);
	ebrys.Push   (make_tuple(-10, "ux",    0.0));
	ebrys.Push   (make_tuple(-20, "ux",    0.0));
	ebrys.Push   (make_tuple(-30, "uy",    0.0));
	ebrys.Push   (make_tuple(-40, "fy",  -10.0));
	ebrys.Push   (make_tuple(-40, "pwp",   0.0));
	FEM::SetBrys (&mesh, NULL, &ebrys, NULL, &g);
	CallSolve    (1, sol);

	// Stage # 2
	sol->SetNumDiv(10)->SetDeltaTime(16666.0);
	ebrys.Resize (0);
	ebrys.Push   (make_tuple(-10, "ux",    0.0));
	ebrys.Push   (make_tuple(-20, "ux",    0.0));
	ebrys.Push   (make_tuple(-30, "uy",    0.0));
	ebrys.Push   (make_tuple(-40, "pwp",   0.0));
	FEM::SetBrys (&mesh, NULL, &ebrys, NULL, &g);
	CallSolve    (2, sol);

	delete sol;

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

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
