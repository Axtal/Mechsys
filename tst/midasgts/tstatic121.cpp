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

/*                      
      ///@-----------@ |
      ///|           | | F
      ///|           | |
      ///@-----------@ V
*/

// STL
#include <iostream>

// MechSys
#include "fem/geometry.h"
#include "fem/functions.h"
#include "fem/elems/quad4pstress.h"
#include "fem/elems/quad8pstress.h"
#include "fem/elems/quad4pstrain.h"
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "fem/output.h"
#include "util/exception.h"
#include "linalg/matrix.h"
#include "mesh/unstructured.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_8s;
using boost::make_tuple;

int main(int argc, char **argv) try
{
	// Constants
	double E       = 10.7e+6; // Young (psi)
	double nu      = 0.3;     // Poisson
	double F       = 120.0;   // Load
	size_t nx      = 20;       // ndivs along x
	size_t ny      = 10;       // ndivs along y
	bool   is_o2   = false;   // use high order elements?

	// Input
	cout << "Input: " << argv[0] << "  is_o2  maxarea\n";
	if (argc>=2) is_o2 = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) nx    =  atof(argv[2]);
	if (argc>=4) ny    =  atof(argv[3]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Dimensions
	double L = 3.0;
	double H = 0.6;
	double B = 0.1;
	// Blocks
	Mesh::Block b;
	b.SetTag    (-1); // tag to be replicated to all generated elements inside this block
	b.SetCoords (false, 4,            // Is3D, NNodes
	             0.0, 3.0, 3.0, 0.0,  // x coordinates
	             0.0, 0.0, 0.6, 0.6); // y coordinates
	b.SetNx     (nx);                 // x weights and num of divisions along x
	b.SetNy     (ny);                 // y weights and num of divisions along y
	b.SetETags  (4, -10, -20, 0, 0);  // edge tags
	Array<Mesh::Block*> blocks;
	blocks.Push (&b);

	// Generate
	Mesh::Structured mesh(/*Is3D*/false);
	if (is_o2) mesh.SetO2();
	mesh.SetBlocks (blocks);
	mesh.Generate  (true);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Geom g(2); // 2D

	// Edges brys
	FEM::EBrys_T ebrys;
	ebrys.Push (make_tuple(-10, "ux", 0.0));
	ebrys.Push (make_tuple(-10, "uy", 0.0));
	ebrys.Push (make_tuple(-20, "fy", -F/0.06));

	// Elements attributes
	String prms; prms.Printf("E=%f nu=%f",E,nu);
	FEM::EAtts_T eatts;
	if (is_o2) eatts.Push (make_tuple(-1, "Quad8PStrain", "LinElastic", prms.CStr(), "ZERO", "", true));
	else       eatts.Push (make_tuple(-1, "Quad4PStress", "LinElastic", prms.CStr(), "ZERO", "", true));

	// Set geometry: nodes, elements, attributes, and boundaries
	FEM::SetNodesElems (&mesh, &eatts, &g);
	FEM::SetBrys       (&mesh, NULL,   &ebrys, NULL, &g);

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol->SetGeom(&g);
	sol->SolveWithInfo();
	delete sol;

	// Output: VTU
	Output o; o.VTU (&g, "tstatic121.vtu");
	cout << "[1;34mFile <tstatic121.vtu> saved.[0m\n\n";

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Check
    Array<double> err_dis_z;
	double I = B*H*H*H/12.0; // Section momentum of inertia

	for (size_t i=0; i<g.NNodes(); i++)
	{
		double x   = g.Nod(i)->X();
		double uy  = g.Nod(i)->Val("uy");
		double auy = F/(E*I)*(x*x*x/6.0 - x*x*L/2.0); // Analitic solution for the vertical displacement
		err_dis_z.Push(fabs(uy - auy));
		cout << "iNode: " << i << "  Uy: " << uy << "    Analytic Uy: " <<  auy << endl;
	}

	// Error summary
	double tol_dis     = 1.0e-3;
	double min_err_dis = err_dis_z[err_dis_z.Min()];
	double max_err_dis = err_dis_z[err_dis_z.Max()];
	cout << _4<< ""    << _8s<<"Min"       << _8s<<"Mean"                                                          << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	cout << _4<< "Dis" << _8s<<min_err_dis << _8s<<err_dis_z.Mean() << (max_err_dis>tol_dis?"[1;31m":"[1;32m") << _8s<<max_err_dis << "[0m" << _8s<<err_dis_z.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err_dis>tol_dis) return 1;
	else return 0;
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
