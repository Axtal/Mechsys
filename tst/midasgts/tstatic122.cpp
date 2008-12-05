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
#include "fem/elems/hex8equilib.h"
#include "fem/elems/hex20equilib.h"
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
	size_t nx      = 1;       // ndivs along x
	size_t ny      = 20;    // ndivs along y
	size_t nz      = 10;    // ndivs along z
	bool   is_o2   = false;   // use high order elements?

	// Input
	cout << "Input: " << argv[0] << "  is_o2  maxarea\n";
	if (argc>=2) is_o2 = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) nx    =  atof(argv[2]);
	if (argc>=4) ny    =  atof(argv[3]);
	if (argc>=5) nz    =  atof(argv[4]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Dimensions
	double L = 3.0;
	double H = 0.6;
	double B = 0.1;

	// Blocks
	Mesh::Block b;
	b.SetTag    (-1); // tag to be replicated to all generated elements inside this block
	b.SetCoords (true, 8,             // Is3D, NNodes
	             0.0, B, B, 0.0, 0.0, B, B, 0.0,  // x coordinates
	             0.0, 0.0, L, L, 0.0, 0.0, L, L,  // y coordinates
	             0.0, 0.0, 0.0, 0.0, H, H, H, H); // z coordinates
	b.SetNx     (nx);                 // x weights and num of divisions along x
	b.SetNy     (ny);                 // y weights and num of divisions along y
	b.SetNz     (nz);                 // z weights and num of divisions along z
	b.SetFTags  (6, -10, -20, -30, -40, 0, 0);  // face tags
	Array<Mesh::Block*> blocks;
	blocks.Push (&b);

	// Generate
	Mesh::Structured mesh(/*Is3D*/true);
	if (is_o2) mesh.SetO2();
	mesh.SetBlocks (blocks);
	mesh.Generate  (true);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Geom g(3); // 2D

	// Face brys
	FEM::FBrys_T fbrys;
	fbrys.Push (make_tuple(-30, "ux", 0.0));
	fbrys.Push (make_tuple(-30, "uy", 0.0));
	fbrys.Push (make_tuple(-30, "uz", 0.0));
	fbrys.Push (make_tuple(-40, "fz", -F/0.06));

	// Elements attributes
	String prms; prms.Printf("E=%f nu=%f",E,nu);
	FEM::EAtts_T eatts;
	if (is_o2) eatts.Push (make_tuple(-1, "Hex20Equilib", "LinElastic", prms.CStr(), "ZERO", "", true));
	else       eatts.Push (make_tuple(-1, "Hex8Equilib" , "LinElastic", prms.CStr(), "ZERO", "", true));

	// Set geometry: nodes, elements, attributes, and boundaries
	FEM::SetNodesElems (&mesh, &eatts, &g);
	FEM::SetBrys       (&mesh, NULL, NULL, &fbrys, &g);

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol->SetGeom(&g);
	sol->SolveWithInfo();
	delete sol;

	// Output: VTU
	Output o; o.VTU (&g, "tstatic122.vtu");
	cout << "[1;34mFile <tstatic122.vtu> saved.[0m\n\n";

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Check
    Array<double> err_dis_z;
	double I = B*H*H*H/12.0; // Section momentum of inertia

	for (size_t i=0; i<g.NNodes(); i++)
	{
		double y   = g.Nod(i)->Y();
		double uz  = g.Nod(i)->Val("uz");
		double auz = F/(E*I)*(y*y*y/6.0 - y*y*L/2.0); // Analitic solution for the vertical displacement
		err_dis_z.Push(fabs(uz - auz));
		cout << "iNode: " << i << "  Uz: " << uz << "    Analytic Uz: " <<  auz << endl;
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
