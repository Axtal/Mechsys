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

/* J N Reddy's FEM book -- Example 8.5.1, page 464
 *  
 * Steady-state heat conduction in an isotropic region
 *
 *               y ^
 *                 | 
 *                    u=T=T0*cos(Pi*x/6a)
 *                #@------@------@------@
 *                #|      |      |      |
 *                #|      |      |      |
 *                #|      |  a   |      |
 *    (insulated) #@------@------@------@  u=T=0
 *       q=0.0    #|      |      |      |
 *                #|      |a     |      |  edgetag=-20
 *   edgetag=-10  #|      |      |      |
 *                #@------@------@------@  --> x
 *                #######################
 *                      (insulated)
 *                         q=0.0
 *
 *                      edgetag=-10
 */

// STL
#include <iostream>
#include <cmath>

// MechSys
#include "fem/geometry.h"
#include "fem/functions.h"
#include "fem/elems/quad4diffusion.h"
#include "models/diffusions/lindiffusion.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "fem/output.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "util/util.h"

using boost::make_tuple;
using std::cout;
using std::endl;
using Util::PI;
using Util::_4;
using Util::_6;
using Util::_8s;

int main(int argc, char **argv) try
{
	// Constants
	double k  = 1.0;     // isotropic conductivity
	double a  = 1.0;     // cell size
	double T0 = 1.0;     // initial temperature (at top)
	double L  = 3.0*a;   // total length
	double H  = 2.0*a;   // total height
	int    nx = 3;       // number of divisions along x
	int    ny = 2;       // number of divisions along y
	String linsol("UM"); // linear solver: UMFPACK

	// Input
	cout << "Input: " << argv[0] << "  linsol(LA,UM,SLU)\n";
	if (argc==2) linsol.Printf("%s",argv[1]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Blocks
	Mesh::Block b;
	b.SetTag    (-1); // tag to be replicated to all generated elements inside this block
	b.SetCoords (false, 4,              // Is3D, NNodes
	             0.,  L, L, 0.,         // x coordinates
	             0., 0., H,  H);        // y coordinates
	b.SetNx     (nx);                   // x weights and num of divisions along x
	b.SetNy     (ny);                   // y weights and num of divisions along y
	b.SetETags  (4,  -10, -20, -10, 0); // edge tags
	Array<Mesh::Block*> blocks;
	blocks.Push (&b);

	// Generate
	Mesh::Structured mesh(/*Is3D*/false);
	clock_t start = std::clock(); // Initial time
	size_t  ne    = mesh.Generate (blocks);
	clock_t total = std::clock() - start; // Time elapsed
	cout << "\nNumber of quadrangles   = " << ne << endl;
	cout << "Time elapsed (mesh)     = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";

	////////////////////////////////////////////////////////////////////////////////////////// FEA /////

	// Geometry
	FEM::Geom g(2); // 2D

	// Edges brys
	FEM::EBrys_T ebrys;
	ebrys.Push (make_tuple(-10, "q", 0.0)); // tag, key, val
	ebrys.Push (make_tuple(-20, "u", 0.0)); // tag, key, val

	// Elements attributes
	String prms; prms.Printf("k=%f", k);
	FEM::EAtts_T eatts;
	eatts.Push (make_tuple(-1, "Quad4Diffusion", "LinDiffusion", prms.CStr(), "", "", true)); // tag, type, model, prms, inis, props

	// Set geometry: nodes, elements, attributes, and boundaries
	FEM::SetNodesElems (&mesh, &eatts, &g);
	FEM::SetBrys       (&mesh, NULL, &ebrys, NULL, &g);

	// Set upper nodes boundary condition
	for (size_t i=0; i<g.NNodes(); ++i)
	{
		double x = g.Nod(i)->X();
		double y = g.Nod(i)->Y();
		if (fabs(y-H)<1.0e-5) // top node
			g.Nod(i)->Bry ("u", T0*cos(PI*x/(6.0*a)));
	}

	// Check conductivity matrices
	double max_err_ke = 0.0;
	LinAlg::Matrix<double> Ke_correct;  Ke_correct.Resize(4,4);
	Ke_correct =  4.0, -1.0, -2.0, -1.0,
				 -1.0,  4.0, -1.0, -2.0,
				 -2.0, -1.0,  4.0, -1.0,
				 -1.0, -2.0, -1.0,  4.0;
	Ke_correct = (k/6.0)*Ke_correct;
	for (size_t i=0; i<g.NElems(); ++i)
	{
		LinAlg::Matrix<double> Ke;
		g.Ele(i)->Order1Matrix(0,Ke);
		double err_ke = 0.0;
		for (int i=0; i<4; ++i)
		for (int j=0; j<4; ++j)
			err_ke += fabs(Ke(i,j)-Ke_correct(i,j));
		if (err_ke>max_err_ke) max_err_ke = err_ke;
	}
	if (max_err_ke>3.4e-5) throw new Fatal("tex831: max_err_ke==%e for quadrangular mesh is bigger than %e.",max_err_ke,3.4e-5);

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol->SetGeom(&g)->SetLinSol(linsol.CStr());
	sol->SolveWithInfo();
	delete sol;

	// Output: Nodes
	cout << _6<<"Node #" << _8s<<"u" << _8s<<"q" << endl;
	for (size_t i=0; i<g.NNodes(); ++i)
		cout << _6<<i << _8s<<g.Nod(i)->Val("u") << _8s<<g.Nod(i)->Val("q") << endl;
	cout << endl;

	//////////////////////////////////////////////////////////////////////////////////////// Check /////
	
	// Check
	Array<double> err_u;
	for (size_t i=0; i<g.NNodes(); ++i)	
	{
		double x     = g.Nod(i)->X();
		double y     = g.Nod(i)->Y();
		double u     = g.Nod(i)->Val("u");
		double ucorr = T0*cosh(PI*y/(6.0*a))*cos(PI*x/(6.0*a))/cosh(PI/3.0);
		err_u.Push ( fabs(u-ucorr) / (1.0+fabs(ucorr)) );
	}

	// Error summary
	double tol_u     = 7.5e-3;
	double min_err_u = err_u[err_u.Min()];
	double max_err_u = err_u[err_u.Max()];
	cout << _4<< ""  << _8s<<"Min"     << _8s<<"Mean"                                                  << _8s<<"Max"                << _8s<<"Norm"       << endl;
	cout << _4<< "u" << _8s<<min_err_u << _8s<<err_u.Mean() << (max_err_u>tol_u?"[1;31m":"[1;32m") << _8s<<max_err_u << "[0m" << _8s<<err_u.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err_u>tol_u) return 1;
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
