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
      u=100*sin(PI*x/10)
       //@-----------@       ---
       //|           |        |
       //|           |        |
  u=0  //|           |q=0    H=10
       //|           |        |
       //|    L=5    |        |
       //@-----------@       ---
       ///////////////
              u=0     
*/

// STL
#include <iostream>
#include <fstream>
#include <cmath>

// MechSys
#include "fem/geometry.h"
#include "fem/functions.h"
#include "fem/elems/tri3diffusion.h"
#include "fem/elems/tri6diffusion.h"
#include "fem/elems/quad4diffusion.h"
#include "fem/elems/quad8diffusion.h"
#include "models/diffusions/lindiffusion.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "fem/output.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "mesh/structured.h"
#include "mesh/unstructured.h"

using std::cout;
using std::endl;
using LinAlg::Vector;
using boost::make_tuple;
using Util::_4;
using Util::_6;
using Util::_8s;
using Util::PI;

int main(int argc, char **argv) try
{
	// Constants
	bool   is_o2   = false; // use high order elements?
	double maxarea = 1.0;   // maximum area for each triangle
	bool   is_quad = false; // use quadrilateral mesh?
	size_t ndiv    = 4;     // number of divisions
	String linsol("UM");    // UMFPACK

	// Input
	cout << "Input: " << argv[0] << "  is_o2  maxarea(1.0)  is_quad  ndiv  linsol(LA,UM,SLU)\n";
	if (argc>=2) is_o2      = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) maxarea    =  atof(argv[2]);
	if (argc>=4) is_quad    = (atoi(argv[3])>0 ? true : false);
	if (argc>=5) ndiv       =  atof(argv[4]);
	if (argc>=6) linsol.Printf("%s",argv[5]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Mesh
	Mesh::Generic     * msh;
	Mesh::Structured    mss(/*Is3D*/false);
	Mesh::Unstructured  msu(/*Is3D*/false);
	if (is_quad) // quadrangles
	{
		// Blocks
		Mesh::Block b;
		b.SetTag    (-1); // tag to be replicated to all generated elements inside this block
		b.SetCoords (false, 4,           // Is3D, NNodes
					 0., 5.,  5.,  0.,   // x coordinates
					 0., 0., 10., 10.);  // y coordinates
		b.SetNx     (ndiv);              // x weights and num of divisions along x
		b.SetNy     (ndiv);              // y weights and num of divisions along y
		b.SetETags  (4,-10,-20,-10,0);   // edge tags
		Array<Mesh::Block*> blocks;
		blocks.Push (&b);

		// Generate
		if (is_o2) mss.SetO2();
		mss.SetBlocks (blocks);
		mss.Generate  (true);
		msh = &mss;
	}
	else // triangles
	{
		// Polygon
		msu.SetPolySize    (/*NPoints*/4, /*NSegments*/4, /*NRegions*/1);
		msu.SetPolyPoint   (0, /*X*/ 0.0, /*Y*/ 0.0);
		msu.SetPolyPoint   (1, /*X*/ 5.0, /*Y*/ 0.0);
		msu.SetPolyPoint   (2, /*X*/ 5.0, /*Y*/10.0);
		msu.SetPolyPoint   (3, /*X*/ 0.0, /*Y*/10.0);
		msu.SetPolySegment (0, /*iPointLeft*/0, /*iPointRight*/1, /*Tag*/-10);
		msu.SetPolySegment (1, /*iPointLeft*/1, /*iPointRight*/2, /*Tag*/-20);
		msu.SetPolySegment (2, /*iPointLeft*/2, /*iPointRight*/3);
		msu.SetPolySegment (3, /*iPointLeft*/3, /*iPointRight*/0, /*Tag*/-10);
		msu.SetPolyRegion  (0, /*Tag*/-1, maxarea, /*X*/2.5, /*Y*/5.0);

		// Generate
		if (is_o2) msu.SetO2();
		msu.Generate (true);
		msh = &msu;
	}

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Geom g(2); // 2D

	// Edges brys (the order matters!)
	FEM::EBrys_T ebrys;
	ebrys.Push (make_tuple(-20, "q", 0.0));
	ebrys.Push (make_tuple(-10, "u", 0.0));

	// Elements attributes
	FEM::EAtts_T eatts;
	if (is_quad)
	{
		if (is_o2) eatts.Push (make_tuple(-1, "Quad8Diffusion", "LinDiffusion", "k=1.0", "", "", true));
		else       eatts.Push (make_tuple(-1, "Quad4Diffusion", "LinDiffusion", "k=1.0", "", "", true));
	}
	else
	{
		if (is_o2) eatts.Push (make_tuple(-1, "Tri6Diffusion", "LinDiffusion", "k=1.0", "", "", true));
		else       eatts.Push (make_tuple(-1, "Tri3Diffusion", "LinDiffusion", "k=1.0", "", "", true));
	}

	// Set geometry: nodes, elements, attributes, and boundaries
	FEM::SetNodesElems (msh, &eatts, &g);
	FEM::SetBrys       (msh, NULL, &ebrys, NULL, &g);

	// Set prescribed u = 100*sin(PI*x/10) on top nodes
	for (size_t i=0; i<g.NNodes(); ++i)
	{
		double x = g.Nod(i)->X();
		double y = g.Nod(i)->Y();
		if (y>10.0-1e-7) // top node
			g.Nod(i)->Bry("u", 100.0*sin(PI*x/10.0));
	}

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol->SetGeom(&g)->SetLinSol(linsol.CStr());
	sol->SolveWithInfo();
	delete sol;

	// Output: VTU
	Output o; o.VTU (&g, "tlaplace01.vtu");
	cout << "[1;34mFile <tlaplace01.vtu> saved.[0m\n\n";

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Check
	Array<double> err_u;
	//cout << _6<<"Node #" << _8s<<"u" << _8s<<"ucorrect" << _8s<<"q" << endl;
	for (size_t i=0; i<g.NNodes(); ++i)	
	{
		double x     = g.Nod(i)->X();
		double y     = g.Nod(i)->Y();
		double u     = g.Nod(i)->Val("u");
		double ucorr = 100.0*sinh(PI*y/10.0)*sin(PI*x/10.0)/sinh(PI);
		err_u.Push ( fabs(u-ucorr) / (1.0+fabs(ucorr)) );
		//cout << _6<<i << _8s<<g.Nod(i)->Val("u") << _8s<<ucorr << _8s<<g.Nod(i)->Val("q") << endl;
	}
	//cout << endl;

	// Error summary
	double tol_u     = 3.8e-2;
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
