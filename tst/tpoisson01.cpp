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
#include <fstream>
#include <cmath>

// MechSys
#include "fem/geometry.h"
#include "fem/functions.h"
#include "fem/elems/tri3diffusion.h"
#include "fem/elems/tri6diffusion.h"
#include "models/diffusions/lindiffusion.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "fem/output.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "mesh/unstructured.h"

using std::cout;
using std::endl;
using LinAlg::Vector;
using boost::make_tuple;
using Util::_8s;

int main(int argc, char **argv) try
{
	// Constants
	size_t ndiv     = 36;                // number of divisions on the circumference of the unit circle
	double dalpha   = 2.0*Util::PI/ndiv; // delta alpha (increment angle on the circumference)
	double maxarea1 = 0.01;              // maximum area for each triangle in the inner circle
	double maxarea2 = 0.1;               // maximum area for each triangle in the outer circle

	// Input
	String linsol("LA");
	if (argc>=2) maxarea1 = atof(argv[1]);
	if (argc>=3) maxarea2 = atof(argv[2]);
	if (argc>=4) linsol.Printf("%s",argv[3]);
	else cout << "[1;32mYou can call this program as in:\t " << argv[0] << " LinSol\n  where LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n [0m[1;34m Now using LA (LAPACK)\n[0m" << endl;

	// Generate mesh
	Mesh::Unstructured mesh;
	mesh.SetPolySize (/*NPoints*/2*ndiv, /*NSegments*/2*ndiv, /*NRegions*/2);

	// Inner circle
	double ri = 0.7;
	for (size_t i=0; i<ndiv; ++i)
	{
		double alpha = i*dalpha;
		mesh.SetPolyPoint   (i, /*X*/ri*cos(alpha), /*Y*/ri*sin(alpha));
		mesh.SetPolySegment (i, /*iPointLeft*/i, /*iPointRight*/(i==ndiv-1 ? 0 : i+1));
	}

	// Outer circle r=1
	for (size_t i=0; i<ndiv; ++i)
	{
		double alpha = i*dalpha;
		mesh.SetPolyPoint   (i+ndiv, /*X*/cos(alpha), /*Y*/sin(alpha));
		mesh.SetPolySegment (i+ndiv, /*iPointLeft*/i+ndiv, /*iPointRight*/(i==ndiv-1 ? ndiv : i+ndiv+1), /*Tag*/-10);
	}

	// Generate
	mesh.SetPolyRegion (0, /*Tag*/-1, maxarea1, /*X*/0.0, /*Y*/0.0);
	mesh.SetPolyRegion (1, /*Tag*/-1, maxarea2, /*X*/0.9, /*Y*/0.0);
	//mesh.SetO2    ();
	mesh.Generate ();

	// Geometry
	FEM::Geom g(2);

	// Edges brys (the order matters!)
	FEM::EBrys_T ebrys;
	ebrys.Push (make_tuple(-10, "u", 0.0));

	// Elements attributes
	FEM::EAtts_T eatts;
	//eatts.Push (make_tuple(-1, "Tri6Diffusion", "LinDiffusion", "k=1.0", ""));
	eatts.Push (make_tuple(-1, "Tri3Diffusion", "LinDiffusion", "k=1.0", ""));

	// Set geometry: nodes, elements, attributes, and boundaries
	FEM::SetNodesElems (&mesh, &eatts, &g);
	FEM::SetBrys       (&mesh, NULL, &ebrys, NULL, &g);

	// Set heat source
	Array<double> source(1); source.SetValues(1.0);
	for (size_t i=0; i<g.NElems(); ++i)
		g.Ele(i)->SetProps(source);

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol -> SetGeom(&g) -> SetLinSol(linsol.CStr()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();
	double norm_resid = LinAlg::Norm(sol->Resid());
	cout << "Norm(Resid=DFext-DFint) = " << norm_resid << "\n\n";

	// Check
	double errors = 0.0;
	for (size_t i=0; i<g.NNodes(); ++i)	
	{
		double x     = g.Nod(i)->X();
		double y     = g.Nod(i)->Y();
		double u     = g.Nod(i)->Val("u");
		double ucorr = (1.0-x*x-y*y)/4.0;
		errors += fabs(u-ucorr);
	}

	// Output: VTU
	Output o; o.VTU (&g, "tpoisson01.vtu");

	if (fabs(errors)>1.0e-14) cout << "[1;31mErrors(" << linsol << ") = " << errors << "[0m\n" << endl;
	else                      cout << "[1;32mErrors(" << linsol << ") = " << errors << "[0m\n" << endl;

	if (fabs(errors)>1.0e-14) return 1;

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
