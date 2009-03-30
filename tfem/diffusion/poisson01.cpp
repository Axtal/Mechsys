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
#include "fem/data.h"
#include "fem/solver.h"
#include "fem/elems/tri3.h"
#include "fem/elems/tri6.h"
#include "fem/diffusionelem.h"
#include "models/diffusions/lindiffusion.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "mesh/unstructured.h"

using std::cout;
using std::endl;
using LinAlg::Vector;
using Util::_4;
using Util::_8s;

double SourceFun(double u)
{
	return 1.0;
}

#define T boost::make_tuple

int main(int argc, char **argv) try
{
	// Constants
	size_t ndiv     = 360;               // number of divisions on the circumference of the unit circle
	double dalpha   = 2.0*Util::PI/ndiv; // delta alpha (increment angle on the circumference)
	double maxarea1 = 0.01;              // maximum area for each triangle in the inner circle
	double maxarea2 = 0.1;               // maximum area for each triangle in the outer circle
	bool   is_o2    = false;             // use high order elements?
	String linsol("UM");                 // UMFPACK

	// Input
	cout << "Input: " << argv[0] << "  is_o2  maxarea1_inner(0.01)  maxarea2_outer(0.1)  linsol(LA,UM,SLU)\n";
	if (argc>=2) is_o2      = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) maxarea1   =  atof(argv[2]);
	if (argc>=4) maxarea2   =  atof(argv[3]);
	if (argc>=5) linsol.Printf("%s",argv[4]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Polygon
	Mesh::Unstructured mesh(/*Is3D*/false);
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
	if (is_o2) mesh.SetO2();
	mesh.Generate (true);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Data and Solver
	FEM::Data   dat (2);
	FEM::Solver sol (dat,"tpoisson01");

	// Elements attributes
	FEM::EAtts_T eatts(1);
	String geom; geom = (is_o2 ? "Tri6" : "Tri3");
	eatts = T(-1, geom.CStr(), "Diffusion", "LinDiffusion", "k=1.0", "", "s=0", &SourceFun, true);

	// Set geometry: nodes and elements
	dat.SetNodesElems (&mesh, &eatts);

	// Stage # 1 -----------------------------------------------------------
	FEM::EBrys_T ebrys;
	ebrys.Push        (T(-10, "u", 0.0));
	dat.SetBrys       (&mesh, NULL, &ebrys, NULL);
	dat.AddVolForces  ();
	sol.SolveWithInfo ();

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Check
	Array<double> err_u;
	for (size_t i=0; i<dat.NNodes(); ++i)	
	{
		double x     = dat.Nod(i)->X();
		double y     = dat.Nod(i)->Y();
		double u     = dat.Nod(i)->Val("u");
		double ucorr = (1.0-x*x-y*y)/4.0;
		err_u.Push ( fabs(u-ucorr) / (1.0+fabs(ucorr)) );
	}

	// Error summary
	double tol_u     = 1.0e-2;
	double min_err_u = err_u[err_u.Min()];
	double max_err_u = err_u[err_u.Max()];
	cout << _4<< ""  << _8s<<"Min"     << _8s<<"Mean"                                                  << _8s<<"Max"                << _8s<<"Norm"       << endl;
	cout << _4<< "u" << _8s<<min_err_u << _8s<<err_u.Mean() << (max_err_u>tol_u?"[1;31m":"[1;32m") << _8s<<max_err_u << "[0m" << _8s<<err_u.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err_u>tol_u) return 1;
	else return 0;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }