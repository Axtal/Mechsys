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
#include "fem/elems/quad4pstrain.h"
#include "fem/elems/quad4diffusion.h"
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "util/exception.h"
#include "util/numstreams.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	// Geometry
	FEM::Geom g(2);

	// Nodes
	g.SetNNodes (4);
	g.SetNode   (0, 0.0, 0.0);
	g.SetNode   (1, 1.0, 0.0);
	g.SetNode   (2, 1.0, 1.0);
	g.SetNode   (3, 0.0, 1.0);

	// Elements
	g.SetNElems (1);
	g.SetElem   (0, "Quad4PStrain", /*IsActive*/true);

	// Set connectivity
	g.Ele(0)->Connect(0, g.Nod(0))
	        ->Connect(1, g.Nod(1))
	        ->Connect(2, g.Nod(2))
	        ->Connect(3, g.Nod(3));

	// Parameters and initial values
	g.Ele(0)->SetModel("LinElastic", "E=100.0 nu=0.25", "Sx=0.0 Sy=0.0 Sxy=0.0");

	// Boundary conditions (must be after connectivity)
	g.Nod(0)->Bry("ux",0.0)->Bry("uy",0.0);
	g.Nod(1)->Bry("uy",0.0);

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol -> SetGeom(&g) -> SetLinSol("UM") -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();
	delete sol;

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
