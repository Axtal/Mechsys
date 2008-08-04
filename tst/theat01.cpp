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
#include "fem/elems/quad4heat.h"
#include "models/heats/linheat.h"
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

int main(int argc, char **argv) try
{
	// Constants
	double k  = 1.0;
	double a  = 1.0;
	double T0 = 1.0;
	double L  = 3.0*a;
	double H  = 2.0*a;
	int    nx = 3;
	int    ny = 2;

	/* J N Reddy's FEM book -- Example 8.5.1, page 464
	 *
	 *               y ^
	 *                 | 
	 *                    T=T0*cos(Pi*x/6a)
	 *                #@------@------@------@
	 *                #|      |      |      |
	 *                #|      |      |      |
	 *                #|      |  a   |      |
	 *    (insulated) #@------@------@------@  T=0
	 *       F=0.0    #|      |      |      |
	 *                #|      |a     |      |  -20
	 *        -10     #|      |      |      |
	 *                #@------@------@------@  --> x
	 *                #######################
	 *                      (insulated)
	 *                         F=0.0
	 *
	 *                         -10
	 */

	// Input
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);
	else cout << "[1;32mYou may call this program as in:\t " << argv[0] << " LinSol\n  where LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n [0m[1;34m Now using LA (LAPACK)\n[0m" << endl;

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Blocks
	String wx; for (int i=0; i<nx; ++i) wx.Printf("%s %f",wx.CStr(),1.0);
	String wy; for (int i=0; i<ny; ++i) wy.Printf("%s %f",wy.CStr(),1.0);
	Mesh::Block b;
	b.SetTag (-1); // tag to be replicated to all generated elements inside this block
	b.Set2D  ();   // 2D
	b.C      () = 0.,  L, L, 0.,    L/2.,    L, L/2.,   0., // x coordinates
	              0., 0., H,  H,      0., H/2.,    H, H/2.; // y coordinates
	b.SetWx  (wx.CStr());                                   // x weights and num of divisions along x
	b.SetWy  (wx.CStr());                                   // y weights and num of divisions along y
	b.ETags  () = -10, -20, -10, 0;                         // edge tags
	Array<Mesh::Block*> blocks;
	blocks.Push (&b);

	// Generate
	cout << "\nMesh Generation: --------------------------------------------------------------" << endl;
	Mesh::Structured ms;
	clock_t start = std::clock(); // Initial time
	size_t  ne    = ms.Generate (blocks);
	clock_t total = std::clock() - start; // Time elapsed
	cout << "[1;33m"<<ne<<" elements[0m. Time elapsed = [1;31m" << static_cast<double>(total)/CLOCKS_PER_SEC << "[0m [1;32mseconds[0m" << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////// FEA /////

	// Geometry
	FEM::Geom g(2); // 2D

	// Edges brys
	FEM::EBrys_T ebrys;
	ebrys.Push (make_tuple(-10, "F", 0.0)); // tag, key, val
	ebrys.Push (make_tuple(-20, "T", 0.0)); // tag, key, val

	// Elements attributes
	String prms; prms.Printf("k=%f", k);
	FEM::EAtts_T eatts;
	eatts.Push (make_tuple(-1, "Quad4Heat", "LinHeat", prms.CStr(), "")); // tag, type, model, prms, inis

	// Set geometry: nodes, elements, attributes, and boundaries
	FEM::SetGeom (&ms, NULL, &ebrys, NULL, &eatts, &g);

	// Set upper nodes boundary condition
	for (size_t i=0; i<g.NNodes(); ++i)
	{
		double x = g.Nod(i)->X();
		double y = g.Nod(i)->Y();
		if (fabs(y-H)<1.0e-5) // top node
			g.Nod(i)->Bry ("T", T0*cos(PI*x/(6.0*a)));
	}

	// Solve
	cout << "\nSolution: ---------------------------------------------------------------------" << endl;
	//FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	FEM::Solver * sol = FEM::AllocSolver("AutoME");
	start = std::clock(); // Initial time
	sol -> SetGeom(&g) -> SetLinSol(linsol.CStr()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> SetCte("DTOL", 1.0e-10);
	sol -> Solve();
	total = std::clock() - start; // Time elapsed
	//cout << "NormResid = "<<sol->GetVar("NormResid")<<". Time elapsed = [1;31m"<<static_cast<double>(total)/CLOCKS_PER_SEC<<"[0m [1;32mseconds[0m"<<std::endl;
	cout << "RelError  = "<<sol->GetVar("RelError" )<<". Time elapsed = [1;31m"<<static_cast<double>(total)/CLOCKS_PER_SEC<<"[0m [1;32mseconds[0m"<<std::endl;

	////////////////////////////////////////////////////////////////////////////////// Output File /////

	// Write file
	cout << "Write output file: ------------------------------------------------------------" << endl;
	start = std::clock(); // Initial time
	Output::VTU o; o.Heat (&g, "theat01.vtu");
	total = std::clock() - start; // Time elapsed
	cout << "[1;34mFile <theat01.vtu> saved.[0m" << endl;
	cout << "Time elapsed = [1;31m"<<static_cast<double>(total)/CLOCKS_PER_SEC<<"[0m [1;32mseconds[0m"<<std::endl;
	cout << endl;

	//////////////////////////////////////////////////////////////////////////////////////// Check /////
	
    double errors = 0.0;

	//if (fabs(errors)>1.0e-10) cout << "[1;31mErrors(" << linsol << ") = " << errors << "[0m\n" << endl;
	//else                      cout << "[1;32mErrors(" << linsol << ") = " << errors << "[0m\n" << endl;

	// Return error flag
	if (fabs(errors)>1.0e-13) return 1;
	else                      return 0;
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
