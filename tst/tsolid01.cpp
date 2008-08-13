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
#include "fem/output.h"
#include "fem/elems/hex8equilib.h"
#include "fem/solvers/autome.h"
#include "fem/solvers/forwardeuler.h"
#include "models/equilibs/linelastic.h"
#include "util/exception.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using boost::make_tuple;

int main(int argc, char **argv) try
{
	// constants
	double E  = 200.0; // Young
	double nu = 0.25;  // Poisson
	double q  = 2.0;   // Downward vertical pressure
	int    nx = 2;     // number of divisions along x
	int    ny = 2;     // number of divisions along y
	int    nz = 2;     // number of divisions along z
	double Lx = 1.0;   // x edge length
	double Ly = 1.0;   // y edge length
	double Lz = 1.0;   // z edge length

	/* Cube with a uniform pressure at the top
	 
	      z
	      |__y      +________________+
	   x,'        ,'|              ,'|
	            ,'               ,'  |
	          ,'    |          ,'    |
	        ,'      .        ,'      | Lz
	      +'_______________+'        |
	      |                |         |
	      |         |      |         |
	      |         + -  - | -  -  - +
	      |       ,        |       ,' 
	      |     ,          |     ,'   
	      |   ,            |   ,'  Lx 
	      | ,       Ly     | ,'       
	      +________________+'         
	*/

	// Input
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);
	else cout << "[1;32mYou may call this program as in:\t " << argv[0] << " LinSol\n  where LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n [0m[1;34m Now using LA (LAPACK)\n[0m" << endl;

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
	
	Mesh::Block b;
	b.SetTag    (-1); // tag to be replicated to all generated elements inside this block
	b.SetCoords (true, 8,                          // Is3D, NNodes
	             0., Lx, Lx, 0.,  0., Lx, Lx, 0.,  // x coordinates
	             0., 0., Ly, Ly,  0., 0., Ly, Ly,  // y coordinates
	             0., 0., 0., 0.,  Lz, Lz, Lz, Lz); // z coordinates
	b.SetFTags  (6, -100,0,-102,0,-104,-105);      // face tags
	b.SetNx     (nx);                              // num of divisions along x
	b.SetNy     (ny);                              // num of divisions along y
	b.SetNz     (nz);                              // num of divisions along z
	Array<Mesh::Block*> blocks;
	blocks.Push (&b);

	// Generate
	cout << "\nMesh Generation: --------------------------------------------------------------" << endl;
	Mesh::Structured ms;
	clock_t start = std::clock(); // Initial time
	size_t  ne    = ms.Generate (blocks);
	clock_t total = std::clock() - start; // Time elapsed
	cout << "[1;33m"<<ne<<" elements[0m. Time elapsed = [1;31m" << static_cast<double>(total)/CLOCKS_PER_SEC << "[0m [1;32mseconds[0m" << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Geom g(3); // 3D

	// Faces brys
	FEM::FBrys_T fbrys;
	fbrys.Push (make_tuple(-100, "ux", 0.0)); // tag, key, val
	fbrys.Push (make_tuple(-102, "uy", 0.0)); // tag, key, val
	fbrys.Push (make_tuple(-104, "uz", 0.0)); // tag, key, val
	fbrys.Push (make_tuple(-105, "fz",  -q)); // tag, key, val

	// Element attributes
	String prms; prms.Printf("E=%f nu=%f",E,nu);
	FEM::EAtts_T eatts;
	eatts.Push (make_tuple(-1, "Hex8Equilib", "LinElastic", prms.CStr(), "ZERO")); // tag, type, model, prms, inis

	// Set geometry
	FEM::SetGeom (&ms, NULL, NULL, &fbrys, &eatts, &g);

	// Solve
	//FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	FEM::Solver * sol = FEM::AllocSolver("AutoME");
	sol -> SetGeom(&g) -> SetLinSol(linsol.CStr()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();

	// Output
	Output o; o.VTU (&g, "tsolid01.vtu");
	cout << "[1;34mFile <tsolid01.vtu> saved.[0m" << endl;

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	double Sx  = 0.0;
	double Sy  = 0.0;
	double Sz  = q;
	double Sxy = 0.0;
	double Syz = 0.0;
	double Szx = 0.0;

	double Ex  = -nu*Sz/E;
	double Ey  = -nu*Sz/E;
	double Ez  = Sz/E;
	double Exy = 0.0;
	double Eyz = 0.0;
	double Ezx = 0.0;

	// Stress and strains
	double err1 = 0.0;
	for (size_t i=0; i<g.NElems(); ++i)
	{
		err1 += fabs(g.Ele(i)->Val("Sx" ) - (Sx ));
		err1 += fabs(g.Ele(i)->Val("Sy" ) - (Sy ));
		err1 += fabs(g.Ele(i)->Val("Sz" ) - (Sz ));
		err1 += fabs(g.Ele(i)->Val("Sxy") - (Sxy));
		err1 += fabs(g.Ele(i)->Val("Syz") - (Syz));
		err1 += fabs(g.Ele(i)->Val("Szx") - (Szx));
		err1 += fabs(g.Ele(i)->Val("Ex" ) - (Ex ));
		err1 += fabs(g.Ele(i)->Val("Ey" ) - (Ey ));
		err1 += fabs(g.Ele(i)->Val("Ez" ) - (Ez ));
		err1 += fabs(g.Ele(i)->Val("Exy") - (Exy));
		err1 += fabs(g.Ele(i)->Val("Eyz") - (Eyz));
		err1 += fabs(g.Ele(i)->Val("Ezx") - (Ezx));
	}

	// Displacements
	double err2 = 0.0;
	for (size_t i=0; i<g.NNodes(); ++i)
	{
		err2 += fabs(g.Nod(i)->Val("ux")-(-Ex*g.Nod(i)->X()));
		err2 += fabs(g.Nod(i)->Val("uy")-(-Ey*g.Nod(i)->Y()));
		err2 += fabs(g.Nod(i)->Val("uz")-(-Ez*g.Nod(i)->Z()));
	}

	if (fabs(err1)>1.0e-13) cout << "[1;31m\nErrors(" << linsol << ") stress/strain = " << err1 << "[0m";
	else                    cout << "[1;32m\nErrors(" << linsol << ") stress/strain = " << err1 << "[0m";
	if (fabs(err2)>1.0e-13) cout << "[1;31m\nErrors(" << linsol << ") displacements = " << err2 << "[0m\n" << endl;
	else                    cout << "[1;32m\nErrors(" << linsol << ") displacements = " << err2 << "[0m\n" << endl;

	// Return error flag
	if (fabs(err1+err2)>1.0e-13) return 1;
	else                         return 0;
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
