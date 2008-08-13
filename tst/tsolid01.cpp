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

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	// constants
	double E  = 20.0; // Young
	double nu = 0.001;  // Poisson
	double q  = 10.0;   // Downward vertical pressure

	/*  Unit cube
	 
	      z
	      |__y      4________________7
	   x,'        ,'|              ,'|
	            ,'               ,'  |
	          ,'    |          ,'    |
	        ,'      .        ,'      |
	      5'_______________6'        |
	      |                |         |
	      |         |      |         |
	      |         0 -  - | -  -  - 3
	      |       ,        |       ,' 
	      |     ,          |     ,'   
	      |   ,            |   ,'     
	      | ,              | ,'       
	      1________________2'         
	*/

	// Input
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);
	else cout << "[1;32mYou may call this program as in:\t " << argv[0] << " LinSol\n  where LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n [0m[1;34m Now using LA (LAPACK)\n[0m" << endl;

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// 0) Geometry
	FEM::Geom g(3); // 3D

	// 1) Nodes
	g.SetNNodes (8);
	g.SetNode   (0, 0.0, 0.0, 0.0);
	g.SetNode   (1, 1.0, 0.0, 0.0);
	g.SetNode   (2, 1.0, 1.0, 0.0);
	g.SetNode   (3, 0.0, 1.0, 0.0);
	g.SetNode   (4, 0.0, 0.0, 1.0);
	g.SetNode   (5, 1.0, 0.0, 1.0);
	g.SetNode   (6, 1.0, 1.0, 1.0);
	g.SetNode   (7, 0.0, 1.0, 1.0);

	// 2) Elements
	g.SetNElems (1);
	g.SetElem   (0, "Hex8Equilib", /*IsActive*/true);

	// 3) Set connectivity
	g.Ele(0)->Connect(0, g.Nod(0))
	        ->Connect(1, g.Nod(1))
	        ->Connect(2, g.Nod(2))
	        ->Connect(3, g.Nod(3))
	        ->Connect(4, g.Nod(4))
	        ->Connect(5, g.Nod(5))
	        ->Connect(6, g.Nod(6))
	        ->Connect(7, g.Nod(7));

	// 4) Boundary conditions (must be after connectivity)
	g.Nod(0)->Bry("ux",0.0)->Bry("uy",0.0)->Bry("uz",0.0);
	g.Nod(1)               ->Bry("uy",0.0)->Bry("uz",0.0);
	g.Nod(2)                              ->Bry("uz",0.0);
	g.Nod(3)->Bry("ux",0.0)               ->Bry("uz",0.0);
	g.Nod(4)->Bry("ux",0.0)->Bry("uy",0.0);
	g.Nod(5)               ->Bry("uy",0.0);
	g.Nod(7)->Bry("ux",0.0);
                                          
	// 5) Parameters and initial values
	String prms; prms.Printf("E=%f nu=%f",E,nu);
	g.Ele(0)->SetModel ("LinElastic", prms.CStr(), "ZERO");
	g.Ele(0)->FaceBry  ("fz", -q, 5); // 5 => top face

	// 6) Solve
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

	if (fabs(err1)>1.0e-14) cout << "[1;31m\nErrors(" << linsol << ") stress/strain = " << err1 << "[0m";
	else                    cout << "[1;32m\nErrors(" << linsol << ") stress/strain = " << err1 << "[0m";
	if (fabs(err2)>1.0e-14) cout << "[1;31m\nErrors(" << linsol << ") displacements = " << err2 << "[0m\n" << endl;
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
