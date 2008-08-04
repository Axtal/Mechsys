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
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "fem/output.h"
#include "util/exception.h"
#include "linalg/matrix.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;

int main(int argc, char **argv) try
{
	double H  = 2.0;   // height
	double L  = 2.0;   // length
	double E  = 207.0; // Young
	double nu = 0.3;   // Poisson
	double q  = 1.0;   // Load

	/*        | | | | | | | | | | | | | | | | |  q
	          V V V V V V V V V V V V V V V V V 
	  -+-     @-------@-------@-------@-------@
	   |      |20     |21     |22     |23     |24
	   |      |   e12 |   e13 |   e14 |   e15 |
	   |      |15     |16     |17     |18     |19
	   |      @-------@-------@-------@-------@
	   |      |       |       |       |       |
	   |      |   e8  |   e9  |   e10 |   e11 |
	          |10     |11     |12     |13     |14
	   H      @-------@-------@-------@-------@
	          |       |       |       |       |
	   |      |   e4  |   e5  |   e6  |   e7  |
	   |      |5      |6      |7      |8      |9
	   |      @-------@-------@-------@-------@
	   |      |       |       |       |       |
	   |      |   e0  |   e1  |   e2  |   e3  |
	   |      |0      |1      |2      |3      |4
	  -+-     @-------@-------@-------@-------@
	         /_\     /_\     /_\     /_\     /_\
	         o o     o o     ///     o o     o o
	
	          |-------------- L --------------|
	 */

	// Input
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);
	else cout << "[1;32mYou may call this program as in:\t " << argv[0] << " LinSol\n  where LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n [0m[1;34m Now using LA (LAPACK)\n[0m" << endl;

	// 0) Geometry
	FEM::Geom g(2); // 2D

	// 1) Nodes
	double xmin  = 0.0;
	double ymin  = 0.0;
	int    ndivx = 20;   // must be even for this problem
	int    ndivy = 30;
	double dx    = L/ndivx;
	double dy    = H/ndivy;
	int    nn    = (ndivx+1)*(ndivy+1);
	size_t k     = 0;
	g.SetNNodes (nn);
	for (int j=0; j<ndivy+1; ++j)
	for (int i=0; i<ndivx+1; ++i)
	{
		g.SetNode (k, xmin+i*dx, ymin+j*dy);
		k++;
	}

	// 2,3) Elements and connectivity
	int ne = ndivx*ndivy;
	g.SetNElems (ne);
	k = 0;
	for (int j=0; j<ndivy; ++j)
	for (int i=0; i<ndivx; ++i)
	{
		int I =  i    +  j   *(ndivx+1);
		int J = (i+1) +  j   *(ndivx+1);
		int K = (i+1) + (j+1)*(ndivx+1);
		int L =  i    + (j+1)*(ndivx+1);
		g.SetElem (k, "Quad4PStrain")->Connect(0, g.Nod(I))
		                             ->Connect(1, g.Nod(J))
		                             ->Connect(2, g.Nod(K))
		                             ->Connect(3, g.Nod(L));
		k++;
	}

	// 4) Boundary conditions (must be after connectivity)
	for (size_t i=0; i<g.NNodes(); ++i)
	{
		if (fabs(g.Nod(i)->Y()-ymin)<1.0e-5) // bottom nodes
		{
			g.Nod(i)->Bry("uy",0.0);
			if (fabs(g.Nod(i)->X()-(xmin+L/2.0))<1.0e-5) // central node
				g.Nod(i)->Bry("ux",0.0);
		}
	}
	for (size_t i=0; i<g.NElems(); ++i)
	{
		Matrix<double> c;
		g.Ele(i)->Coords(c);
		if (fabs(c(2,1)-H)<1.0e-5) // top element
			g.Ele(i)->EdgeBry("fy", -q, 3); // 3=> top edge
	}

	// 5) Parameters and initial values
	String prms; prms.Printf("E=%f  nu=%f",E,nu);
	for (size_t i=0; i<g.NElems(); ++i)
		g.Ele(i)->SetModel("LinElastic", prms.CStr(), "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0");

	// 6) Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	//FEM::Solver * sol = FEM::AllocSolver("AutoME");
	sol -> SetGeom(&g) -> SetLinSol(linsol.CStr()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();
	cout << "NormResid = " << sol->GetVar("NormResid") << endl;

	// Output
	cout << "Node 20: ux = " << g.Nod(20)->Val("ux") << " : uy = " << g.Nod(20)->Val("uy") << " : fx = "  << g.Nod(20)->Val("fx")  << " : fy = "  << g.Nod(20)->Val("fy")  << endl;
	cout << "Node 22: ux = " << g.Nod(22)->Val("ux") << " : uy = " << g.Nod(22)->Val("uy") << " : fx = "  << g.Nod(22)->Val("fx")  << " : fy = "  << g.Nod(22)->Val("fy")  << endl;
	cout << "Node 24: ux = " << g.Nod(24)->Val("ux") << " : uy = " << g.Nod(24)->Val("uy") << " : fx = "  << g.Nod(24)->Val("fx")  << " : fy = "  << g.Nod(24)->Val("fy")  << endl << endl;
	cout << "Node 0:  ux = " << g.Nod( 0)->Val("ux") << " : uy = " << g.Nod( 0)->Val("uy") << " : fx = "  << g.Nod( 0)->Val("fx")  << " : fy = "  << g.Nod( 0)->Val("fy")  << endl;
	cout << "Node 2:  ux = " << g.Nod( 2)->Val("ux") << " : uy = " << g.Nod( 2)->Val("uy") << " : fx = "  << g.Nod( 2)->Val("fx")  << " : fy = "  << g.Nod( 2)->Val("fy")  << endl;
	cout << "Node 4:  ux = " << g.Nod( 4)->Val("ux") << " : uy = " << g.Nod( 4)->Val("uy") << " : fx = "  << g.Nod( 4)->Val("fx")  << " : fy = "  << g.Nod( 4)->Val("fy")  << endl << endl;
	cout << "Elem 9:  Sx = " << g.Ele( 9)->Val("Sx") << " : Sy = " << g.Ele( 9)->Val("Sy") << " : Sxy = " << g.Ele( 9)->Val("Sxy") << endl;
	cout << "Elem 9:  Ex = " << g.Ele( 9)->Val("Ex") << " : Ey = " << g.Ele( 9)->Val("Ey") << " : Exy = " << g.Ele( 9)->Val("Exy") << endl << endl;

	// Check
    double errors = 0.0;

	double Sy = q;
	double Ex = -nu*(1.0+nu)*Sy/E;
	double Ey =  (1.0-nu*nu)*Sy/E;
	double Sz = (E/(1.0+nu))*(nu/(1.0-2.0*nu))*(Ex+Ey);

	// Stress and strains
	for (size_t i=0; i<g.NElems(); ++i)
	{
		errors += fabs(g.Ele(i)->Val("Ex" ) - (Ex));
		errors += fabs(g.Ele(i)->Val("Ey" ) - (Ey));
		errors += fabs(g.Ele(i)->Val("Exy") - (0.0));
		errors += fabs(g.Ele(i)->Val("Sx" ) - (0.0));
		errors += fabs(g.Ele(i)->Val("Sy" ) - (Sy ));
		errors += fabs(g.Ele(i)->Val("Sz" ) - (Sz ));
		errors += fabs(g.Ele(i)->Val("Sxy") - (0.0));
	}

	// Displacements
	for (size_t i=0; i<g.NNodes(); ++i)
	{
		if (fabs(g.Nod(i)->Y()-ymin)<1.0e-5) // bottom nodes
		{
			errors += fabs(g.Nod(i)->Val("uy") - (0.0));
			if (fabs(g.Nod(i)->X()-(xmin+L/2.0))<1.0e-5) // central node
				errors += fabs(g.Nod(i)->Val("ux") - (0.0));
		}
		else if (fabs(g.Nod(i)->Y()-(ymin+H/2.0))<1.0e-5) // mid nodes
			errors += fabs(g.Nod(i)->Val("uy") - (-0.5*H*Ey));
		else if (fabs(g.Nod(i)->Y()-(ymin+H))<1.0e-5) // top nodes
			errors += fabs(g.Nod(i)->Val("uy") - (-H*Ey));
		if (fabs(g.Nod(i)->X()-(xmin))<1.0e-5) // left nodes
			errors += fabs(g.Nod(i)->Val("ux") - (0.5*L*Ex));
		if (fabs(g.Nod(i)->X()-(xmin+L))<1.0e-5) // right nodes
			errors += fabs(g.Nod(i)->Val("ux") - (-0.5*L*Ex));
	}

	if (fabs(errors)>1.0e-10) cout << "[1;31m\nErrors(" << linsol << ") = " << errors << "[0m\n" << endl;
	else                      cout << "[1;32m\nErrors(" << linsol << ") = " << errors << "[0m\n" << endl;

	// Write file
	Output o; o.VTK (&g, "tpstrain01.vtk");
	          o.VTU (&g, "tpstrain01.vtu");
	cout << "[1;34mFile <tpstrain01.vtk> saved.[0m" << endl;
	cout << "[1;34mFile <tpstrain01.vtu> saved.[0m" << endl;

	// Return error flag
	if (fabs(errors)>1.0e-10) return 1;
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
