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
#include "fem/data.h"
#include "fem/node.h"
#include "fem/elems/quad4pstrain.h"
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "util/exception.h"
#include "linalg/matrix.h"

using std::cout;
using std::endl;
using FEM::Nodes;
using FEM::Elems;
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

	// 0) Problem dimension
	FEM::Dim = 2; // 2D

	// 1) Nodes
	double xmin  = 0.0;
	double ymin  = 0.0;
	int    ndivx = 4;
	int    ndivy = 4;
	double dx    = L/ndivx;
	double dy    = H/ndivy;
	for (int j=0; j<ndivy+1; ++j)
	for (int i=0; i<ndivx+1; ++i)
		FEM::AddNode (xmin+i*dx, ymin+j*dy);

	// 2,3) Elements and connectivity
	for (int j=0; j<ndivy; ++j)
	for (int i=0; i<ndivx; ++i)
	{
		int I =  i    +  j   *(ndivx+1);
		int J = (i+1) +  j   *(ndivx+1);
		int K = (i+1) + (j+1)*(ndivx+1);
		int L =  i    + (j+1)*(ndivx+1);
		FEM::AddElem("Quad4PStrain")->SetNode(0,I)->SetNode(1,J)->SetNode(2,K)->SetNode(3,L);
	}

	// 4) Boundary conditions (must be after connectivity)
	for (size_t i=0; i<Nodes.Size(); ++i)
	{
		if (fabs(Nodes[i]->Y()-ymin)<1.0e-5) // bottom nodes
		{
			Nodes[i]->Bry("uy",0.0);
			if (fabs(Nodes[i]->X()-L/2.0)<1.0e-5) // central node
				Nodes[i]->Bry("ux",0.0);
		}
	}
	for (size_t i=0; i<Elems.Size(); ++i)
	{
		Matrix<double> c;
		Elems[i]->Coords(c);
		if (fabs(c(2,1)-H)<1.0e-5) // top element
			Elems[i]->BryL("fy", -q, 2, 3,2); // Actually, fy is traction == ty (list of nodes is LOCAL)
	}

	// 5) Parameters and initial values
	String prms; prms.Printf("E=%f  nu=%f",E,nu);
	for (size_t i=0; i<Elems.Size(); ++i)
		Elems[i]->SetModel("LinElastic", prms.GetSTL().c_str(), "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0");

	// 6) Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	//FEM::Solver * sol = FEM::AllocSolver("AutoME");
	sol -> SetLinSol(linsol.GetSTL().c_str()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();

	// Output
	cout << "Node 20: ux = " << Nodes[20]->Val("ux") << " : uy = " << Nodes[20]->Val("uy") << " : fx = "  << Nodes[20]->Val("fx")  << " : fy = "  << Nodes[20]->Val("fy")  << endl;
	cout << "Node 22: ux = " << Nodes[22]->Val("ux") << " : uy = " << Nodes[22]->Val("uy") << " : fx = "  << Nodes[22]->Val("fx")  << " : fy = "  << Nodes[22]->Val("fy")  << endl;
	cout << "Node 24: ux = " << Nodes[24]->Val("ux") << " : uy = " << Nodes[24]->Val("uy") << " : fx = "  << Nodes[24]->Val("fx")  << " : fy = "  << Nodes[24]->Val("fy")  << endl << endl;
	cout << "Node 0:  ux = " << Nodes[ 0]->Val("ux") << " : uy = " << Nodes[ 0]->Val("uy") << " : fx = "  << Nodes[ 0]->Val("fx")  << " : fy = "  << Nodes[ 0]->Val("fy")  << endl;
	cout << "Node 2:  ux = " << Nodes[ 2]->Val("ux") << " : uy = " << Nodes[ 2]->Val("uy") << " : fx = "  << Nodes[ 2]->Val("fx")  << " : fy = "  << Nodes[ 2]->Val("fy")  << endl;
	cout << "Node 4:  ux = " << Nodes[ 4]->Val("ux") << " : uy = " << Nodes[ 4]->Val("uy") << " : fx = "  << Nodes[ 4]->Val("fx")  << " : fy = "  << Nodes[ 4]->Val("fy")  << endl << endl;
	cout << "Elem 9:  Sx = " << Elems[ 9]->Val("Sx") << " : Sy = " << Elems[ 9]->Val("Sy") << " : Sxy = " << Elems[ 9]->Val("Sxy") << endl;
	cout << "Elem 9:  Ex = " << Elems[ 9]->Val("Ex") << " : Ey = " << Elems[ 9]->Val("Ey") << " : Exy = " << Elems[ 9]->Val("Exy") << endl << endl;

	// Check
    double errors = 0.0;

	double Sy = q;
	double Ex = -nu*(1.0+nu)*Sy/E;
	double Ey =  (1.0-nu*nu)*Sy/E;

	errors += fabs(Elems[9]->Val("Ex" ) - (Ex));
	errors += fabs(Elems[9]->Val("Ey" ) - (Ey));
	errors += fabs(Elems[9]->Val("Exy") - (0.0));
	errors += fabs(Elems[9]->Val("Sx" ) - (0.0));
	errors += fabs(Elems[9]->Val("Sy" ) - (Sy ));
	errors += fabs(Elems[9]->Val("Sxy") - (0.0));

	//errors += fabs(Nodes[0]->Val("ux") - ( 0.5*L*Ex));
	//errors += fabs(Nodes[1]->Val("ux") - (-0.5*L*Ex));
	//errors += fabs(Nodes[2]->Val("ux") - (-0.5*L*Ex));
	//errors += fabs(Nodes[3]->Val("ux") - ( 0.5*L*Ex));
	//errors += fabs(Nodes[4]->Val("ux") - (      0.0));
	//errors += fabs(Nodes[5]->Val("ux") - (-0.5*L*Ex));
	//errors += fabs(Nodes[6]->Val("ux") - (      0.0));
	//errors += fabs(Nodes[7]->Val("ux") - ( 0.5*L*Ex));
//
	//errors += fabs(Nodes[0]->Val("uy") - (      0.0));
	//errors += fabs(Nodes[1]->Val("uy") - (      0.0));
	//errors += fabs(Nodes[2]->Val("uy") - (    -H*Ey));
	//errors += fabs(Nodes[3]->Val("uy") - (    -H*Ey));
	//errors += fabs(Nodes[4]->Val("uy") - (      0.0));
	//errors += fabs(Nodes[5]->Val("uy") - (-0.5*H*Ey));
	//errors += fabs(Nodes[6]->Val("uy") - (    -H*Ey));
	//errors += fabs(Nodes[7]->Val("uy") - (-0.5*H*Ey));
//
	//errors += fabs(Nodes[3]->Val("fy") - (    -q*L/6.0));
	//errors += fabs(Nodes[6]->Val("fy") - (-2.0*q*L/3.0));
	//errors += fabs(Nodes[2]->Val("fy") - (    -q*L/6.0));
//
	//errors += fabs(Nodes[3]->Val("fy")+Nodes[6]->Val("fy")+Nodes[2]->Val("fy")-(-q*L));

	if (fabs(errors)>1.0e-13) cout << "[1;31m\nErrors(" << linsol << ") = " << errors << "[0m\n" << endl;
	else                      cout << "[1;32m\nErrors(" << linsol << ") = " << errors << "[0m\n" << endl;

	// Write file
	FEM::WriteVTU ("tpstrain01.vtu");
	cout << "[1;34mFile <tpstrain01.vtu> saved.[0m" << endl;

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
