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

// MechSys
#include "fem/geometry.h"
#include "fem/functions.h"
#include "fem/elems/hex8equilib.h"
#include "models/equilibs/camclay.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "fem/output.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "linalg/matrix.h"
#include "tensors/tensors.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using boost::make_tuple;
using Tensors::Tensor1;
using Util::_8s;

int main(int argc, char **argv) try
{
	// constants
	double lam   = 0.0891;  // lambda
	double kap   = 0.0196;  // kappa
	double phics = 31.6;    // shear angle at CS
	double G     = 18130.0; // shear modulus (kPa)
	double p_ini = 198.0;   // initial mean stress (kPa)
	double v_ini = 1.6910;  // initial specific volume

	/* Cube -- true triaxial test
	 
	      z
	      |__y      +________________+
	   x,'        ,'|              ,'|
	            ,'               ,'  |
	          ,'    |          ,'    | 1.0
	        ,'      .        ,'      | 
	      +'_______________+'        |
	      |                |         |    13 = # elem in centre
	      |         |      |         |
	      |         + -  - | -  -  - +
	      |       ,        |       ,' 
	      |     ,          |     ,'   
	      |   ,            |   ,'  1.0
	      | ,      1.0     | ,'       
	      +________________+'         
	*/

	// Input
	String linsol("UM");
	if (argc==2) linsol.Printf("%s",argv[1]);
	else cout << "[1;32mYou may call this program as in:\t " << argv[0] << " LinSol\n  where LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n [0m[1;34m Now using LA (LAPACK)\n[0m" << endl;

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
	
	Mesh::Block b;
	b.SetTag    (-1); // tag to be replicated to all generated elements inside this block
	b.SetCoords (true, 8,                           // Is3D, NNodes
	             0., 1., 1., 0.,  0., 1., 1., 0.,   // x coordinates
	             0., 0., 1., 1.,  0., 0., 1., 1.,   // y coordinates
	             0., 0., 0., 0.,  1., 1., 1., 1.);  // z coordinates
	b.SetFTags  (6, -100,-101,-102,-103,-104,-105); // face tags
	b.SetNx     (3);                                // num of divisions along x
	b.SetNy     (3);                                // num of divisions along y
	b.SetNz     (3);                                // num of divisions along z
	Array<Mesh::Block*> blocks;
	blocks.Push (&b);

	// Generate
	cout << "\nMesh Generation: --------------------------------------------------------------" << endl;
	Mesh::Structured ms(/*Is3D*/true);
	clock_t start = std::clock(); // Initial time
	size_t  ne    = ms.Generate (blocks);
	clock_t total = std::clock() - start; // Time elapsed
	cout << "[1;33m"<<ne<<" elements[0m. Time elapsed = [1;31m" << static_cast<double>(total)/CLOCKS_PER_SEC << "[0m [1;32mseconds[0m" << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Geom g(3); // 3D

	// Solver
	FEM::Solver * sol = FEM::AllocSolver("AutoME");
	sol -> SetGeom(&g) -> SetLinSol(linsol.CStr()) -> SetNumDiv(1) -> SetDeltaTime(0.0);

	// Element attributes
	String prms; prms.Printf("lam=%f kap=%f phics=%f G=%f",lam,kap,phics,G);
	String inis; inis.Printf("Sx=%f Sy=%f Sz=%f Sxy=0 Syz=0 Szx=0 v=%f",p_ini,p_ini,p_ini,v_ini);
	FEM::EAtts_T eatts;
	eatts.Push (make_tuple(-1, "Hex8Equilib", "CamClay", prms.CStr(), inis.CStr(), "")); // tag, type, model, prms, inis, props

	// Set Nodes and Elements
	FEM::SetNodesElems (&ms, &eatts, &g);

	// Open file
	std::ofstream res("tcamclay01.cal", std::ios::out);
	res << _8s<<"Sx" << _8s<<"Sy" << _8s<<"Sz" << _8s<<"Ex" << _8s<<"Ey" << _8s<<"Ez" << _8s<<"p" << _8s<<"q" << _8s<<"Ev" << _8s<<"Ed" << "\n";

	// Solve each stage -- stress path
	Array<Tensor1> dtrac;  dtrac.Resize(2); // stress path (delta traction)
	dtrac[0] = 0.0, 0.0, -2.0;              // dtracx, dtracy, dtracz
	dtrac[1] = 0.0, 0.0,  3.0;
	for (size_t i=0; i<dtrac.Size(); ++i)
	{
		size_t ndiv = 10;
		for (size_t j=0; j<ndiv; ++j)
		{
			// Faces brys
			FEM::FBrys_T fbrys;
			fbrys.Push (make_tuple(-100, "ux", 0.0));
			fbrys.Push (make_tuple(-101, "fx", dtrac[i](0)/ndiv));
			fbrys.Push (make_tuple(-102, "uy", 0.0));
			fbrys.Push (make_tuple(-103, "fy", dtrac[i](1)/ndiv));
			fbrys.Push (make_tuple(-104, "uz", 0.0));
			fbrys.Push (make_tuple(-105, "fz", dtrac[i](2)/ndiv));

			// Set boundary conditions
			FEM::SetBrys (&ms, NULL, NULL, &fbrys, &g);

			// Solve
			sol -> Solve();

			// Output
			res << _8s<<g.Ele(13)->Val("Sx") << _8s<<g.Ele(13)->Val("Sy") << _8s<<g.Ele(13)->Val("Sz");
			res << _8s<<g.Ele(13)->Val("Ex") << _8s<<g.Ele(13)->Val("Ey") << _8s<<g.Ele(13)->Val("Ez");
			res << _8s<<g.Ele(13)->Val("p")  << _8s<<g.Ele(13)->Val("q")  << _8s<<g.Ele(13)->Val("Ev") << _8s<<g.Ele(13)->Val("Ed") << "\n";
		}
	}
	delete sol;

	// Close file
	res.close();
	cout << "[1;34mFile <tcamclay01.cal> saved.[0m" << endl;

	// Output final state
	Output o; o.VTU (&g, "tcamclay01.vtu");
	cout << "[1;34mFile <tcamclay01.vtu> saved.[0m" << endl;
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
