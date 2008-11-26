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

/* Unit cube with a uniform pressure at the top
 
      z
      |__y      +________________+
   x,'        ,'|              ,'|
            ,'               ,'  |
          ,'    |          ,'    |
        ,'      .        ,'      | 1.0
      +'_______________+'        |
      |                |         |
      |         |      |         |
      |         + -  - | -  -  - +
      |       ,        |       ,' 
      |     ,          |     ,'   
      |   ,            |   ,'  1.0
      | ,      1.0     | ,'       
      +________________+'         
*/

// STL
#include <iostream>

// MechSys
#include "fem/geometry.h"
#include "fem/functions.h"
#include "fem/output.h"
#include "fem/elems/hex8equilib.h"
#include "fem/elems/hex20equilib.h"
#include "fem/solvers/autome.h"
#include "fem/solvers/forwardeuler.h"
#include "models/equilibs/linelastic.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using boost::make_tuple;
using Util::_4;
using Util::_8s;

int main(int argc, char **argv) try
{
	// constants
	double E     = 200.0; // Young
	double nu    = 0.25;  // Poisson
	double q     = -2.0;  // Downward vertical pressure
	int    ndiv  = 4;     // number of divisions along x, y, and z
	bool   is_o2 = false; // use high order elements?
	String linsol("UM");  // UMFPACK

	// Input
	cout << "Input: " << argv[0] << "  is_o2  ndiv  linsol(LA,UM,SLU)\n";
	if (argc>=2) is_o2      = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) ndiv       =  atof(argv[2]);
	if (argc>=4) linsol.Printf("%s",argv[3]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
	
	Mesh::Block b;
	b.SetTag    (-1); // tag to be replicated to all generated elements inside this block
	b.SetCoords (true, 8,                          // Is3D, NNodes
	             0., 1., 1., 0.,  0., 1., 1., 0.,  // x coordinates
	             0., 0., 1., 1.,  0., 0., 1., 1.,  // y coordinates
	             0., 0., 0., 0.,  1., 1., 1., 1.); // z coordinates
	b.SetFTags  (6, -100,0,-102,0,-104,-105);      // face tags
	b.SetNx     (ndiv);                            // num of divisions along x
	b.SetNy     (ndiv);                            // num of divisions along y
	b.SetNz     (ndiv);                            // num of divisions along z
	Array<Mesh::Block*> blocks;
	blocks.Push (&b);

	// Generate
	Mesh::Structured mesh(/*Is3D*/true);
	if (is_o2) mesh.SetO2();                // Non-linear elements
	clock_t start = std::clock();           // Initial time
	size_t  ne    = mesh.Generate (blocks); // Discretize domain
	clock_t total = std::clock() - start;   // Time elapsed
	if (is_o2) cout << "\nNumber of hexahedra(o2) = " << ne << endl;
	else       cout << "\nNumber of hexahedra     = " << ne << endl;
	cout << "Time elapsed (mesh)     = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Geom g(3); // 3D

	// Faces brys
	FEM::FBrys_T fbrys;
	fbrys.Push (make_tuple(-100, "ux", 0.0)); // tag, key, val
	fbrys.Push (make_tuple(-102, "uy", 0.0)); // tag, key, val
	fbrys.Push (make_tuple(-104, "uz", 0.0)); // tag, key, val
	fbrys.Push (make_tuple(-105, "fz",   q)); // tag, key, val

	// Element attributes
	String prms; prms.Printf("E=%f nu=%f",E,nu);
	FEM::EAtts_T eatts;
	if (is_o2) eatts.Push (make_tuple(-1, "Hex20Equilib", "LinElastic", prms.CStr(), "ZERO", "")); // tag, type, model, prms, inis, props
	else       eatts.Push (make_tuple(-1, "Hex8Equilib",  "LinElastic", prms.CStr(), "ZERO", "")); // tag, type, model, prms, inis, props

	// Set geometry
	FEM::SetNodesElems (&mesh, &eatts, &g);
	FEM::SetBrys       (&mesh, NULL, NULL, &fbrys, &g);

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol -> SetGeom(&g) -> SetLinSol(linsol.CStr()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	start = std::clock();
	sol -> Solve();
	total = std::clock() - start;
	double norm_resid = LinAlg::Norm(sol->Resid());
	cout << "Time elapsed (solution) = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";
	cout << "[1;35mNorm(Resid=DFext-DFint) = " << norm_resid << "[0m\n";
	cout << "[1;32mNumber of DOFs          = " << sol->nDOF() << "[0m\n";
	delete sol;

	// Output
	Output o; o.VTU (&g, "tsolid01.vtu");
	cout << "[1;34mFile <tsolid01.vtu> saved.[0m\n\n";

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Check
    Array<double> err_eps;
    Array<double> err_sig;
    Array<double> err_dis;

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
	for (size_t i=0; i<g.NElems(); ++i)
	{
		for (size_t j=0; j<g.Ele(i)->NNodes(); ++j)
		{
			// eps
			err_eps.Push( fabs(g.Ele(i)->Val(j,"Ex" ) - Ex ) / (1.0+fabs(Ex )) );
			err_eps.Push( fabs(g.Ele(i)->Val(j,"Ey" ) - Ey ) / (1.0+fabs(Ey )) );
			err_eps.Push( fabs(g.Ele(i)->Val(j,"Ez" ) - Ez ) / (1.0+fabs(Ez )) );
			err_eps.Push( fabs(g.Ele(i)->Val(j,"Exy") - Exy) / (1.0+fabs(Exy)) );
			err_eps.Push( fabs(g.Ele(i)->Val(j,"Eyz") - Eyz) / (1.0+fabs(Eyz)) );
			err_eps.Push( fabs(g.Ele(i)->Val(j,"Ezx") - Ezx) / (1.0+fabs(Ezx)) );
			// sig
			err_sig.Push( fabs(g.Ele(i)->Val(j,"Sx" ) - Sx ) / (1.0+fabs(Sx )) );
			err_sig.Push( fabs(g.Ele(i)->Val(j,"Sy" ) - Sy ) / (1.0+fabs(Sy )) );
			err_sig.Push( fabs(g.Ele(i)->Val(j,"Sz" ) - Sz ) / (1.0+fabs(Sz )) );
			err_sig.Push( fabs(g.Ele(i)->Val(j,"Sxy") - Sxy) / (1.0+fabs(Sxy)) );
			err_sig.Push( fabs(g.Ele(i)->Val(j,"Syz") - Syz) / (1.0+fabs(Syz)) );
			err_sig.Push( fabs(g.Ele(i)->Val(j,"Szx") - Szx) / (1.0+fabs(Szx)) );
		}
	}

	// Displacements
	for (size_t i=0; i<g.NNodes(); ++i)
	{
		double ux_correct = Ex*g.Nod(i)->X();
		double uy_correct = Ey*g.Nod(i)->Y();
		double uz_correct = Ez*g.Nod(i)->Z();
		err_dis.Push ( fabs(g.Nod(i)->Val("ux") - ux_correct) / (1.0+fabs(ux_correct)) );
		err_dis.Push ( fabs(g.Nod(i)->Val("uy") - uy_correct) / (1.0+fabs(uy_correct)) );
		err_dis.Push ( fabs(g.Nod(i)->Val("uz") - uz_correct) / (1.0+fabs(uz_correct)) );
	}

	// Error summary
	double tol_eps     = 1.0e-16;
	double tol_sig     = 1.0e-14;
	double tol_dis     = 1.0e-16;
	double min_err_eps = err_eps[err_eps.Min()];
	double min_err_sig = err_sig[err_sig.Min()];
	double min_err_dis = err_dis[err_dis.Min()];
	double max_err_eps = err_eps[err_eps.Max()];
	double max_err_sig = err_sig[err_sig.Max()];
	double max_err_dis = err_dis[err_dis.Max()];
	cout << _4<< ""    << _8s<<"Min"       << _8s<<"Mean"                                                        << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	cout << _4<< "Eps" << _8s<<min_err_eps << _8s<<err_eps.Mean() << (max_err_eps>tol_eps?"[1;31m":"[1;32m") << _8s<<max_err_eps << "[0m" << _8s<<err_eps.Norm() << endl;
	cout << _4<< "Sig" << _8s<<min_err_sig << _8s<<err_sig.Mean() << (max_err_sig>tol_sig?"[1;31m":"[1;32m") << _8s<<max_err_sig << "[0m" << _8s<<err_sig.Norm() << endl;
	cout << _4<< "Dis" << _8s<<min_err_dis << _8s<<err_dis.Mean() << (max_err_dis>tol_dis?"[1;31m":"[1;32m") << _8s<<max_err_dis << "[0m" << _8s<<err_dis.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err_eps>tol_eps || max_err_sig>tol_sig || max_err_dis>tol_dis) return 1;
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
