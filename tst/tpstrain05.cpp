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
#include "fem/elems/quad8pstrain.h"
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "fem/output.h"
#include "util/exception.h"
#include "linalg/matrix.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_8s;
using Util::PI;
using boost::make_tuple;

void Analyitical(double r, double R, double p0, double p1, double E, double nu, double L, double & SigR, double & SigT, double & UR)
{
	double c1 = 0.5*(r*r*p1-R*R*p0)/(R*R-r*r);
	double c2 = r*r*R*R*(p0-p1)/(R*R-r*r);
	SigR = -(2.0*c1 + c2/(L*L));
	SigT = -(2.0*c1 - c2/(L*L));
	double _E  =  E/(1.0-nu*nu);
	double _nu = nu/(1.0-nu);
	UR = -(1.0/_E*(2*c1*(1-_nu)*L-c2*(1+_nu)/L));
}

int main(int argc, char **argv) try
{
	// Constants
	double E     = 207.0; // Young
	double nu    = 0.3;   // Poisson
	double p0    = -2.0;  // Load
	double p1    = -1.0;  // Load
	int    ndiv  = 30;    // number of divisions along x and y
	bool   is_o2 = false; // use high order elements?
	String linsol("UM");  // UMFPACK

	// Input
	cout << "Input: " << argv[0] << "  is_o2  ndiv  linsol(LA,UM,SLU)\n";
	if (argc>=2) is_o2      = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) ndiv       =  atof(argv[2]);
	if (argc>=4) linsol.Printf("%s",argv[3]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
	
	/*            |---------- R ----------|
	 *            @--,__ 
	 *            |     '--,__    
	 *        -40 |           '--,   -20
	 *            @               '@
	 *            |                 ',
	 *            |                .  \
	 *            @-,_                 ',
	 *                '-,          .     \
	 *                   '@               \
	 *               -10  .',      .      | 
	 *                    .  '  y^        |
	 *                    .   \  | .      |
	 *                    .   |  +-->x    |
	 *            +----r----->@-----@-----@
	 *                    .   .    .  -30
	 *                    .   |---- a ----|
	 *            |-- b --|        .
	 *            |------- c ------|
	 */

	// Geometry
	double R = 10.0;
	double r = 6.0;
	double a = R-r;
	double b = r*cos(PI/4.);
	double c = R*cos(PI/4.);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Blocks
	Mesh::Block bl;
	bl.SetTag    (-1); // tag to be replicated to all generated elements inside this block
	bl.SetCoords (false, 8,                              // Is3D, NNodes
	               r,  R, 0., 0., r+a/2., c,     0., b,  // x coordinates
	              0., 0.,  R,  r,     0., c, r+a/2., b); // y coordinates
	bl.SetNx      (ndiv,/*Ax*/1.,/*Nonlinear*/false);    // x weights and num of divisions along x
	bl.SetNy      (ndiv);                                // y weights and num of divisions along y
	bl.SetETags   (4, -10, -20, -30, -40);               // edge tags
	Array<Mesh::Block*> blocks;
	blocks.Push (&bl);

	// Generate
	Mesh::Structured mesh(/*Is3D*/false);
	if (is_o2) mesh.SetO2();                // Non-linear elements
	clock_t start = std::clock();           // Initial time
	size_t  ne    = mesh.Generate (blocks); // Discretize domain
	clock_t total = std::clock() - start;   // Time elapsed
	if (is_o2) cout << "\nNum of quadrangles (o2) = " << ne << endl;
	else       cout << "\nNumber of quadrangles   = " << ne << endl;
	cout << "Time elapsed (mesh)     = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Geom g(2); // 2D

	// Edges brys
	FEM::EBrys_T ebrys;
	ebrys.Push (make_tuple(-30, "uy", 0.0));
	ebrys.Push (make_tuple(-40, "ux", 0.0));
	ebrys.Push (make_tuple(-20, "Q",   p0));
	ebrys.Push (make_tuple(-10, "Q",   p1));

	// Elements attributes
	String prms; prms.Printf("E=%f nu=%f",E,nu);
	FEM::EAtts_T eatts;
	if (is_o2) eatts.Push (make_tuple(-1, "Quad8PStrain", "LinElastic", prms.CStr(), "ZERO", "", true));
	else       eatts.Push (make_tuple(-1, "Quad4PStrain", "LinElastic", prms.CStr(), "ZERO", "", true));

	// Set geometry: nodes, elements, attributes, and boundaries
	FEM::SetNodesElems (&mesh, &eatts, &g);
	FEM::SetBrys       (&mesh, NULL, &ebrys, NULL, &g);

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol->SetGeom(&g)->SetLinSol(linsol.CStr());
	sol->SolveWithInfo(/*NDiv*/1, /*DTime*/0.0);
	delete sol;

	// Output: VTU
	Output o; o.VTU (&g, "tpstrain05.vtu");
	cout << "[1;34mFile <tpstrain05.vtu> saved.[0m\n\n";

	//////////////////////////////////////////////////////////////////////////////////////// Check /////


	// Stress
	Array <double> err_sR;
	Array <double> err_sT;
	Array <double> err_uR;
	for (size_t i=0; i<g.NElems(); ++i)
	{
		for (size_t j=0; j<g.Ele(i)->NNodes(); ++j)
		{
			// Analytical
			double x = g.Ele(i)->Nod(j)->X();
			double y = g.Ele(i)->Nod(j)->Y();
			double L = sqrt(x*x+y*y);
			double sigRc, sigTc, uRc; // correct stress components
			Analyitical(r, R, p0, p1, E, nu, L, sigRc, sigTc, uRc);

			// Numerical
			double c     = x/L;
			double s     = y/L;
			double cc    = c*c;
			double ss    = s*s;
			double sc    = s*c;
			double Sx    = g.Ele(i)->Val(j,"Sx");
			double Sy    = g.Ele(i)->Val(j,"Sy");
			double Sxy   = g.Ele(i)->Val(j,"Sxy");
			double ux    = g.Ele(i)->Val(j,"ux");
			double uy    = g.Ele(i)->Val(j,"uy");
			double sigR  = Sx*cc + Sy*ss + 2.0*Sxy*sc;
			double sigT  = Sx*ss + Sy*cc - 2.0*Sxy*sc;
			double uR    = ux*c  + uy*s;

			// Error
			err_sR.Push (fabs(sigR - sigRc));
			err_sT.Push (fabs(sigT - sigTc));
			err_uR.Push (fabs(  uR - uRc  ));

			if (i==400 && j==0)
			{
				cout << "sigR = " << sigR << " sigRc = " << sigRc << endl;
				cout << "sigT = " << sigT << " sigTc = " << sigTc << endl;
				cout << "uR = "   << uR   << " uRc = "   << uRc   << endl;
			}
		}
	}

	// Error summary
	double tol_sR      = 3.0e0;
	double tol_sT      = 3.0e0;
	double tol_uR      = 1.0e-1;
	double min_err_sR  = err_sR [err_sR .Min()];   double max_err_sR  = err_sR [err_sR .Max()];
	double min_err_sT  = err_sT [err_sT .Min()];   double max_err_sT  = err_sT [err_sT .Max()];
	double min_err_uR  = err_uR [err_uR .Min()];   double max_err_uR  = err_uR [err_uR .Max()];
	cout << _4<< ""    << _8s<<"Min"       << _8s<<"Mean"                                                        << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	cout << _4<< "sR"  << _8s<<min_err_sR  << _8s<<err_sR .Mean() << (max_err_sR >tol_sR ?"[1;31m":"[1;32m") << _8s<<max_err_sR  << "[0m" << _8s<<err_sR.Norm()  << endl;
	cout << _4<< "sT"  << _8s<<min_err_sT  << _8s<<err_sT .Mean() << (max_err_sT >tol_sT ?"[1;31m":"[1;32m") << _8s<<max_err_sT  << "[0m" << _8s<<err_sT.Norm()  << endl;
	cout << _4<< "uR"  << _8s<<min_err_uR  << _8s<<err_uR .Mean() << (max_err_uR >tol_uR ?"[1;31m":"[1;32m") << _8s<<max_err_uR  << "[0m" << _8s<<err_uR.Norm()  << endl;

	total = std::clock() - start;
	cout << "Time elapsed (error check) = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";

	// Return error flag
	if (max_err_sR>tol_sR || max_err_sT>tol_sT || max_err_uR>tol_uR) return 1;
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
