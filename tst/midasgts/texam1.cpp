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
#include "fem/elems/beam.h"
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/output.h"
#include "util/exception.h"
#include "linalg/matrix.h"
#include "linalg/laexpr.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_8s;
using Util::PI;
using boost::make_tuple;

// Calculate Kirsch's solution for a cylindrical hole (stresses)
inline void KirschStress(double ph, double pv, double r, double R, double T, double & SigR, double & SigT, double & SigRT)
{
	double pm = (ph+pv)/2.0;
	double pd = (ph-pv)/2.0;
	double c1 = r*r/(R*R);
	SigR  =  pm*(1.0-c1) + pd*(1.0-4.0*c1+3.0*c1*c1)*cos(2.0*T);
	SigT  =  pm*(1.0+c1) - pd*(1.0+3.0*c1*c1)*cos(2.0*T);
	SigRT = -pd*(1.0+2.0*c1-3.0*c1*c1)*sin(2.0*T);
}

// Calculate Kirsch's solution for a cylindrical hole (displacements)
void KirschDisp(double ph, double pv, double r, double R, double T, double E, double nu, double & uR, double & uT)
{
	double G  = E/(2.0*(1.0+nu)); // Shear modulus
	double c1 = r*r/R;
	double qm = (ph+pv)/(4.0*G);
	double qd = (ph-pv)/(4.0*G);
	uR =  qm*c1 + qd*c1*(4.0*(1.0-nu)-c1/R)*cos(2.0*T);
	uT = -qd*c1*(2.0*(1.0-2.0*nu)+c1/R)*sin(2.0*T);
}

int main(int argc, char **argv) try
{
	// Constants
	double ph      = -30.0;    // horizontal pressure
	double pv      = -30.0;    // vertical pressure
	double E_soil  = 6000.0;   // Young [MPa]
	double nu_soil = 0.2;      // Poisson
	double r       = 1.0;      // radius
	double L       = 10.0;     // length
	double H       = 10.0;     // height
	bool   is_o2   = false;    // use high order elements?
	int    ndivy   = 15;       // ndivy
	double Ax      = 2.0;      // rate of increase of X divisions
	double NonLinX = false;    // nonlinear divisions along X?

	// Input
	cout << "Input: " << argv[0] << "  is_o2  ndivy(ndivx=2*ndivy)\n";
	if (argc>=2) is_o2 = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) ndivy =  atoi(argv[2]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
	
	/*                     | | | pv | | |
	 *                     V V V    V V V
	 *
	 *                |---------- L ----------|
	 *   -+- -+- -+-  o___________o_-30_______o
	 *    |   |   |   |                      ,|
	 *    |   |   |   |-40    [b1]         ,' |
	 *    |   |   d   o                  ,'   |
	 *    |   f   |   |    y       x   ,'     |      <--
	 *    |   |   |   |     ',   ,'  ,o       |      <--
	 *    |   |  -+-  o-,_    '+'  ,'         |      <--
	 *    H   |           o-,    ,'           o -20  <-- ph
	 *    |  -+- . . . . . . 'o '      [b0]   |      <--
	 *    |   |               .',             |      <--
	 *    |   e               .  o  y^        |      <--
	 *    |   |               .   \  |        |
	 *    |   |               .   |  +-->x    |
	 *   -+- -+-      +----r----->o-----o-----o
	 *                        .       -10
	 *                        .   |---- a ----|
	 *                |-- b --|------ c ------|
	 */

	// Geometry
	double a = L-r;
	double b = r*cos(2.*PI/8.);
	double c = L-b;
	double d = H-r;
	double e = r*sin(2.*PI/8.);
	double f = H-e;

	// Lower block -- coordinates
	Mesh::Block b0;
	b0.SetTag    (-1);
	b0.SetCoords (false, 8, // Is3D, NNodes
	               r,  L, L, b,    r+a/2.,    L, b+c/2., r*cos(PI/8.),
	              0., 0., H, e,        0., H/2., e+f/2., r*sin(PI/8.));
	b0.SetNx     (2*ndivy, /*Ax*/Ax, /*Nonlinear*/NonLinX);
	b0.SetNy     (ndivy);
	b0.SetETags  (4, 0, -20, -10, 0);

	// Upper block -- coordinates
	Mesh::Block b1;
	b1.SetTag    (-1);
	b1.SetCoords (false, 8,
	              b, L, 0., 0.,   b+c/2., L/2.,     0., r*cos(3.*PI/8.),
	              e, H,  H,  r,   e+f/2.,    H, r+d/2., r*sin(3.*PI/8.));
	b1.SetNx     (2*ndivy, /*Ax*/Ax, /*Nonlinear*/NonLinX);
	b1.SetNy     (ndivy);
	b1.SetETags  (4, 0, -30,  0, -40);

	// Blocks
	Array<Mesh::Block*> blocks;  blocks.Resize(2);
	blocks[0] = &b0;
	blocks[1] = &b1;

	// Generate
	Mesh::Structured mesh(/*Is3D*/false);
	if (is_o2) mesh.SetO2();               // Non-linear elements
	clock_t start = std::clock();          // Initial time
	size_t  ne    = mesh.Generate(blocks); // Discretize domain
	clock_t total = std::clock() - start;  // Time elapsed
	cout << "\nNumber of quads     = " << ne << endl;
	cout << "Time elapsed (mesh) = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Geom g(2); // 2D

	// Edges brys
	FEM::EBrys_T ebrys;
	ebrys.Push (make_tuple(-10, "uy", 0.));
	ebrys.Push (make_tuple(-20, "fx", ph));
	ebrys.Push (make_tuple(-30, "fy", pv));
	ebrys.Push (make_tuple(-40, "ux", 0.));

	// Elements attributes
	String prms; prms.Printf("E=%f nu=%f",E_soil,nu_soil);
	FEM::EAtts_T eatts;
	if (is_o2) eatts.Push (make_tuple(-1, "Quad8PStrain", "LinElastic", prms.CStr(), "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0"));
	else       eatts.Push (make_tuple(-1, "Quad4PStrain", "LinElastic", prms.CStr(), "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0"));

	// Set geometry: nodes, elements, attributes, and boundaries
	FEM::SetNodesElems (&mesh, &eatts, &g, 1.0e-5);
	FEM::SetBrys       (&mesh, NULL, &ebrys, NULL, &g);

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol -> SetGeom(&g) -> SetLinSol("UM") -> SetNumDiv(1) -> SetDeltaTime(0.0);
	start = std::clock();
	sol -> Solve();
	total = std::clock() - start;
	double norm_resid = LinAlg::Norm(sol->Resid());
	cout << "Time elapsed (solution) = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";
	cout << "[1;35mNorm(Resid=DFext-DFint) = " << norm_resid << "[0m\n";
	cout << "[1;32mNumber of DOFs          = " << sol->nDOF() << "[0m\n";

	// Output: VTU
	start = std::clock();
	Output o; o.VTU (&g, "texam1.vtu");
	total = std::clock() - start;
	cout << "Time elapsed (output file) = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";
	cout << "[1;34mFile <texam1.vtu> saved.[0m\n\n";

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	start = std::clock();

	// Stress
	Array <double> err_sR;
	Array <double> err_sT;
	Array <double> err_sRT;
	for (size_t i=0; i<g.NElems(); ++i)
	{
		for (size_t j=0; j<g.Ele(i)->NNodes(); ++j)
		{
			// Analytical
			double x = g.Ele(i)->Nod(j)->X();
			double y = g.Ele(i)->Nod(j)->Y();
			double t = atan(y/x);
			double R = sqrt(x*x+y*y);
			double sigRc, sigTc, sigRTc; // correct stress components
			KirschStress (ph,pv,r,R,t, sigRc,sigTc,sigRTc);

			// Numerical
			double c     = x/R;
			double s     = y/R;
			double cc    = c*c;
			double ss    = s*s;
			double sc    = s*c;
			double Sx    = g.Ele(i)->Val(j,"Sx");
			double Sy    = g.Ele(i)->Val(j,"Sy");
			double Sxy   = g.Ele(i)->Val(j,"Sxy");
			double sigR  = Sx*cc + Sy*ss + 2.0*Sxy*sc;
			double sigT  = Sx*ss + Sy*cc - 2.0*Sxy*sc;
			double sigRT = (Sy-Sx)*sc + (cc-ss)*Sxy;

			// Error
			err_sR .Push (fabs(sigR  - sigRc ));
			err_sT .Push (fabs(sigT  - sigTc ));
			err_sRT.Push (fabs(sigRT - sigRTc));
		}
	}

	// Displacements
	Array <double> err_uR;
	Array <double> err_uT;
	for (size_t i=0; i<g.NNodes(); ++i)
	{
		// Analytical
		double x  = g.Nod(i)->X();
		double y  = g.Nod(i)->Y();
		double t  = atan(y/x);
		double R  = sqrt(x*x+y*y);
		double uRc, uTc; // correct displacement components
		KirschDisp (ph,pv,r,R,t, E_soil,nu_soil, uRc,uTc);

		// Numerical
		double c  = x/R;
		double s  = y/R;
		double ux = g.Nod(i)->Val("ux");
		double uy = g.Nod(i)->Val("uy");
		double uR = ux*c + uy*s;
		double uT = uy*c - ux*s;

		// Error
		err_uR.Push (fabs(uR-uRc));
		err_uT.Push (fabs(uT-uTc));
	}

	// Error summary
	double tol_sR      = 3.0e0;
	double tol_sT      = 3.0e0;
	double tol_sRT     = 3.0e0;
	double tol_uR      = 1.0e-1;
	double tol_uT      = 1.0e-3;
	double min_err_sR  = err_sR [err_sR .Min()];   double max_err_sR  = err_sR [err_sR .Max()];
	double min_err_sT  = err_sT [err_sT .Min()];   double max_err_sT  = err_sT [err_sT .Max()];
	double min_err_sRT = err_sRT[err_sRT.Min()];   double max_err_sRT = err_sRT[err_sRT.Max()];
	double min_err_uR  = err_uR [err_uR .Min()];   double max_err_uR  = err_uR [err_uR .Max()];
	double min_err_uT  = err_uT [err_uT .Min()];   double max_err_uT  = err_uT [err_uT .Max()];
	cout << _4<< ""    << _8s<<"Min"       << _8s<<"Mean"                                                        << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	cout << _4<< "sR"  << _8s<<min_err_sR  << _8s<<err_sR .Mean() << (max_err_sR >tol_sR ?"[1;31m":"[1;32m") << _8s<<max_err_sR  << "[0m" << _8s<<err_sR.Norm()  << endl;
	cout << _4<< "sT"  << _8s<<min_err_sT  << _8s<<err_sT .Mean() << (max_err_sT >tol_sT ?"[1;31m":"[1;32m") << _8s<<max_err_sT  << "[0m" << _8s<<err_sT.Norm()  << endl;
	cout << _4<< "sRT" << _8s<<min_err_sRT << _8s<<err_sRT.Mean() << (max_err_sRT>tol_sRT?"[1;31m":"[1;32m") << _8s<<max_err_sRT << "[0m" << _8s<<err_sRT.Norm() << endl;
	cout << _4<< "uR"  << _8s<<min_err_uR  << _8s<<err_uR .Mean() << (max_err_uR >tol_uR ?"[1;31m":"[1;32m") << _8s<<max_err_uR  << "[0m" << _8s<<err_uR.Norm()  << endl;
	cout << _4<< "uT"  << _8s<<min_err_uT  << _8s<<err_uT .Mean() << (max_err_uT >tol_uT ?"[1;31m":"[1;32m") << _8s<<max_err_uT  << "[0m" << _8s<<err_uT.Norm()  << endl;

	total = std::clock() - start;
	cout << "Time elapsed (error check) = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";

	// Return error flag
	if (max_err_sR>tol_sR || max_err_sT>tol_sT || max_err_sRT>tol_sRT || max_err_uR>tol_uR || max_err_uT>tol_uT) return 1;
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
