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
#include "mesh/structured.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_8s;
using Util::PI;
using boost::make_tuple;

void Einsten_Schwartz(
		double p0,   // initial stress load
		double K,    // ratio beetween horizontal and vertical initial stress load
		double R,    // radius: polar coordinate
		double th,   // angle:  polar coordinate
		double r,    // radius of the hole
		double E,    // rock Young modulus
		double nu,   // rock Poisson ratio
		double E_s,  // liner Young modulus
		double nu_s, // liner Poisson ratio
		double I,    // liner moment of inertia
		double A,    // cross sectional area of the liner for a 1m long section
		double & M,  // bending moment in the liner
		double & N,  // axial force in the liner
		double & uT, // tangential displacement
		double & uR  // radial displacement
		)
{
	double Nu   = 1.0-nu;
	double Nu2  = 1.0-nu*nu;
	double Nus2 = 1.0-nu_s*nu_s;
	double C    = E*R*Nus2 / (E_s*A*Nu2);
	double F    = E*pow(R,3)*Nus2 / ( E_s*I*(1.0-nu*nu) );
	double beta = ((6.0+F)*C*Nu + 2*F*nu) / (3*F + 3*C + 2*C*F*Nu);
	double b2   = C*Nu * 0.5 / (C*Nu + 4*nu - 6*beta - 3*beta*C*Nu);
	double a2   = beta*b2;
	double a0   = C*F*Nu / (C+F+C*F*Nu);
	double G    = E/2.0/(1+nu);

	M = p0*r*r*0.25*(1.0-K)*(1.0 - 2.0*a2 + 2.0*b2)*cos(2.0*th);
	N = p0*r*(0.5*(1.0+K)*(1.0 - a0) + 0.5*(1.0-K)*(1.0+2*a2)*cos(2.0*th) );
	uT = 0.5*p0*r/G*(1.0-K)*(a2 + (1.0 - 2*nu)*b2)*sin(2.0*th);
	uR = 0.5*p0*r/G*(-0.5*(1.0+K)*a0 + 0.5*(1.0-K)*( 4.0*(1.0-nu)*b2 - 2.0*a2 )*cos(2*th) );
}

int main(int argc, char **argv) try
{
	// Constants
	double E_soil   = 6000.0;   // Young [MPa]
	double nu_soil  = 0.2;      // Poisson
	double E_beam   = 20000.0;  // Young [MPa]
	double Izz_beam = 0.01042;  // Inertia m^4
	double A_beam   = 0.5;      // Area m^2
	double r        = 1.0;      // radius
	double L        = 10.0;     // length
	double H        = 10.0;     // height
	bool   is_o2    = false;    // use high order elements?

	// Input
	cout << "Input: " << argv[0] << "  is_o2\n";
	if (argc>=2) is_o2 = (atoi(argv[1])>0 ? true : false);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
	
	/*                |---------- L ----------|
	 *   -+- -+- -+-  o___________o_-30_______o
	 *    |   |   |   |                      ,|
	 *    |   |   |   |-40    [b1]         ,' |
	 *    |   |   d   o                  ,'   |
	 *    |   f   |   |    y       x   ,'     |
	 *    |   |   |   |     ',   ,'  ,o       |
	 *    |   |  -+-  o-,_    '+'  ,'         |
	 *    H   |       -55 o-,    ,'           o -10
	 *    |  -+- . . . . . . 'o '      [b0]   |
	 *    |   |               .',             |
	 *    |   e               .  o  y^        |
	 *    |   |               .-55\  |        |
	 *    |   |               .   |  +-->x    |
	 *   -+- -+-      +----r----->o-----o-----o
	 *                        .       -20
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
	b0.SetNx     (4, 2.0, true);
	b0.SetNy     (4);
	b0.SetETags  (4, -55, -10, -20, 0);

	// Upper block -- coordinates
	Mesh::Block b1;
	b1.SetTag    (-1);
	b1.SetCoords (false, 8,
	              b, L, 0., 0.,   b+c/2., L/2.,     0., r*cos(3.*PI/8.),
	              e, H,  H,  r,   e+f/2.,    H, r+d/2., r*sin(3.*PI/8.));
	b1.SetNx     (4, 2.0, true);
	b1.SetNy     (4);
	b1.SetETags  (4, -55, -30,  0, -40);

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

	// Nodes brys
	FEM::NBrys_T nbrys;
	nbrys.Push (make_tuple(r,0.0,0.0,"wz",0.0));
	nbrys.Push (make_tuple(0.0,r,0.0,"wz",0.0));

	// Edges brys
	FEM::EBrys_T ebrys;
	double p0 = -15.0;
	double  K =  2.0;
	ebrys.Push (make_tuple(-10, "fx", p0*K));
	ebrys.Push (make_tuple(-20, "uy",  0.0));
	ebrys.Push (make_tuple(-30, "fy",   p0));
	ebrys.Push (make_tuple(-40, "ux",  0.0));

	// Elements attributes
	String prms; prms.Printf("E=%f nu=%f",E_soil,nu_soil);
	FEM::EAtts_T eatts;
	if (is_o2) eatts.Push (make_tuple(-1, "Quad8PStrain", "LinElastic", prms.CStr(), "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0"));
	else       eatts.Push (make_tuple(-1, "Quad4PStrain", "LinElastic", prms.CStr(), "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0"));
	//eatts.Push (make_tuple(-55, "Quad4PStrain", "LinElastic", prms.CStr(), "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0"));

	// Beam attributes
	String beamprms;
	beamprms.Printf("E=%f A=%f Izz=%f",E_beam,A_beam,Izz_beam);
	eatts.Push (make_tuple(-55, "Beam", "", beamprms.CStr(), "ZERO"));

	// Set geometry: nodes, elements, attributes, and boundaries
	FEM::SetNodesElems (&mesh, &eatts, &g, 1.0e-5);
	FEM::SetBrys       (&mesh, &nbrys, &ebrys, NULL, &g);

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

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

    Vector<double> sig(3);
    Vector<double> sig_k(3);
    Array <double> err_sig_x;
    Array <double> err_sig_y;
    Array <double> err_sig_xy;
    Array <double> err_disp;

	// Stress 
	//for (size_t i=0; i<g.NElems(); ++i)
	//{
	//	for (size_t j=0; j<g.Ele(i)->NNodes(); ++j)
	//	{
	//		// Analytical
	//		double x  = g.Ele(i)->Nod(j)->X();
	//		double y  = g.Ele(i)->Nod(j)->Y();
	//		double th = atan(y/x);
	//		double R  = sqrt(x*x + y*y);
	//		//Kirsch_stress(p1, p2, r, R, th, sig_k); // Calculate the Kirsch solution for a cylindrical hole



	//		// Analysis
	//		sig = g.Ele(i)->Val(j,"Sx"), g.Ele(i)->Val(j,"Sy"), g.Ele(i)->Val(j,"Sxy");
	//		Vector<double> sig_polar(3); // Vector to transform to polar coordinates
	//		sig_polar = sig(0)*pow(cos(th),2.0) + sig(1)*pow(sin(th),2.0) + 2.0*sig(2)*sin(th)*cos(th),
	//		            sig(0)*pow(sin(th),2.0) + sig(1)*pow(cos(th),2.0) - 2.0*sig(2)*sin(th)*cos(th),
	//		            (sig(1)-sig(0))*sin(th)*cos(th) + sig(2)*(pow(cos(th),2.0) - pow(sin(th),2.0));

	//		err_sig_x .Push ( fabs(sig_k(0) - sig_polar(0)) );
	//		err_sig_y .Push ( fabs(sig_k(1) - sig_polar(1)) );
	//		err_sig_xy.Push ( fabs(sig_k(2) - sig_polar(2)) );
	//	}
	//}

	// Displacements
	for (size_t i=0; i<g.NNodes(); ++i)
	{
		// Analytical
		double x  = g.Nod(i)->X();
		double y  = g.Nod(i)->Y();
		double th = atan(y/x);
		double R  = sqrt(x*x + y*y);
		double uR_correct;
		double uT_correct;
		double nu_beam = 0.25;
		double M, N;

		Einsten_Schwartz( p0, K, R, th, r, E_soil, nu_soil, E_beam, nu_beam, Izz_beam, A_beam, M, N, uT_correct, uR_correct);

		// Analysis
		double ux = g.Nod(i)->Val("ux");
		double uy = g.Nod(i)->Val("uy");
		double uR = ux*cos(th) + uy*sin(th);
		double uT = uy*cos(th) - ux*sin(th);

		err_disp.Push ( uR - uR_correct) ;
		err_disp.Push ( uT - uT_correct) ;
	}
	// Error summary
	double tol_disp       = 1.0e-1;
	double min_err_disp   = err_disp  [err_disp.Min()];
	double max_err_disp   = err_disp  [err_disp.Max()];
	cout << _4 << ""      << _8s <<"Min"       << _8s << "Mean"                                                        << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	cout << _4 << "Disp"  << _8s <<min_err_disp << _8s<<err_disp  .Mean() << (max_err_disp  >tol_disp?"[1;31m":"[1;32m") << _8s<<max_err_disp   << "[0m" << _8s<<err_disp  .Norm() << endl;
	cout << endl;




	// Output: VTU
	Output o; o.VTU (&g, "texam5.vtu");
	cout << "[1;34mFile <texam5.vtu> saved.[0m\n\n";

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
