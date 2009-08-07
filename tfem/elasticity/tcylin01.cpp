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
#include "fem/solver.h"
#include "fem/elems/quad4.h"
#include "fem/elems/quad8.h"
#include "fem/equilibelem.h"
#include "models/equilibs/linelastic.h"
#include "util/fatal.h"
#include "linalg/matrix.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_8s;
using Util::PI;

#define T boost::make_tuple

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
	 *            2--,__ 
	 *            |     '--,__    
	 *        -40 |           '--,   -20
	 *            6               '5
	 *            |                 ',
	 *            |                .  \
	 *            3-,_                 ',
	 *                '-,          .     \
	 *                   '7               \
	 *               -10  .',      .      | 
	 *                    .  '  y^        |
	 *                    .   \  | .      |
	 *                    .   |  +-->x    |
	 *            +----r----->0-----4-----1
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
	Array<Mesh::Block> bks(1);

	// Block # 0 --------------------------------
    Mesh::Verts_T ve0(8);
    Mesh::Edges_T ed0(8);
    Mesh::ETags_T et0(8);
	ve0 = T(0,  r      , 0.     , 0.),
	      T(1,  R      , 0.     , 0.),
	      T(2,  0.     , R      , 0.),
	      T(3,  0.     , r      , 0.),
	      T(4,  r+a/2. , 0.     , 0.),
	      T(5,  c      , c      , 0.),
	      T(6,  0.     , r+a/2. , 0.),
	      T(7,  b      , b      , 0.);
    ed0 = T(0,4), T(4,1), T(1,5), T(5,2), T(2,6), T(6,3), T(3,7), T(7,0);
    et0 = T(0,4,-30), T(4,1,-30), T(1,5,-20), T(5,2,-20), T(2,6,-40), T(6,3,-40), T(3,7,-10), T(7,0,-10);
    bks[0].Set   (-1, ve0, ed0, &et0, NULL, /*orig*/0, /*xplus*/4, /*yplus*/7);
	bks[0].SetNx (ndiv, /*Ax*/1.0, /*NonLin*/false);
	bks[0].SetNy (ndiv);

	// Generate
	Mesh::Structured mesh(/*Is3D*/false);
	if (is_o2) mesh.SetO2();
	mesh.SetBlocks (bks);
	mesh.Generate  (true);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Data and solver
	FEM::Data   dat (2); // 2D
	FEM::Solver sol (dat,"tpstrain05");

	// Elements attributes
	String prms; prms.Printf("E=%f nu=%f",E,nu);
	String geom; geom = (is_o2 ? "Quad8" : "Quad4");
	FEM::EAtts_T eatts(1);
	eatts = T(-1, geom.CStr(), "PStrain", "LinElastic", prms.CStr(), "ZERO", "gam=20", FNULL, true);

	// Set geometry: nodes and elements
	dat.SetNodesElems (&mesh, &eatts);

	// Stage # 1 -----------------------------------------------------------
	FEM::EBrys_T ebrys;
	ebrys.Push  (T(-30, "uy", 0.0));
	ebrys.Push  (T(-40, "ux", 0.0));
	ebrys.Push  (T(-20, "Q",   p0));
	ebrys.Push  (T(-10, "Q",   p1));
	dat.SetBrys (&mesh, NULL, &ebrys, NULL);
	sol.SolveWithInfo(/*NDiv*/1, /*DTime*/0.0);

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Stress
	Array <double> err_sR;
	Array <double> err_sT;
	Array <double> err_uR;
	for (size_t i=0; i<dat.NElems(); ++i)
	{
		for (size_t j=0; j<dat.Ele(i)->NNodes(); ++j)
		{
			// Analytical
			double x = dat.Ele(i)->Nod(j)->X();
			double y = dat.Ele(i)->Nod(j)->Y();
			double L = sqrt(x*x+y*y);
			double sigRc, sigTc, uRc; // correct stress components
			Analyitical(r, R, p0, p1, E, nu, L, sigRc, sigTc, uRc);

			// Numerical
			double c     = x/L;
			double s     = y/L;
			double cc    = c*c;
			double ss    = s*s;
			double sc    = s*c;
			double Sx    = dat.Ele(i)->Val(j,"Sx");
			double Sy    = dat.Ele(i)->Val(j,"Sy");
			double Sxy   = dat.Ele(i)->Val(j,"Sxy");
			double ux    = dat.Ele(i)->Val(j,"ux");
			double uy    = dat.Ele(i)->Val(j,"uy");
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

	// Return error flag
	if (max_err_sR>tol_sR || max_err_sT>tol_sT || max_err_uR>tol_uR) return 1;
	else return 0;

}
MECHSYS_CATCH
