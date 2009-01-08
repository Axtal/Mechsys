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
#include "fem/elems/hex8.h"
#include "fem/equilibelem.h"
#include "models/equilibs/linelastic.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using Util::_4;
using Util::_8s;

#define T boost::make_tuple

int main(int argc, char **argv) try
{
	// Description:
	// Test to evaluate the rotation of a loaded cube with just vertical constraints

	// constants
	double E     = 200.0; // Young
	double nu    = 0.25;  // Poisson
	double q     = 1.0;   // Downward vertical pressure
	int    ndiv  = 2;     // number of divisions along x, y, and z
	bool   is_o2 = false; // use high order elements?
	String linsol("UM");  // UMFPACK

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
	cout << "Input: " << argv[0] << "  is_o2  ndiv  linsol(LA,UM,SLU)\n";
	if (argc>=2) is_o2      = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) ndiv       =  atof(argv[2]);
	if (argc>=4) linsol.Printf("%s",argv[3]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
	
    /*
                      4----------------7  
                    ,'|              ,'| 
                  ,'  |            ,'  | 
                ,'    | -6    -1 ,'    | 
              ,'      |        ,'      | 
            5'===============6'        | 
            |         |      |    -4   | 
            |    -3   |      |         | 
            |         0- - - | -  - - -3  
            |       ,'       |       ,'  
            |     ,' -2      |     ,'    
            |   ,'        -5 |   ,'      
            | ,'             | ,'        
            1----------------2'          
    */

	// Blocks
	Array<Mesh::Block> bks(1);

	// Block # 0 --------------------------------
    Mesh::Verts_T ve0( 8);
    Mesh::Edges_T ed0(12);
    Mesh::FTags_T ft0( 6);
    ve0 = T(0, 0., 0., 0.), T(1, 1., 0., 0.), T(2, 1., 1., 0.), T(3, 0., 1., 0.),
          T(4, 0., 0., 1.), T(5, 1., 0., 1.), T(6, 1., 1., 1.), T(7, 0., 1., 1.);
    ed0 = T(0,1), T(1,2), T(2,3), T(3,0),
          T(4,5), T(5,6), T(6,7), T(7,4),
          T(0,4), T(1,5), T(2,6), T(3,7);
    ft0 = T(0,3,7,4,-1), T(1,2,6,5,-2), T(1,0,4,5,-3), T(2,3,7,6,-4), T(0,1,2,3,-5), T(4,5,6,7,-6);
    bks[0].Set   (-1, ve0, ed0, NULL, &ft0, /*orig*/0, /*xplus*/1, /*yplus*/3, /*zplus*/4);
	bks[0].SetNx (ndiv);
	bks[0].SetNy (ndiv);
	bks[0].SetNz (ndiv);

	// Generate
	Mesh::Structured mesh(/*Is3D*/true);
	if (is_o2) mesh.SetO2();
	mesh.SetBlocks (bks);
	mesh.Generate  (true);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Data and solver
	FEM::Data   dat (3); // 3D
	FEM::Solver sol (dat,"tsolid02");

	// Element attributes
	String prms; prms.Printf("E=%f nu=%f",E,nu);
	String geom; geom = (is_o2 ? "Hex20" : "Hex8");
	FEM::EAtts_T eatts(1);
	eatts = T(-1, geom.CStr(), "Equilib", "LinElastic", prms.CStr(), "ZERO", "gam=20", true);

	// Set geometry: nodes and elements
	dat.SetNodesElems (&mesh, &eatts);

	// Stage # 1 --------------------------------------------
	FEM::NBrys_T nbrys;
	FEM::FBrys_T fbrys;
	nbrys.Push        (T(0.5, 0.5, 0., "ux", 0.0));
	nbrys.Push        (T(0.5, 0.5, 0., "uy", 0.0));
	fbrys.Push        (T(-5, "uz", 0.0));
	fbrys.Push        (T(-6, "fz",  -q));
	dat.SetBrys       (&mesh, NULL, NULL, &fbrys);
	sol.SolveWithInfo (/*NDiv*/1, /*DTime*/0.0);

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
	for (size_t i=0; i<dat.NElems(); ++i)
	{
		for (size_t j=0; j<dat.Ele(i)->NNodes(); ++j)
		{
			// eps
			err_eps.Push( fabs(dat.Ele(i)->Val(j,"Ex" ) - Ex ) / (1.0+fabs(Ex )) );
			err_eps.Push( fabs(dat.Ele(i)->Val(j,"Ey" ) - Ey ) / (1.0+fabs(Ey )) );
			err_eps.Push( fabs(dat.Ele(i)->Val(j,"Ez" ) - Ez ) / (1.0+fabs(Ez )) );
			err_eps.Push( fabs(dat.Ele(i)->Val(j,"Exy") - Exy) / (1.0+fabs(Exy)) );
			err_eps.Push( fabs(dat.Ele(i)->Val(j,"Eyz") - Eyz) / (1.0+fabs(Eyz)) );
			err_eps.Push( fabs(dat.Ele(i)->Val(j,"Ezx") - Ezx) / (1.0+fabs(Ezx)) );
			// sig
			err_sig.Push( fabs(dat.Ele(i)->Val(j,"Sx" ) - Sx ) / (1.0+fabs(Sx )) );
			err_sig.Push( fabs(dat.Ele(i)->Val(j,"Sy" ) - Sy ) / (1.0+fabs(Sy )) );
			err_sig.Push( fabs(dat.Ele(i)->Val(j,"Sz" ) - Sz ) / (1.0+fabs(Sz )) );
			err_sig.Push( fabs(dat.Ele(i)->Val(j,"Sxy") - Sxy) / (1.0+fabs(Sxy)) );
			err_sig.Push( fabs(dat.Ele(i)->Val(j,"Syz") - Syz) / (1.0+fabs(Syz)) );
			err_sig.Push( fabs(dat.Ele(i)->Val(j,"Szx") - Szx) / (1.0+fabs(Szx)) );
		}
	}

	// Displacements
	for (size_t i=0; i<dat.NNodes(); ++i)
	{
		double ux_correct = -Ex*dat.Nod(i)->X();
		double uy_correct = -Ey*dat.Nod(i)->Y();
		double uz_correct = -Ez*dat.Nod(i)->Z();
		err_dis.Push ( fabs(dat.Nod(i)->Val("ux") - ux_correct) / (1.0+fabs(ux_correct)) );
		err_dis.Push ( fabs(dat.Nod(i)->Val("uy") - uy_correct) / (1.0+fabs(uy_correct)) );
		err_dis.Push ( fabs(dat.Nod(i)->Val("uz") - uz_correct) / (1.0+fabs(uz_correct)) );
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
	cout << endl;
	cout << _4<< ""    << _8s<<"Min"       << _8s<<"Mean"                                                        << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	cout << _4<< "Eps" << _8s<<min_err_eps << _8s<<err_eps.Mean() << (max_err_eps>tol_eps?"[1;31m":"[1;32m") << _8s<<max_err_eps << "[0m" << _8s<<err_eps.Norm() << endl;
	cout << _4<< "Sig" << _8s<<min_err_sig << _8s<<err_sig.Mean() << (max_err_sig>tol_sig?"[1;31m":"[1;32m") << _8s<<max_err_sig << "[0m" << _8s<<err_sig.Norm() << endl;
	cout << _4<< "Dis" << _8s<<min_err_dis << _8s<<err_dis.Mean() << (max_err_dis>tol_dis?"[1;31m":"[1;32m") << _8s<<max_err_dis << "[0m" << _8s<<err_dis.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err_eps>tol_eps || max_err_sig>tol_sig || max_err_dis>tol_dis) return 1;
	else return 0;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
