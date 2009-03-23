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

/*                      
      ///@-----------@ |
      ///|           | | F
      ///|           | |
      ///@-----------@ V
*/

// STL
#include <iostream>

// MechSys
#include "fem/data.h"
#include "fem/elems/hex8.h"
#include "fem/elems/hex20.h"
#include "fem/equilibelem.h"
#include "models/equilibs/linelastic.h"
#include "fem/solver.h"
#include "fem/output.h"
#include "util/exception.h"
#include "linalg/matrix.h"
#include "mesh/unstructured.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_8s;

#define T boost::make_tuple

int main(int argc, char **argv) try
{
	// Constants
	double E       = 10.7e+6; // Young (psi)
	double nu      = 0.3;     // Poisson
	double F       = 120.0;   // Load
	size_t nx      = 1;       // ndivs along x
	size_t ny      = 20;    // ndivs along y
	size_t nz      = 10;    // ndivs along z
	bool   is_o2   = false;   // use high order elements?

	// Input
	cout << "Input: " << argv[0] << "  is_o2  maxarea\n";
	if (argc>=2) is_o2 = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) nx    =  atof(argv[2]);
	if (argc>=4) ny    =  atof(argv[3]);
	if (argc>=5) nz    =  atof(argv[4]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Dimensions
	double L = 3.0;
	double H = 0.6;
	double B = 0.1;

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
           H|       ,'       |       ,'  
            |     ,' -2      |     ,'    
            |   ,'        -5 |   ,' B    
            | ,'             | ,'        
            1----------------2'          
	                L
    */

	// Blocks
	Array<Mesh::Block> bks(1);

	// Block # 0 --------------------------------
    Mesh::Verts_T ve0( 8);
    Mesh::Edges_T ed0(12);
    Mesh::FTags_T ft0( 6);
    ve0 = T(0, 0., 0., 0.), T(1, B, 0., 0.), T(2, B, L, 0.), T(3, 0., L, 0.),
          T(4, 0., 0., H),  T(5, B, 0., H),  T(6, B, L, H),  T(7, 0., L, H),
    ed0 = T(0,1), T(1,2), T(2,3), T(3,0),
          T(4,5), T(5,6), T(6,7), T(7,4),
          T(0,4), T(1,5), T(2,6), T(3,7);
    ft0 = T(0,3,7,4,-1), T(1,2,6,5,-2), T(1,0,4,5,-3), T(2,3,7,6,-4), T(0,1,2,3,-5), T(4,5,6,7,-6);
    bks[0].Set   (-1, ve0, ed0, NULL, &ft0, /*orig*/0, /*xplus*/1, /*yplus*/3, /*zplus*/4);
	bks[0].SetNx (nx);
	bks[0].SetNy (ny);
	bks[0].SetNz (nz);

	// Generate
	Mesh::Structured mesh(/*Is3D*/true);
	if (is_o2) mesh.SetO2();
	mesh.SetBlocks (bks);
	mesh.Generate  (true);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Data and Solver
	FEM::Data   dat (3);
	FEM::Solver sol (dat, "tstatic122");

	// Elements attributes
	String prms; prms.Printf("E=%f nu=%f",E,nu);
	String geom; geom = (is_o2 ? "Hex20" : "Hex8");
	FEM::EAtts_T eatts(1);
	eatts = T(-1, geom.CStr(), "Equilib", "LinElastic", prms.CStr(), "ZERO", "gam=20", FNULL, true);

	// Set geometry: nodes and elements
	dat.SetNodesElems (&mesh, &eatts);

	// State # 1 -----------------------------------------------
	FEM::FBrys_T fbrys;
	fbrys.Push  (T(-3, "ux", 0.0));
	fbrys.Push  (T(-3, "uy", 0.0));
	fbrys.Push  (T(-3, "uz", 0.0));
	fbrys.Push  (T(-4, "fz", -F/0.06));
	dat.SetBrys (&mesh, NULL, NULL, &fbrys);
	sol.SolveWithInfo();

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Check
    Array<double> err_dis_z;
	double I = B*H*H*H/12.0; // Section momentum of inertia

	for (size_t i=0; i<dat.NNodes(); i++)
	{
		double y   = dat.Nod(i)->Y();
		double uz  = dat.Nod(i)->Val("uz");
		double auz = F/(E*I)*(y*y*y/6.0 - y*y*L/2.0); // Analitic solution for the vertical displacement
		err_dis_z.Push(fabs(uz - auz));
		cout << "iNode: " << i << "  Uz: " << uz << "    Analytic Uz: " <<  auz << endl;
	}

	// Error summary
	double tol_dis     = 1.0e-3;
	double min_err_dis = err_dis_z[err_dis_z.Min()];
	double max_err_dis = err_dis_z[err_dis_z.Max()];
	cout << _4<< ""    << _8s<<"Min"       << _8s<<"Mean"                                                          << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	cout << _4<< "Dis" << _8s<<min_err_dis << _8s<<err_dis_z.Mean() << (max_err_dis>tol_dis?"[1;31m":"[1;32m") << _8s<<max_err_dis << "[0m" << _8s<<err_dis_z.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err_dis>tol_dis) return 1;
	else return 0;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
