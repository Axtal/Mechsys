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

// STL
#include <iostream>
#include <fstream>

// MechSys
#include "fem/data.h"
#include "fem/solver.h"
#include "fem/elems/hex8.h"
#include "fem/equilibelem.h"
#include "models/equilibs/camclay.h"
#include "fem/output.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "linalg/matrix.h"
#include "tensors/tensors.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Tensors::Tensor1;
using Util::_8s;

#define T boost::make_tuple

int main(int argc, char **argv) try
{
	// constants
	double lam   = 0.0891;  // lambda
	double kap   = 0.0196;  // kappa
	double phics = 31.6;    // shear angle at CS
	double nu    = 0.2;     // Poisson ratio
	double p_ini = 198.0;   // initial mean stress (kPa)
	double v_ini = 1.6910;  // initial specific volume
	bool   is_o2 = false;   // use high order elements?
	int    ndiv  = 1;       // number of divisions along x, y and z
	String linsol("UM");    // UMFPACK

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
	bks[0].SetNx (3);
	bks[0].SetNy (3);
	bks[0].SetNz (3);

	// Generate
	Mesh::Structured mesh(/*Is3D*/false);
	mesh.SetBlocks (bks);              // Set Blocks
	if (is_o2) mesh.SetO2();           // Non-linear elements
	mesh.Generate (/*WithInfo*/ true); // Discretize domain

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Data and Solver
	FEM::Data   dat (3); // 3D
	FEM::Solver sol (dat, "tcamclay01");
	sol.SetLinSol   (linsol.CStr());

	// Element attributes
	String prms; prms.Printf("lam=%f kap=%f nu=%f phics=%f",lam,kap,nu,phics);
	String inis; inis.Printf("Sx=%f Sy=%f Sz=%f Sxy=0 Syz=0 Szx=0 v=%f",p_ini,p_ini,p_ini,v_ini);
	String geom; geom = (is_o2 ? "Hex20" : "Hex8");
	FEM::EAtts_T eatts(1);
	eatts = T(-1, geom.CStr(), "Equilib", "CamClay", prms.CStr(), inis.CStr(), "gam=20", FNULL, true);

	// Set geometry: nodes and elements
	dat.SetNodesElems (&mesh, &eatts);

	// Open file
	std::ofstream res("tcamclay01.cal", std::ios::out);
	res << _8s<<"Sx" << _8s<<"Sy" << _8s<<"Sz" << _8s<<"Ex" << _8s<<"Ey" << _8s<<"Ez" << _8s<<"p" << _8s<<"q" << _8s<<"Ev" << _8s<<"Ed" << "\n";

    // Output
    res << _8s<<dat.Ele(13)->Val("Sx") << _8s<<dat.Ele(13)->Val("Sy") << _8s<<dat.Ele(13)->Val("Sz");
    res << _8s<<dat.Ele(13)->Val("Ex") << _8s<<dat.Ele(13)->Val("Ey") << _8s<<dat.Ele(13)->Val("Ez");
    res << _8s<<dat.Ele(13)->Val("p")  << _8s<<dat.Ele(13)->Val("q")  << _8s<<dat.Ele(13)->Val("Ev") << _8s<<dat.Ele(13)->Val("Ed") << "\n";

	// Solve each stage -- stress path
	Array<Tensor1> dtrac;  dtrac.Resize(2); // stress path (delta traction)
	dtrac[0] = 0.0, 0.0, -2.0;              // dtracx, dtracy, dtracz
	dtrac[1] = 0.0, 0.0,  3.0;
	for (size_t i=0; i<dtrac.Size(); ++i)
	{
        // Faces brys
        FEM::FBrys_T fbrys;
        fbrys.Push (T(-1, "ux", 0.0));
        fbrys.Push (T(-2, "fx", dtrac[i](0)));
        fbrys.Push (T(-3, "uy", 0.0));
        fbrys.Push (T(-4, "fy", dtrac[i](1)));
        fbrys.Push (T(-5, "uz", 0.0));
        fbrys.Push (T(-6, "fz", dtrac[i](2)));

        // Set boundary conditions
        dat.SetBrys (&mesh, NULL, NULL, &fbrys);

        // Solve
        sol.SolveWithInfo();

        // Output
        res << _8s<<dat.Ele(13)->Val("Sx") << _8s<<dat.Ele(13)->Val("Sy") << _8s<<dat.Ele(13)->Val("Sz");
        res << _8s<<dat.Ele(13)->Val("Ex") << _8s<<dat.Ele(13)->Val("Ey") << _8s<<dat.Ele(13)->Val("Ez");
        res << _8s<<dat.Ele(13)->Val("p")  << _8s<<dat.Ele(13)->Val("q")  << _8s<<dat.Ele(13)->Val("Ev") << _8s<<dat.Ele(13)->Val("Ed") << "\n";
	}

	// Close file
	res.close();
	cout << "[1;34mFile <tcamclay01.cal> saved.[0m" << endl;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
