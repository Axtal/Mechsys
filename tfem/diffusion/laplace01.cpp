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
      u=100*sin(PI*x/10)
       //@-----------@       ---
       //|           |        |
       //|           |        |
  u=0  //|           |f=0    H=10
       //|           |        |
       //|    L=5    |        |
       //@-----------@       ---
       ///////////////
              u=0     
*/

// STL
#include <iostream>
#include <fstream>
#include <cmath>

// MechSys
#include "fem/data.h"
#include "fem/solver.h"
#include "fem/elems/tri3.h"
#include "fem/elems/tri6.h"
#include "fem/elems/quad4.h"
#include "fem/elems/quad8.h"
#include "fem/diffusionelem.h"
#include "models/diffusions/lindiffusion.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "mesh/structured.h"
#include "mesh/unstructured.h"

using std::cout;
using std::endl;
using LinAlg::Vector;
using Util::_4;
using Util::_6;
using Util::_8s;
using Util::PI;

#define T boost::make_tuple

int main(int argc, char **argv) try
{
	// Constants
	bool   is_o2   = false; // use high order elements?
	double maxarea = 1.0;   // maximum area for each triangle
	bool   is_quad = false; // use quadrilateral mesh?
	size_t ndiv    = 4;     // number of divisions
	String linsol("UM");    // UMFPACK

	// Input
	cout << "Input: " << argv[0] << "  is_o2  maxarea(1.0)  is_quad  ndiv  linsol(LA,UM,SLU)\n";
	if (argc>=2) is_o2      = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) maxarea    =  atof(argv[2]);
	if (argc>=4) is_quad    = (atoi(argv[3])>0 ? true : false);
	if (argc>=5) ndiv       =  atof(argv[4]);
	if (argc>=6) linsol.Printf("%s",argv[5]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Mesh
	Mesh::Generic     * msh;
	Mesh::Structured    mss(/*Is3D*/false);
	Mesh::Unstructured  msu(/*Is3D*/false);
	if (is_quad) // quadrangles
	{
		// Blocks
		Array<Mesh::Block> bks(1);

		// Block # 0 --------------------------------
		Mesh::Verts_T ve0(4);
		Mesh::Edges_T ed0(4);
		Mesh::ETags_T et0(3);
		ve0 = T(0,0.0,0.0,0.), T(1,5.0,0.0,0.), T(2,5.0,10.0,0.), T(3,0.0,10.0,0.);
		ed0 = T(0,1), T(1,2), T(2,3), T(3,0);
		et0 = T(3,0,-10), T(1,2,-20), T(0,1,-10);
		bks[0].Set (-1, ve0, ed0, &et0, NULL, /*orig*/0, /*xplus*/1, /*yplus*/3);
		bks[0].SetNx (ndiv);
		bks[0].SetNy (ndiv);

		// Generate
		if (is_o2) mss.SetO2();
		mss.SetBlocks (bks);
		mss.Generate  (true);
		msh = &mss;
	}
	else // triangles
	{
		// Polygon
		msu.SetPolySize    (/*NPoints*/4, /*NSegments*/4, /*NRegions*/1);
		msu.SetPolyPoint   (0, /*X*/ 0.0, /*Y*/ 0.0);
		msu.SetPolyPoint   (1, /*X*/ 5.0, /*Y*/ 0.0);
		msu.SetPolyPoint   (2, /*X*/ 5.0, /*Y*/10.0);
		msu.SetPolyPoint   (3, /*X*/ 0.0, /*Y*/10.0);
		msu.SetPolySegment (0, /*iPointLeft*/0, /*iPointRight*/1, /*Tag*/-10);
		msu.SetPolySegment (1, /*iPointLeft*/1, /*iPointRight*/2, /*Tag*/-20);
		msu.SetPolySegment (2, /*iPointLeft*/2, /*iPointRight*/3);
		msu.SetPolySegment (3, /*iPointLeft*/3, /*iPointRight*/0, /*Tag*/-10);
		msu.SetPolyRegion  (0, /*Tag*/-1, maxarea, /*X*/2.5, /*Y*/5.0);

		// Generate
		if (is_o2) msu.SetO2();
		msu.Generate (true);
		msh = &msu;
	}

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Data dat(2); // 2D

	// Edges brys (the order matters!)
	FEM::EBrys_T ebrys(2);
	ebrys = T(-20, "f", 0.0),
	        T(-10, "u", 0.0);

	// Elements attributes
	FEM::EAtts_T eatts(1);
	String geom;
	if (is_quad) geom = (is_o2 ? "Quad8" : "Quad4");
	else         geom = (is_o2 ? "Tri6"  : "Tri3" );
	eatts = T(-1, geom.CStr(), "Diffusion", "LinDiffusion", "k=1.0", "", "s=0.0", FNULL, true);

	// Set geometry: nodes, elements, attributes, and boundaries
	dat.SetNodesElems (msh, &eatts);
	dat.SetBrys       (msh, NULL, &ebrys, NULL);

	// Set prescribed u = 100*sin(PI*x/10) on top nodes
	for (size_t i=0; i<dat.NNodes(); ++i)
	{
		double x = dat.Nod(i)->X();
		double y = dat.Nod(i)->Y();
		if (y>10.0-1e-7) // top node
			dat.Nod(i)->Bry("u", 100.0*sin(PI*x/10.0));
	}

	// Solve
	FEM::Solver sol(dat, "tlaplace01");
	sol.SolveWithInfo();

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Check
	Array<double> err_u;
	//cout << _6<<"Node #" << _8s<<"u" << _8s<<"fcorrect" << _8s<<"f" << endl;
	for (size_t i=0; i<dat.NNodes(); ++i)	
	{
		double x     = dat.Nod(i)->X();
		double y     = dat.Nod(i)->Y();
		double u     = dat.Nod(i)->Val("u");
		double ucorr = 100.0*sinh(PI*y/10.0)*sin(PI*x/10.0)/sinh(PI);
		err_u.Push ( fabs(u-ucorr) / (1.0+fabs(ucorr)) );
	}

	// Error summary
	double tol_u     = 3.8e-2;
	double min_err_u = err_u[err_u.Min()];
	double max_err_u = err_u[err_u.Max()];
	cout << _4<< ""  << _8s<<"Min"     << _8s<<"Mean"                                                  << _8s<<"Max"                << _8s<<"Norm"       << endl;
	cout << _4<< "u" << _8s<<min_err_u << _8s<<err_u.Mean() << (max_err_u>tol_u?"[1;31m":"[1;32m") << _8s<<max_err_u << "[0m" << _8s<<err_u.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err_u>tol_u) return 1;
	else return 0;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
