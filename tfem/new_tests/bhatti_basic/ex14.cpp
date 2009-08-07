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

/*  Bhatti (2005): Example 1.4, p25
 *  ===============================
 */

// STL
#include <iostream>

// boost::python
#include <boost/python.hpp>
namespace BPy = boost::python;

#define USE_BOOST_PYTHON

// MechSys
#include "fem/data.h"
#include "fem/solver.h"
#include "fem/elems/rod.h"
#include "models/equilibs/rodelastic.h"
#include "util/fatal.h"
#include "util/table.h"
#include "mesh/mesh.h"

using std::cout;
using std::endl;

#define T boost::make_tuple

int main(int argc, char **argv) try
{
	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
	
	Mesh::Generic mesh(/*Is3D*/false);
	mesh.SetNVerts  (4);
	mesh.SetNElems  (5);
	mesh.SetVert    (0, true,  0.0,  0.0, 0, -100); // true => OnBry
	mesh.SetVert    (1, true, 1500, 3500, 0, -200);
	mesh.SetVert    (2, true,  0.0, 5000, 0,    0);
	mesh.SetVert    (3, true, 5000, 5000, 0, -100);
	mesh.SetElem    (0, -1, true, VTK_LINE); // true => OnBry
	mesh.SetElem    (1, -1, true, VTK_LINE);
	mesh.SetElem    (2, -2, true, VTK_LINE);
	mesh.SetElem    (3, -2, true, VTK_LINE);
	mesh.SetElem    (4, -3, true, VTK_LINE);
	mesh.SetElemCon (0, 0, 0);  mesh.SetElemCon(0, 1, 1);
	mesh.SetElemCon (1, 0, 1);  mesh.SetElemCon(1, 1, 3);
	mesh.SetElemCon (2, 0, 0);  mesh.SetElemCon(2, 1, 2);
	mesh.SetElemCon (3, 0, 2);  mesh.SetElemCon(3, 1, 3);
	mesh.SetElemCon (4, 0, 1);  mesh.SetElemCon(4, 1, 2);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////
	
	// Data and Solver
	FEM::Data   dat (2); // 2D
	FEM::Solver sol (dat, "ex14");

	// Parameters and initial value
	FEM::EAtts_T eatts(3);
	String p1; p1.Printf("E=%f A=%f", 200000.0, 4000.0);
	String p2; p2.Printf("E=%f A=%f", 200000.0, 3000.0);
	String p3; p3.Printf("E=%f A=%f",  70000.0, 2000.0);
	eatts = T(-1, "Lin2", "Rod", "RodElastic", p1.CStr(), "ZERO", "gam=0", FNULL, true),
	        T(-2, "Lin2", "Rod", "RodElastic", p2.CStr(), "ZERO", "gam=0", FNULL, true),
	        T(-3, "Lin2", "Rod", "RodElastic", p3.CStr(), "ZERO", "gam=0", FNULL, true);

	// Set geometry: nodes and elements
	dat.SetOnlyFrame  (true);
	dat.SetNodesElems (&mesh, &eatts);

	// Stage # 1 -----------------------------------------------------------
	FEM::NBrysT_T nbrys(3);
	nbrys = T(-100, "ux", 0.0),
	        T(-100, "uy", 0.0),
	        T(-200, "fy", -150000);
	dat.SetNBrys      (mesh, nbrys);
	sol.SolveWithInfo (/*NDiv*/1, /*DTime*/0.0);

	//////////////////////////////////////////////////////////////////////////////////////// Output ////

	dat.PrintResults();

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	/*
	Table nod_sol("                  ux                      uy                     Rux                     Ruy", 4,
	              0.000000000000000e+00,  0.000000000000000e+00,  5.492667465378455e+04,  1.599266746537845e+05,
	              5.389536380057676e-01, -9.530613006371175e-01,  0.000000000000000e+00,  0.000000000000000e+00,
	              2.647036149579491e-01, -2.647036149579490e-01,  0.000000000000000e+00,  0.000000000000000e+00,
	              0.000000000000000e+00,  0.000000000000000e+00, -5.492667465378455e+04, -9.926674653784567e+03);
	*/

	// correct solution
	Table nod_sol("                  ux                      uy", /*NRows*/4,
	              0.000000000000000e+00,  0.000000000000000e+00,
	              5.389536380057676e-01, -9.530613006371175e-01,
	              2.647036149579491e-01, -2.647036149579490e-01,
	              0.000000000000000e+00,  0.000000000000000e+00);

	Table ele_sol("                   Ea                      Sa                       N", /*NRows*/5,
	              -1.742954548428455e-04, -3.485909096856910e+01, -1.394363638742764e+05,
	              -3.149970910789726e-05, -6.299941821579453e+00, -2.519976728631781e+04,
	              -5.294072299158981e-05, -1.058814459831796e+01, -3.176443379495388e+04,
	              -5.294072299158982e-05, -1.058814459831796e+01, -3.176443379495389e+04,
	               3.208692362423290e-04,  2.246084653696303e+01,  4.492169307392606e+04);

	// error tolerance
	FEM::StrDbl_t nod_tol, ele_tol;
	nod_tol["ux"] = 1.0e-15;
	nod_tol["uy"] = 1.0e-15;
	ele_tol["Ea"] = 1.0e-15;
	ele_tol["Sa"] = 1.0e-14;
	ele_tol["N" ] = 1.0e-10;

	// Return error flag
	return dat.CheckError(nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH
