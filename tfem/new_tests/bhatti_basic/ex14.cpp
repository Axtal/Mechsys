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
using LinAlg::Matrix;

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

	/*
	// Output: Nodes
	cout << _6<<"Node #" << _8s<<"ux" << _8s<<"uy" << _8s<<"fx"<< _8s<<"fy" << endl;
	for (size_t i=0; i<dat.NNodes(); ++i)
		cout << _6<<i << _8s<<dat.Nod(i)->Val("ux") <<  _8s<<dat.Nod(i)->Val("uy") << _8s<<dat.Nod(i)->Val("fx") << _8s<<dat.Nod(i)->Val("fy") << endl;
	cout << endl;

	// Output: Elements
	cout << _6<<"Elem #" << _8s<<"Sa(left)" << _8s<<"Sa(right)" << _8s<<"Ea(left)" << _8s<<"Ea(right)" << endl;
	for (size_t i=0; i<dat.NElems(); ++i)
	{
		dat.Ele(i)->CalcDeps();
		cout << _6<<i;
		for (size_t j=0; j<dat.Ele(i)->NNodes(); ++j) cout << _8s<<dat.Ele(i)->Val(j, "Sa");
		for (size_t j=0; j<dat.Ele(i)->NNodes(); ++j) cout << _8s<<dat.Ele(i)->Val(j, "Ea");
		cout << endl;
	}
	cout << endl;
	*/

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	/*
	Py_Initialize();
	BPy::object main_module((BPy::handle<>(BPy::borrowed(PyImport_AddModule("__main__")))));
	BPy::object main_namespace = main_module.attr("__dict__");
	BPy::handle<> ignored((PyRun_String("print \"Hello, World\"", Py_file_input, main_namespace.ptr(), main_namespace.ptr())));
	*/

	Table nod_sol("                  ux                      uy                     Rux                     Ruy", /*NRows*/4,
	              0.000000000000000e+00,  0.000000000000000e+00,  5.492667465378455e+04,  1.599266746537845e+05,
	              5.389536380057676e-01, -9.530613006371175e-01,  0.000000000000000e+00,  0.000000000000000e+00,
	              2.647036149579491e-01, -2.647036149579490e-01,  0.000000000000000e+00,  0.000000000000000e+00,
	              0.000000000000000e+00,  0.000000000000000e+00, -5.492667465378455e+04, -9.926674653784567e+03);

	Table ele_sol("                   ea                      sa                      Fa", /*NRows*/5,
	              -1.742954548428455e-04, -3.485909096856910e+01, -1.394363638742764e+05,
	              -3.149970910789726e-05, -6.299941821579453e+00, -2.519976728631781e+04,
	              -5.294072299158981e-05, -1.058814459831796e+01, -3.176443379495388e+04,
	              -5.294072299158982e-05, -1.058814459831796e+01, -3.176443379495389e+04,
	               3.208692362423290e-04,  2.246084653696303e+01,  4.492169307392606e+04);

	/*
	Table nod_sol;
	nod_sol["ux"].Resize(dat.NNodes());
	nod_sol["ux"] = 0.000000000000000e+00, 5.389536380057676e-01, 2.647036149579491e-01, 0.000000000000000e+00;
	cout << nod_sol["ux"].Size() << endl;
	*/

	/*
	// Displacements
	Array<double> err_u(6);
	err_u[0] =  fabs(dat.Nod(0)->Val("ux") - ( 0.0));
	err_u[1] =  fabs(dat.Nod(0)->Val("uy") - (-0.5));
	err_u[2] =  fabs(dat.Nod(1)->Val("ux") - ( 0.0));
	err_u[3] =  fabs(dat.Nod(1)->Val("uy") - ( 0.4));
	err_u[4] =  fabs(dat.Nod(2)->Val("ux") - (-0.5));
	err_u[5] =  fabs(dat.Nod(2)->Val("uy") - ( 0.2));

	// Forces
	Array<double> err_f(6);
	err_f[0] = fabs(dat.Nod(0)->Val("fx") - (-2.0));
	err_f[1] = fabs(dat.Nod(0)->Val("fy") - (-2.0));
	err_f[2] = fabs(dat.Nod(1)->Val("fx") - ( 0.0));
	err_f[3] = fabs(dat.Nod(1)->Val("fy") - ( 1.0));
	err_f[4] = fabs(dat.Nod(2)->Val("fx") - ( 2.0));
	err_f[5] = fabs(dat.Nod(2)->Val("fy") - ( 1.0));

	// Correct axial normal stresses
	Array<double> err_s(dat.NElems());
	err_s.SetValues(0.0);
	err_s[1] = fabs(dat.Ele(1)->Val(0, "N") - (-1.0));

	// Error summary
	double tol_u     = 1.0e-7;
	double tol_f     = 1.0e-7;
	double tol_s     = 1.0e-7;
	double min_err_u = err_u[err_u.Min()];
	double max_err_u = err_u[err_u.Max()];
	double min_err_f = err_f[err_f.Min()];
	double max_err_f = err_f[err_f.Max()];
	double min_err_s = err_s[err_s.Min()];
	double max_err_s = err_s[err_s.Max()];
	cout << _4<< ""  << _8s<<"Min"     << _8s<<"Mean"                                                  << _8s<<"Max"                << _8s<<"Norm"       << endl;
	cout << _4<< "u" << _8s<<min_err_u << _8s<<err_u.Mean() << (max_err_u>tol_u?"[1;31m":"[1;32m") << _8s<<max_err_u << "[0m" << _8s<<err_u.Norm() << endl;
	cout << _4<< "f" << _8s<<min_err_f << _8s<<err_f.Mean() << (max_err_f>tol_f?"[1;31m":"[1;32m") << _8s<<max_err_f << "[0m" << _8s<<err_f.Norm() << endl;
	cout << _4<< "N" << _8s<<min_err_s << _8s<<err_s.Mean() << (max_err_s>tol_s?"[1;31m":"[1;32m") << _8s<<max_err_s << "[0m" << _8s<<err_s.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err_u>tol_u || max_err_f>tol_f || max_err_s>tol_s) return 1;
	else return 0;
	*/
	return 1;
}
MECHSYS_CATCH
