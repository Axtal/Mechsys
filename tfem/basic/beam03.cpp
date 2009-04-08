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

/* Yuan-Yu Hsieh & S. T. Mau, Elementary Theory of Structured
 * Example 3-5, Fig 3-12, p46 */

// STL
#include <iostream>
#include <fstream>

// MechSys
#include "fem/data.h"
#include "fem/solver.h"
#include "fem/elems/beam.h"
#include "models/equilibs/beamelastic.h"
#include "util/exception.h"
#include "mesh/mesh.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_6;
using Util::_8s;
using boost::make_tuple;

int main(int argc, char **argv) try
{
	// Data
	double E   = 1.1;
	double A   = 2.2;
	double Izz = 1.8;
	double Qb  = -1.3;
	double P   = -4.0;
	double l   = 2.1;

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
	
	Mesh::Generic mesh(/*Is3D*/false);
	mesh.SetNVerts  (5);
	mesh.SetNElems  (4);
	mesh.SetVert    (0, true,    0.0, 0.0); // true => OnBry
	mesh.SetVert    (1, true, 0.99*l, 0.0);
	mesh.SetVert    (2, true,      l, 0.0);
	mesh.SetVert    (3, true,  1.5*l, 0.0);
	mesh.SetVert    (4, true,  2.0*l, 0.0);
	mesh.SetElem    (0, -1, true, VTK_LINE);
	mesh.SetElem    (1, -2, true, VTK_LINE);
	mesh.SetElem    (2, -3, true, VTK_LINE);
	mesh.SetElem    (3, -3, true, VTK_LINE);
	mesh.SetElemCon (0, 0, 0);  mesh.SetElemCon(0, 1, 1);
	mesh.SetElemCon (1, 0, 1);  mesh.SetElemCon(1, 1, 2);
	mesh.SetElemCon (2, 0, 2);  mesh.SetElemCon(2, 1, 3);
	mesh.SetElemCon (3, 0, 3);  mesh.SetElemCon(3, 1, 4);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////
	
	// Data and solver
	FEM::Data   dat (2); // 2D
	FEM::Solver sol (dat, "beam03");

	// Elements attributes
	FEM::EAtts_T eatts;
	String prms; prms.Printf("E=%f A=%f Izz=%f",E,A,Izz);
	eatts.Push (make_tuple(-1, "Lin2", "Beam", "BeamElastic", prms.CStr(), "ZERO", "cq=1",        FNULL, true));
	eatts.Push (make_tuple(-2, "Lin2", "Beam", "BeamElastic", prms.CStr(), "ZERO", "cq=1 rpin=1", FNULL, true));
	eatts.Push (make_tuple(-3, "Lin2", "Beam", "BeamElastic", prms.CStr(), "ZERO", "cq=1",        FNULL, true));

	// Set geometry: nodes and elements
	dat.SetOnlyFrame  (true);
	dat.SetNodesElems (&mesh, &eatts);

	// Stage # 1 -----------------------------------------------------------
	dat.Nod(0)->Bry("ux", 0.0)->Bry("uy", 0.0);
	dat.Nod(3)->Bry("fy", P);
	dat.Nod(4)->Bry("ux", 0.0)->Bry("uy", 0.0)->Bry("wz", 0.0);
	dat.Ele(0)->EdgeBry("Qb",Qb,Qb,0);
	dat.Ele(1)->EdgeBry("Qb",Qb,Qb,0);
	sol.SolveWithInfo();

	//////////////////////////////////////////////////////////////////////////////////////// Output ////
	
	// Extra
	Matrix<double> coords;
	Vector<double> norm;
	Matrix<double> values;
	Array<String>  labels;

	// Elements
	std::ofstream of("beam03_elems.cal", std::ios::out);
	of << _8s<<"X" << _8s<<"Y" << _8s<<"M" << _8s<<"N" << _8s<<"V" << endl;
	for (size_t k=0; k<dat.NElems(); ++k)
	{
		dat.Ele(k)->OutExtra(coords, norm, values, labels);
		for (int i=0; i<coords.Rows(); ++i)
			of << _8s<<coords(i,0) << _8s<<coords(i,1) << _8s<<values(i,0) << _8s<<values(i,1) << _8s<<values(i,2) << endl;
	}
	cout << "\nFile <[1;31mbeam03_elems.cal[0m> created" << endl;
	of.close();

	// Elements
	std::ofstream pf("beam03_nodes.cal", std::ios::out);
	pf << _8s<<"X" << _8s<<"Y" << _8s<<"ux" << _8s<<"uy" << _8s<<"wz" << _8s<<"fx" << _8s<<"fy" << _8s<<"mz" << endl;
	for (size_t k=0; k<dat.NNodes(); ++k)
	{
		pf << _8s<<dat.Nod(k)->X() << _8s<<dat.Nod(k)->Y() << _8s<<dat.Nod(k)->Val("ux") << _8s<<dat.Nod(k)->Val("uy") << _8s<<dat.Nod(k)->Val("wz") << _8s<<dat.Nod(k)->Val("fx") << _8s<<dat.Nod(k)->Val("fy") << _8s<<dat.Nod(k)->Val("mz") << endl;
	}
	cout << "File <[1;31mbeam03_nodes.cal[0m> created" << endl;
	pf.close();

	//////////////////////////////////////////////////////////////////////////////////////// Check /////
	
	// Bending moments and shear forces
	Array<double> err_M;
	Array<double> err_V;
	dat.Ele(0)->OutExtra(coords, norm, values, labels);
    err_M.Push(fabs(values(0,0) -  0.0));
    err_M.Push(fabs(values(5,0) -  -Qb*l*l/8.0));
	err_V.Push(fabs(values(0,2) - -Qb*l/2.0));
	dat.Ele(2)->OutExtra(coords, norm, values, labels);
    err_M.Push(fabs(values(0,0) -  0.0));
	err_V.Push(fabs(values(0,2) - Qb*l/2.0));
	err_V.Push(fabs(values(10,2) - Qb*l/2.0));
	dat.Ele(3)->OutExtra(coords, norm, values, labels);
    err_M.Push(fabs(values(0,0) -  Qb*l*l/4.0));
    err_M.Push(fabs(values(10,0) - (P*l/2.0+Qb*l*l/2.0)));
	err_V.Push(fabs(values(10,2) - (P+Qb*l/2.0)));

	// Error summary
	double tol_M     = 1.0e-4;
	double tol_V     = 1.0e-4;
	double min_err_M = err_M[err_M.Min()];
	double max_err_M = err_M[err_M.Max()];
	double min_err_V = err_V[err_V.Min()];
	double max_err_V = err_V[err_V.Max()];
	cout << endl;
	cout << _4<< ""  << _8s<<"Min"     << _8s<<"Mean"                                                  << _8s<<"Max"                << _8s<<"Norm"       << endl;
	cout << _4<< "M" << _8s<<min_err_M << _8s<<err_M.Mean() << (max_err_M>tol_M?"[1;31m":"[1;32m") << _8s<<max_err_M << "[0m" << _8s<<err_M.Norm() << endl;
	cout << _4<< "V" << _8s<<min_err_V << _8s<<err_V.Mean() << (max_err_V>tol_V?"[1;31m":"[1;32m") << _8s<<max_err_V << "[0m" << _8s<<err_V.Norm() << endl;

	// Return error flag
	if (max_err_M>tol_M || max_err_V>tol_V) return 1;
	else return 0;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
