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

// Std lib
#include <iostream>

// MechSys
#include "linalg/lawrap.h"
#include "linalg/laexpr.h"
#include "linalg/matrix.h"
#include "util/exception.h"
#include "util/numstreams.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_8s;

int main(int argc, char **argv) try
{
	Matrix<double> A(5,5);
	A = 12, 28, 22, 20,  8,
	     0,  3,  5, 17, 28,
	    56,  0, 23,  1,  0,
	    12, 29, 27, 10,  1,
	     9,  4, 13,  8, 22;
	LinAlg::Geinv(A);

	Matrix<double> invA(5,5);
	Matrix<double> b(5,5);
	invA =  115723/2.0 ,  -62128.0,  -8266.0,  -57806.0,   60659.0,
	          133387.0 , -146329.0, -23691.0, -130580.0,  143668.0,
	       -273625/2.0 ,  146442.0,  22992.0,  136648.0, -142842.0,
	       -187113/2.0 ,  111002.0,  17781.0,   94232.0, -111538.0,
	        133883/2.0 ,  -74877.0, -12363.0,  -67623.0,   77834.0;
	invA = invA / 83701.0;

	double tol   = 1.0e-13;
	double error = 0.0;
	for (size_t i=0; i<5; ++i)
	for (size_t j=0; j<5; ++j)
		error += fabs(A(i,j)-invA(i,j));

	cout << "Error = " << (error>tol ? "[1;31m":"[1;32m") << _8s<<error << "[0m" << endl;

	if (error>tol) return 1;
	else return 0;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
