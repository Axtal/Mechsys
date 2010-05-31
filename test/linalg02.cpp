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
#include <cmath>

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    /////////////////////////////////////////////////////////////////////////////////////////// Correct values /////

    Vec3_t corLM, corV0M, corV1M, corV2M;
    corLM  =  3.2062308395514684e+01,  9.0230601363860448e+00, -2.1085368531900723e+01;
    corV0M =  2.4142214133881629e-01,  1.2938771667600124e-01,  9.6175577380369870e-01;
    corV1M =  9.6413836337045689e-01,  8.0603801105791670e-02, -2.5286408112785591e-01;
    corV2M =  1.1023867719052860e-01, -9.8831262565074463e-01,  1.0528811913323206e-01;

    double sq2h = sqrt(2.0)/2.0;
    Vec3_t corLM2, corV0M2, corV1M2, corV2M2;
    corLM2  =  1.0,   -3.0, -1.0;
    corV0M2 =  sq2h,  sq2h,  0.0;
    corV1M2 = -sq2h,  sq2h,  0.0;
    corV2M2 =   0.0,   0.0,  1.0;
    
    Vec3_t corxA;
    corxA = 4.0, -2.0, -2.0;

    double tol   = 1.0e-13;
    double error = 0.0;

    /////////////////////////////////////////////////////////////////////////////////////////// Eigen-problem /////

	cout << "\n--- eigen problem --------------------------------\n";
    Mat3_t M;
    M = 10.0,   4.0,  5.0,
         4.0, -20.0,  6.0,
         5.0,   6.0, 30.0;
    cout << "M  =\n" << PrintMatrix(M);
    Vec3_t L, V0, V1, V2;
    Eig (M, L, V0, V1, V2);
    cout << "M  =\n" << PrintMatrix(M);
    cout << "L  ="   << PrintVector(L);
    cout << "V0 ="   << PrintVector(V0);
    cout << "V1 ="   << PrintVector(V1);
    cout << "V2 ="   << PrintVector(V2);
    error += CompareVectors (L,  corLM);
    error += CompareVectors (V0, corV0M);
    error += CompareVectors (V1, corV1M);
    error += CompareVectors (V2, corV2M);

    cout << endl;
    M = -1.0,   2.0,   0.0,
         2.0,  -1.0,   0.0,
         0.0,   0.0,  -1.0;
    cout << "M  =\n" << PrintMatrix(M);
    Eig (M, L, V0, V1, V2);
    cout << "M  =\n" << PrintMatrix(M);
    cout << "L  ="   << PrintVector(L);
    cout << "V0 ="   << PrintVector(V0);
    cout << "V1 ="   << PrintVector(V1);
    cout << "V2 ="   << PrintVector(V2);
    error += CompareVectors (L,  corLM2);
    error += CompareVectors (V0, corV0M2);
    error += CompareVectors (V1, corV1M2);
    error += CompareVectors (V2, corV2M2);

    /////////////////////////////////////////////////////////////////////////////////////////// Solver /////
    
    Mat3_t A;
    A = 1.0,  1.0,  1.0,
        1.0, -2.0,  2.0,
        1.0,  2.0, -1.0;
    Vec3_t b, x;
    b = 0.0, 4.0, 2.0;
    Sol (A, b, x); // x = inv(A)*b

	cout << "\n--- solve A*b = x --------------------------------\n";
    cout << "A =\n" << PrintMatrix(A);
    cout << "b ="   << PrintVector(b);
    cout << "x ="   << PrintVector(x);
    error += CompareVectors (x, corxA);

	cout << "\n--- error ----------------------------------------\n";
	cout << "error = " << (error>tol ? "[1;31m" : "[1;32m") << error << "[0m" << endl;

    if (error>tol) return 1;
    else return 0;
}
MECHSYS_CATCH
