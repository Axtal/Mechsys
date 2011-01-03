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
using Util::SQ2;

int main(int argc, char **argv) try
{
    /////////////////////////////////////////////////////////////////////////////////////////// Correct values /////

    double corDetM = -6100.0;
    Mat3_t corMi;
    corMi = 1.0426229508196722e-01,   1.4754098360655736e-02,  -2.0327868852459016e-02,
            1.4754098360655738e-02,  -4.5081967213114749e-02,   6.5573770491803270e-03,
           -2.0327868852459016e-02,   6.5573770491803270e-03,   3.5409836065573769e-02;
    Vec_t corTi(6);
    corTi = corMi(0,0), corMi(1,1), corMi(2,2), corMi(0,1)*SQ2, corMi(1,2)*SQ2, corMi(0,2)*SQ2;

    Vec3_t corL, corV0, corV1, corV2;
    corL  =  3.2062308395514684e+01,  9.0230601363860448e+00, -2.1085368531900723e+01;
    corV0 =  2.4142214133881629e-01,  1.2938771667600124e-01,  9.6175577380369870e-01;
    corV1 =  9.6413836337045689e-01,  8.0603801105791670e-02, -2.5286408112785591e-01;
    corV2 =  1.1023867719052860e-01, -9.8831262565074463e-01,  1.0528811913323206e-01;
    Mat3_t corP0m, corP1m, corP2m;
    Vec_t  corP0,  corP1,  corP2;
    Dyad    (corV0, corV0, corP0m);
    Dyad    (corV1, corV1, corP1m);
    Dyad    (corV2, corV2, corP2m);
    Mat2Ten (corP0m, corP0, /*ncp*/6);
    Mat2Ten (corP1m, corP1, /*ncp*/6);
    Mat2Ten (corP2m, corP2, /*ncp*/6);

    double tol   = 1.0e-13;
    double error = 0.0;

    //////////////////////////////////////////////////////////////////////////////////// Determinat and Inverse ///
    
	cout << "\n--- determinat and inverse -----------------------\n";
    Mat3_t M, Mi, MxMi, MxM, MxMblitz;
    M = 10.0,   4.0,  5.0,
         4.0, -20.0,  6.0,
         5.0,   6.0, 30.0;
    Inv  (M, Mi);
    Mult (M,Mi, MxMi);
    Mult (M,M, MxM);
    MxMblitz = blitz::product (M,M);
    cout << "M  =\n"    << PrintMatrix(M);
    cout << "Mi =\n"    << PrintMatrix(Mi);
    cout << "M*Mi = \n" << PrintMatrix(MxMi);
    error += fabs(Det(M)-corDetM);
    error += CompareMatrices (Mi,  corMi);
    error += CompareMatrices (MxM, MxMblitz);

    ////////////////////////////////////////////////////////////////////////////////////////////// Tensors ////////
    
	cout << "\n--- tensors --------------------------------------\n";
    Vec_t T, Ti;
    Mat2Ten (M, T, /*ncp*/6);
    Inv     (T, Ti);
    cout << "T  = " << PrintVector(T);
    cout << "Ti = " << PrintVector(Ti);
    error += CompareVectors (Ti, corTi);

    //////////////////////////////////////////////////////////////////////////////////////// Eigen-projectors /////

	cout << "\n--- eigen-projectors -----------------------------\n";
    Vec3_t L;
    Vec3_t v0,v1,v2;
    Vec_t  P0,P1,P2;
    Mat3_t p0,p1,p2, p0xp0, p0xp1;
    EigenProj (T, L, v0, v1, v2, P0, P1, P2);
    Ten2Mat   (P0, p0);
    Ten2Mat   (P1, p1);
    Ten2Mat   (P2, p2);
    cout << "L     = " << PrintVector(L);
    cout << "P0    = " << PrintVector(P0);
    cout << "P1    = " << PrintVector(P1);
    cout << "P2    = " << PrintVector(P2);
    cout << "detP0 = " << Det(p0) << endl;
    cout << "detP1 = " << Det(p1) << endl;
    cout << "detP2 = " << Det(p2) << endl;
    error += CompareVectors  (L,  corL);
    error += CompareVectors  (P0, corP0);
    error += CompareVectors  (P1, corP1);
    error += CompareVectors  (P2, corP2);
    error += CompareMatrices (p0, corP0m);
    error += CompareMatrices (p1, corP1m);
    error += CompareMatrices (p2, corP2m);
    Mat3_t * P[3] = {&p0,&p1,&p2};
    Mat3_t PxP, zero;
    set_to_zero(zero);
    for (size_t i=0; i<3; ++i)
    for (size_t j=0; j<3; ++j)
    {
        Mult ((*P[i]), (*P[j]), PxP);
        //cout << "P" << i << "*P" << j << " =\n" << PrintMatrix(PxP);
        if (i==j) error += CompareMatrices (PxP, (*P[i]));
        else      error += CompareMatrices (PxP, zero);
    }
    //Mat3_t P0i;
    //Inv (p0, P0i);
    //cout << "P0i = " << PrintMatrix(P0i);

    /////////////////////////////////////////////////////////////////////////////////////////// Error /////////////
    
	cout << "\n--- error ----------------------------------------\n";
	cout << "error = " << (error>tol ? "[1;31m" : "[1;32m") << error << "[0m" << endl;

    if (error>tol) return 1;
    else return 0;
}
MECHSYS_CATCH
