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

/*  Bhatti (2005): Example 4.8, p260  *
 *  ================================  */

// STL
#include <iostream>

// MechSys
#include <mechsys/mesh/mesh.h>
#include <mechsys/fem/beam.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solver.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    double L = 2.0;
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize  (4/*nodes*/, 3/*cells*/);
    mesh.SetVert  (0, -100,  0.0, 0.0);
    mesh.SetVert  (1, -200,    L, 0.0);
    mesh.SetVert  (2, -300, 2.*L, 0.0);
    mesh.SetVert  (3, -400, 3.*L, 0.0);
    mesh.SetCell  (0,   -1, Array<int>(0, 1));
    mesh.SetCell  (1,   -1, Array<int>(1, 2));
    mesh.SetCell  (2,   -2, Array<int>(2, 3));

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob E A Izz fra", PROB("Beam"), 420.0e+6, 0.0001, 4.0e-4, TRUE);
    prps.Set(-2, "prob E A Izz fra", PROB("Beam"), 210.0e+6, 0.0001, 4.0e-4, TRUE);

    // domain
    FEM::Domain dom(mesh, prps, /*mdls*/Dict(), /*inis*/Dict());

    // check matrices
    bool err_flag = false;
    {
        double tol   = 1.0e-15;
        double error = 0.0;
        Mat_t K0c(6,6),K1c(6,6),K2c(6,6);
        K0c =
           21000.0,       0.0,       0.0,  -21000.0,       0.0,       0.0,
               0.0,  252000.0,  252000.0,       0.0, -252000.0,  252000.0,
               0.0,  252000.0,  336000.0,       0.0, -252000.0,  168000.0,
          -21000.0,       0.0,       0.0,   21000.0,       0.0,       0.0,
               0.0, -252000.0, -252000.0,       0.0,  252000.0, -252000.0,
               0.0,  252000.0,  168000.0,       0.0, -252000.0,  336000.0;
        K1c =
           21000.0,       0.0,       0.0,  -21000.0,       0.0,       0.0,
               0.0,  252000.0,  252000.0,       0.0, -252000.0,  252000.0,
               0.0,  252000.0,  336000.0,       0.0, -252000.0,  168000.0,
          -21000.0,       0.0,       0.0,   21000.0,       0.0,       0.0,
               0.0, -252000.0, -252000.0,       0.0,  252000.0, -252000.0,
               0.0,  252000.0,  168000.0,       0.0, -252000.0,  336000.0;
        K2c =
           10500.0,       0.0,       0.0,  -10500.0,       0.0,       0.0,
               0.0,  126000.0,  126000.0,       0.0, -126000.0,  126000.0,
               0.0,  126000.0,  168000.0,       0.0, -126000.0,   84000.0,
          -10500.0,       0.0,       0.0,   10500.0,       0.0,       0.0,
               0.0, -126000.0, -126000.0,       0.0,  126000.0, -126000.0,
               0.0,  126000.0,   84000.0,       0.0, -126000.0,  168000.0;
        Mat_t K0,K1,K2;
        dom.Eles[0]->CalcK(K0);
        dom.Eles[1]->CalcK(K1);
        dom.Eles[2]->CalcK(K2);
        error += CompareMatrices (K0,K0c);
        error += CompareMatrices (K1,K1c);
        error += CompareMatrices (K2,K2c);
        cout << "\n[1;37m--- Matrices: Error ----------------------------------------------------------[0m\n";
        cout << "error (K) = " << (error>tol ? "[1;31m" : "[1;32m") << error << "[0m" << endl;
        if (error>tol) err_flag = true;
    }

    // solver
    FEM::Solver sol(dom);
    //sol.CteTg  = true;
    //sol.Scheme = FEM::Solver::FE_t;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set(-100, "ux uy wz", 0.0, 0.0, 0.0);
    bcs.Set(-200, "fy",     -18.0);
    bcs.Set(-400, "ux uy",    0.0, 0.0);
    bcs.Set(  -2, "qn",     -10.0);
    dom.SetBCs (bcs);
    sol.Solve  (1);

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults ("%11.6g");
    dom.WriteMPY     ("ex48_res",/*sf*/0.01);

    //////////////////////////////////////////////////////////////////////////////////////// Check /////

    // correct solution
    Table nod_sol;
    nod_sol.Set("ux  uy  wz  Rux  Ruy  Rwz", dom.Nods.Size(),
                 0.0,   0.000000000000000e+00,   0.000000000000000e+00,   0.0,  2.060714285714285e+01,  3.164285714285713e+01,
                 0.0,  -2.131519274376416e-04,  -1.313775510204081e-04,   0.0,  0.000000000000000e+00,  0.000000000000000e+00,
                 0.0,  -3.412698412698411e-04,   1.360544217687071e-05,   0.0,  0.000000000000000e+00,  0.000000000000000e+00,
                 0.0,   0.000000000000000e+00,   2.689909297052154e-04,   0.0,  1.739285714285714e+01,  0.000000000000000e+00);

    // error tolerance
    Table ele_sol;
    SDPair nod_tol, ele_tol;
    nod_tol.Set("ux uy wz  Rux Ruy Rwz", 1.0e-15,1.0e-15,1.0e-15, 1.0e-15,1.0e-13,1.0e-14);
    err_flag = err_flag || dom.CheckError (nod_sol, ele_sol, nod_tol, ele_tol);

    // bending moment
    Mat_t Mmat(7,2);
    Mmat = 0.0, -3.1642857142857132e+01,
           1.0, -1.1035714285714279e+01,
           2.0,  9.5714285714285747e+00,
           3.0,  1.2178571428571418e+01,
           4.0,  1.4785714285714272e+01,
           5.0,  1.2392857142857141e+01,
           6.0,  0.0000000000000000e+00;
    Array<double> err_M;
    double tol_M = 1.0e-13;
    for (size_t i=0; i<7; ++i)
    {
        double x  = Mmat(i,0);
        double Mc = Mmat(i,1);
        double P,V,M;
        int    idx;
        double r;
        if      (x<L)    { idx = 0;  r = x/L;        }
        else if (x<2.*L) { idx = 1;  r = (x-L)/L;    }
        else             { idx = 2;  r = (x-2.*L)/L; }
        FEM::Beam const * ele = static_cast<FEM::Beam const *>(dom.Eles[idx]);
        SDPair res;
        ele->CalcRes (r, P,V,M);
        err_M.Push (fabs(Mc-M));
    }
    double max_err_M = err_M[err_M.Max()];
    cout << Util::_4<< "M" << Util::_8s<<err_M[err_M.Min()] << Util::_8s<<err_M.Mean();
    cout << (max_err_M>tol_M ? "[1;31m" : "[1;32m") << Util::_8s<<max_err_M << "[0m" << Util::_8s<<err_M.Norm() << "\n";
    if (max_err_M>tol_M) err_flag = false;

    // return error flag
    return err_flag;
}
MECHSYS_CATCH
