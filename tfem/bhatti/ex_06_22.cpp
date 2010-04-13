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

/*  Bhatti (2005): Example 6.22 p449 *
 *  ================================ */

// STL
#include <iostream>

// MechSys
#include <mechsys/mesh/mesh.h>
#include <mechsys/fem/elems/quad8.h>
#include <mechsys/fem/flowelem.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solver.h>
#include <mechsys/models/linflow.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize   (13/*verts*/, 2/*cells*/);
    mesh.SetVert   ( 0,    0,  0.0,   0.03);
    mesh.SetVert   ( 1,    0,  0.015, 0.03);
    mesh.SetVert   ( 2,    0,  0.03,  0.03);
    mesh.SetVert   ( 3,    0,  0.03,  0.0225);
    mesh.SetVert   ( 4,    0,  0.03,  0.015);
    mesh.SetVert   ( 5,    0,  0.045, 0.015);
    mesh.SetVert   ( 6,    0,  0.06,  0.015);
    mesh.SetVert   ( 7,    0,  0.06,  0.0075);
    mesh.SetVert   ( 8, -100,  0.06,  0.0);
    mesh.SetVert   ( 9, -100,  0.03,  0.0);
    mesh.SetVert   (10, -100,  0.0,   0.0);
    mesh.SetVert   (11,    0,  0.0,   0.015);
    mesh.SetVert   (12,    0,  0.015, 0.0075);
    mesh.SetCell   ( 0,   -1, Array<int>(10,4,2,0, 12,3,1,11));
    mesh.SetCell   ( 1,   -1, Array<int>(10,8,6,4,  9,7,5,12));
    mesh.SetBryTag (0, 1, -30);
    mesh.SetBryTag (0, 2, -30);
    mesh.SetBryTag (0, 3, -40);
    mesh.SetBryTag (1, 1, -20);
    mesh.SetBryTag (1, 2, -30);
    mesh.WriteMPY  ("ex622");

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom", PROB("Flow"), GEOM("Quad8"));

    // models
    Dict mdls;
    mdls.Set(-1, "name k", MODEL("LinFlow"), 45.0);

    // initial values
    Dict inis;
    inis.Set(-1, "vx vy", 0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);

    // solver
    FEM::Solver sol(dom);
    //sol.CteTg = true;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set(-100, "H",           110.0);
    bcs.Set( -20, "flux",        0.0);
    bcs.Set( -30, "conv h Tinf", 1.0, 55.0, 20.0);
    bcs.Set( -40, "flux",        8000.0);
    bcs.Set(  -1, "s",           5.0e+6);
    dom.SetBCs (bcs);
    sol.Solve  ();
    //sol.Solve  (10,1);

    // check matrices
    {
        double tol   = 1.0e-11;
        double error = 0.0;
        Mat_t K0c(8,8),K1c(8,8);
        K0c =
          3.8515873015873034e+01,  2.3194444444444454e+01,  2.1468253968253972e+01,  2.2027777777777782e+01, -2.0317460317460323e+01, -3.0912698412698425e+01, -1.4682539682539684e+01, -3.9293650793650805e+01,
          2.3194444444444454e+01,  8.4197301587301610e+01,  3.3583611111111111e+01,  3.0456349206349213e+01, -2.6984126984127005e+01, -8.2643412698412718e+01, -2.8015873015873023e+01, -3.3650793650793659e+01,
          2.1468253968253975e+01,  3.3583611111111111e+01,  5.8643492063492083e+01,  2.3139444444444450e+01, -1.9365079365079371e+01, -6.8119603174603185e+01, -2.0524920634920647e+01, -2.8412698412698422e+01,
          2.2027777777777782e+01,  3.0456349206349213e+01,  2.3139444444444454e+01,  5.7366825396825405e+01, -1.6507936507936513e+01, -3.6150793650793659e+01, -3.3382063492063487e+01, -4.6674603174603185e+01,
         -2.0317460317460323e+01, -2.6984126984127005e+01, -1.9365079365079371e+01, -1.6507936507936513e+01,  9.5079365079365090e+01,  1.6349206349206362e+01, -5.0793650793650738e+00, -2.3174603174603174e+01,
         -3.0912698412698422e+01, -8.2643412698412718e+01, -6.8119603174603185e+01, -3.6150793650793659e+01,  1.6349206349206362e+01,  1.5631301587301590e+02,  3.6507936507936645e+00,  4.2063492063492077e+01,
         -1.4682539682539685e+01, -2.8015873015873023e+01, -2.0524920634920647e+01, -3.3382063492063494e+01, -5.0793650793650720e+00,  3.6507936507936662e+00,  9.5959365079365071e+01,  3.1746031746031838e+00,
         -3.9293650793650805e+01, -3.3650793650793659e+01, -2.8412698412698422e+01, -4.6674603174603185e+01, -2.3174603174603181e+01,  4.2063492063492077e+01,  3.1746031746031860e+00,  1.2596825396825400e+02;
        K1c =
          4.3051587301587311e+01,  2.5313492063492074e+01,  2.6646825396825406e+01,  4.9623015873015888e+01,  1.6634920634920629e+01, -4.8015873015873034e+01, -2.1269841269841269e+01, -9.1984126984127016e+01,
          2.5313492063492074e+01,  1.1757539682539687e+02,  4.9623015873015888e+01,  6.2599206349206369e+01, -1.2888888888888921e+01, -1.4468253968253967e+02, -4.2222222222222229e+01, -5.5317460317460345e+01,
          2.6646825396825406e+01,  4.9623015873015881e+01,  7.0676349206349244e+01,  5.1948968253968275e+01, -1.1269841269841287e+01, -1.0396825396825398e+02, -2.7350317460317481e+01, -5.6031746031746046e+01,
          4.9623015873015881e+01,  6.2599206349206362e+01,  5.1948968253968275e+01,  1.7377158730158735e+02, -3.2222222222222236e+01, -9.0634920634920647e+01, -8.5445555555555558e+01, -1.2936507936507942e+02,
          1.6634920634920622e+01, -1.2888888888888918e+01, -1.1269841269841285e+01, -3.2222222222222243e+01,  1.5625396825396828e+02,  1.2698412698412733e+01, -3.6507936507936492e+01, -9.2698412698412696e+01,
         -4.8015873015873034e+01, -1.4468253968253967e+02, -1.0396825396825398e+02, -9.0634920634920647e+01,  1.2698412698412730e+01,  2.5888888888888886e+02,  1.4603174603174629e+01,  1.0111111111111113e+02,
         -2.1269841269841269e+01, -4.2222222222222229e+01, -2.7350317460317481e+01, -8.5445555555555558e+01, -3.6507936507936492e+01,  1.4603174603174615e+01,  1.3389587301587301e+02,  6.5396825396825392e+01,
         -9.1984126984126988e+01, -5.5317460317460345e+01, -5.6031746031746046e+01, -1.2936507936507942e+02, -9.2698412698412668e+01,  1.0111111111111113e+02,  6.5396825396825392e+01,  2.5888888888888897e+02;
        Mat_t K0,K1;
        dom.Eles[0]->CalcK(K0);
        dom.Eles[1]->CalcK(K1);
        error += CompareMatrices (K0,K0c);
        error += CompareMatrices (K1,K1c);
        cout << "\n[1;37m--- Matrices: Error ----------------------------------------------------------[0m\n";
        cout << "error (K) = " << (error>tol ? "[1;31m" : "[1;32m") << error << "[0m" << endl;
    }

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults ("%11.6g");

    //////////////////////////////////////////////////////////////////////////////////////// Check /////

    // correct solution
    Table nod_sol;
    nod_sol.Set("H", dom.Nods.Size(),
                 1.564405024662020e+02,
                 1.507560541872985e+02,
                 1.491964629456365e+02,
                 1.442245542836662e+02,
                 1.338432701060946e+02,
                 1.240019529443106e+02,
                 1.217463572762217e+02,
                 1.191481315065258e+02,
                 1.100000000000000e+02,
                 1.100000000000000e+02,
                 1.100000000000000e+02,
                 1.446754222244301e+02,
                 1.291320079882026e+02);

    Table ele_sol;
    ele_sol.Set("gx gy", dom.Eles.Size(),
                -2.552961113487533e+02,  9.610687199598173e+02,
                -2.218639218150420e+02,  1.155327451435747e+03);

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("H", 1.0e-12);
    ele_tol.Set("gx gy", 1.0e-11, 1.0e-10);

    // return error flag
    return dom.CheckError (nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH
