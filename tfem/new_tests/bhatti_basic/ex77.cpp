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

/*  Bhatti (2005): Example 7.7, p510  *
 *  ================================  */

// STL
#include <iostream>

// MechSys
#include "mesh/mesh.h"
#include "fem/elems/quad4.h"
#include "fem/equilibelem.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "models/linelastic.h"
#include "util/maps.h"
#include "util/fatal.h"

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize   (8/*verts*/, 3/*cells*/);
    mesh.SetVert   (0, -100,  0.0,  5.0, 0);
    mesh.SetVert   (1, -100,  0.0, 12.0, 0);
    mesh.SetVert   (2,    0,  6.0,  0.0, 0);
    mesh.SetVert   (3,    0,  6.0,  5.0, 0);
    mesh.SetVert   (4,    0, 20.0,  0.0, 0);
    mesh.SetVert   (5,    0, 20.0, 12.0, 0);
    mesh.SetVert   (6, -200, 54.0,  0.0, 0);
    mesh.SetVert   (7, -200, 54.0, 12.0, 0);
    mesh.SetCell   (0,   -1, /*NVerts*/4, 0,3,5,1);
    mesh.SetCell   (1,   -1, /*NVerts*/4, 2,4,5,3);
    mesh.SetCell   (2,   -1, /*NVerts*/4, 4,6,7,5);
    mesh.SetBryTag (0, 2, -10);
    mesh.SetBryTag (2, 2, -10);
    //mesh.WriteVTU  ("ex77");
    //mesh.WriteMPY  ("ex77");
    //cout << mesh << endl;

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom active  h pse  nip", PROB("Equilib"), GEOM("Quad4"), TRUE,  4.0, TRUE, 4.0);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu pse", MODEL("LinElastic"), 3.0e+6, 0.2, TRUE);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);

    // check matrices
    {
        double tol   = 1.0e-12;
        double error = 0.0;
        Mat_t K0c(8,8),K1c(8,8),K2c(8,8);
        K0c =
          6.2620867124142228e+06, -1.6375545851528419e+04, -4.1192295695570810e+06,  1.2663755458515282e+06, -9.5173112913287571e+05, -1.1299126637554583e+06, -1.1911260137242672e+06, -1.2008733624454145e+05,
         -1.6375545851528419e+04,  9.0354023705552090e+06,  2.5163755458515282e+06, -3.6782595134123508e+06, -1.1299126637554585e+06,  2.2847785402370521e+05, -1.3700873362445417e+06, -5.5856207111665634e+06,
         -4.1192295695570810e+06,  2.5163755458515282e+06,  1.1262086712414224e+07, -3.7663755458515272e+06, -1.1911260137242670e+06, -1.3700873362445412e+06, -5.9517311291328743e+06,  2.6200873362445412e+06,
          1.2663755458515282e+06, -3.6782595134123508e+06, -3.7663755458515272e+06,  2.1535402370555203e+07, -1.2008733624454139e+05, -5.5856207111665625e+06,  2.6200873362445412e+06, -1.2271522145976292e+07,
         -9.5173112913287582e+05, -1.1299126637554585e+06, -1.1911260137242670e+06, -1.2008733624454151e+05,  2.5448378041172801e+06,  7.8602620087336248e+05, -4.0198066126013710e+05,  4.6397379912663752e+05,
         -1.1299126637554585e+06,  2.2847785402370524e+05, -1.3700873362445415e+06, -5.5856207111665625e+06,  7.8602620087336260e+05,  2.5506862133499691e+06,  1.7139737991266376e+06,  2.8064566437928881e+06,
         -1.1911260137242675e+06, -1.3700873362445412e+06, -5.9517311291328752e+06,  2.6200873362445412e+06, -4.0198066126013710e+05,  1.7139737991266372e+06,  7.5448378041172791e+06, -2.9639737991266372e+06,
         -1.2008733624454128e+05, -5.5856207111665634e+06,  2.6200873362445412e+06, -1.2271522145976292e+07,  4.6397379912663740e+05,  2.8064566437928881e+06, -2.9639737991266372e+06,  1.5050686213349968e+07;
        K1c =
          5.7276021655606003e+06,  5.5929095354523242e+05, -7.7191756898358394e+05, -1.1705378973105133e+06, -1.4602252881592731e+06, -1.3294621026894865e+06, -3.4954593084177431e+06,  1.9407090464547675e+06,
          5.5929095354523230e+05,  9.6590115263709389e+06,  7.9462102689486797e+04,  2.7462451973454412e+06, -1.3294621026894865e+06, -3.6391023402025844e+06,  6.9070904645476758e+05, -8.7661543835137952e+06,
         -7.7191756898358418e+05,  7.9462102689486812e+04,  6.3632989870764920e+06, -2.7414425427872855e+06, -1.0061561299336362e+06,  2.4144254278728593e+05, -4.5852252881592717e+06,  2.4205378973105131e+06,
         -1.1705378973105133e+06,  2.7462451973454407e+06, -2.7414425427872855e+06,  7.3973978344393987e+06,  1.4914425427872858e+06, -5.2545406915822560e+06,  2.4205378973105131e+06, -4.8891023402025830e+06,
         -1.4602252881592731e+06, -1.3294621026894867e+06, -1.0061561299336362e+06,  1.4914425427872858e+06,  3.2382989870764930e+06,  1.0085574572127139e+06, -7.7191756898358383e+05, -1.1705378973105131e+06,
         -1.3294621026894865e+06, -3.6391023402025839e+06,  2.4144254278728599e+05, -5.2545406915822560e+06,  1.0085574572127139e+06,  6.1473978344393987e+06,  7.9462102689486521e+04,  2.7462451973454407e+06,
         -3.4954593084177431e+06,  6.9070904645476746e+05, -4.5852252881592717e+06,  2.4205378973105131e+06, -7.7191756898358383e+05,  7.9462102689486695e+04,  8.8526021655605994e+06, -3.1907090464547677e+06,
          1.9407090464547672e+06, -8.7661543835137952e+06,  2.4205378973105131e+06, -4.8891023402025830e+06, -1.1705378973105131e+06,  2.7462451973454417e+06, -3.1907090464547677e+06,  1.0909011526370937e+07;
        K2c =
          6.1928104575163396e+06,  1.8749999999999998e+06,  8.9052287581699330e+05, -6.2499999999999988e+05, -3.0964052287581693e+06, -1.8749999999999995e+06, -3.9869281045751632e+06,  6.2500000000000000e+05,
          1.8749999999999998e+06,  1.2393790849673200e+07,  6.2499999999999977e+05,  5.3145424836601298e+06, -1.8749999999999995e+06, -6.1968954248366002e+06, -6.2499999999999977e+05, -1.1511437908496730e+07,
          8.9052287581699330e+05,  6.2500000000000000e+05,  6.1928104575163387e+06, -1.8749999999999998e+06, -3.9869281045751637e+06, -6.2500000000000000e+05, -3.0964052287581689e+06,  1.8749999999999995e+06,
         -6.2499999999999977e+05,  5.3145424836601298e+06, -1.8749999999999995e+06,  1.2393790849673200e+07,  6.2499999999999977e+05, -1.1511437908496730e+07,  1.8749999999999995e+06, -6.1968954248365993e+06,
         -3.0964052287581689e+06, -1.8749999999999995e+06, -3.9869281045751637e+06,  6.2500000000000000e+05,  6.1928104575163387e+06,  1.8749999999999998e+06,  8.9052287581699318e+05, -6.2500000000000000e+05,
         -1.8749999999999995e+06, -6.1968954248365993e+06, -6.2499999999999988e+05, -1.1511437908496730e+07,  1.8749999999999995e+06,  1.2393790849673200e+07,  6.2499999999999977e+05,  5.3145424836601298e+06,
         -3.9869281045751632e+06, -6.2500000000000000e+05, -3.0964052287581693e+06,  1.8749999999999995e+06,  8.9052287581699342e+05,  6.2500000000000000e+05,  6.1928104575163396e+06, -1.8750000000000000e+06,
          6.2499999999999988e+05, -1.1511437908496730e+07,  1.8749999999999995e+06, -6.1968954248366002e+06, -6.2499999999999988e+05,  5.3145424836601289e+06, -1.8749999999999995e+06,  1.2393790849673202e+07;
        Mat_t K0,K1,K2;
        dom.Eles[0]->CalcK(K0);
        dom.Eles[1]->CalcK(K1);
        dom.Eles[2]->CalcK(K2);
        K0c /= 1.0e+6;  K0 /= 1.0e+6;
        K1c /= 1.0e+6;  K1 /= 1.0e+6;
        K2c /= 1.0e+6;  K2 /= 1.0e+6;
        error += CompareMatrices (K0,K0c);
        error += CompareMatrices (K1,K1c);
        error += CompareMatrices (K2,K2c);
        cout << "\n[1;37m--- Matrices: Error ----------------------------------------------------------[0m\n";
        cout << "error (K) = " << (error>tol ? "[1;31m" : "[1;32m") << error << "[0m" << endl;
    }

    // solver
    FEM::Solver sol(dom);

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set( -10, "qn",  -50.0);
    bcs.Set(-100, "ux",    0.0);
    bcs.Set(-200, "ux uy", 0.0,0.0);
    dom.SetBCs (bcs);
    sol.Solve  ();

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults (cout, Util::_12_6);

    //////////////////////////////////////////////////////////////////////////////////////// Check /////

    // correct solution
    Table nod_sol;
    nod_sol.Set("                   ux                       uy", /*NRows*/8,
                 0.000000000000000e+00,  -1.831554164758489e-02,
                 0.000000000000000e+00,  -1.832042377879234e-02,
                 2.759152276950652e-03,  -1.664860383229998e-02,
                 1.145517124108189e-03,  -1.646340822298832e-02,
                 3.050033425646983e-03,  -1.135663521522087e-02,
                -2.101277392742003e-03,  -1.162541210354794e-02,
                 0.000000000000000e+00,   0.000000000000000e+00,
                 0.000000000000000e+00,   0.000000000000000e+00);

    Table ele_sol;
    ele_sol.Set("                   sx                      sy                      sz                     sxy                       ex                       ey                      ez                     exy", /*NRows*/3,
                -1.045712337020673e+02,  2.854398647071007e+01,  0.000000000000000e+00,  1.669780828315312e+02,  -3.676001033206977e-05,   1.648607773704118e-05,  5.068483148757149e-06,  6.679123313261246e-05,
                -2.208482551524299e+01, -1.916660257635668e+01,  0.000000000000000e+00, -4.365550871990912e+01,  -6.083834999990551e-06,  -4.916545824436026e-06,  2.750095206106645e-06, -1.746220348796365e-05,
                -5.060031876334184e+01, -4.371717479355189e+01,  0.000000000000000e+00,  1.541666666666671e+02,  -1.395229460154382e-05,  -1.119903701362784e-05,  6.287832903792917e-06,  6.166666666666685e-05);

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("ux uy", 1.0e-15, 1.0e-15);
    ele_tol.Set("sx sy sz sxy  ex ey ez exy", 1.0e-11,1.0e-11,1.0e-15,1.0e-11, 1.0e-15,1.0e-15,1.0e-15,1.0e-15);

    // return error flag
    return dom.CheckError (cout, nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH
