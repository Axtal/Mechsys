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

// MechSys
#include <mechsys/mesh/mesh.h>
#include <mechsys/fem/equilibelem.h>
#include <mechsys/fem/beam.h>
#include <mechsys/fem/rod.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solver.h>
#include <mechsys/models/linelastic.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    // elements properties
    double gra = 9.806; // m/s2
    double wid = 1.0;   // m
    double hei = 0.3;   // m
    double rho = 2.5;   // g/m3
    double E   = 38e+6; // kN/m2
    double A   = 0.3;   // m2
    double Izz = wid*pow(hei,3.0)/12.0; // m4
    //cout << "Izz = " << Izz << endl;
    Dict prps;
    prps.Set (-1, "prob fra  rho E A Izz", PROB("Beam"), 1.,  rho, E, A, Izz);

    // boundary conditions
    Dict bcs;

    bool   dead_load  = false;
    double sf         = 1.0;
    double L          = 1.0;
    bool   with_error = false;
    for (size_t tst=1; tst<=7; ++tst)
    {
        bcs.clear();
        printf("\n\n\n%s========================================= Test %zd ===========================================%s\n",TERM_YELLOW_BLUE,tst,TERM_RST);
        Mesh::Generic mesh(/*NDim*/2);
        switch (tst)
        {
            case 1:
            {
                mesh.SetSize (2, 1);
                mesh.SetVert (0, -100, 0.0, 0.0);
                mesh.SetVert (1, -200,   L, 0.0);
                mesh.SetCell (0,    -1, Array<int>(0,1));

                bcs.Set(-100, "ux uy", 0.0, 0.0);
                bcs.Set(-200, "uy",    0.0);
                bcs.Set(  -1, "qn",   -1.0);
                dead_load = false;
                //bcs.Set(  -1, "gravity", gra);
                break;
            }
            case 2:
            {
                mesh.SetSize (2, 1);
                mesh.SetVert (0, -100, 0.0, 0.0);
                mesh.SetVert (1, -200,   L, 0.0);
                mesh.SetCell (0,    -1, Array<int>(1,0));

                bcs.Set(-100, "ux uy", 0.0, 0.0);
                bcs.Set(-200, "uy",    0.0);
                bcs.Set(  -1, "qn",   -1.0);
                dead_load = false;
                //bcs.Set(  -1, "gravity", gra);
                break;
            }
            case 3:
            {
                mesh.SetSize (2, 1);
                mesh.SetVert (0, -100, 0.0, 0.0);
                mesh.SetVert (1, -200,   L, 0.0);
                mesh.SetCell (0,    -1, Array<int>(0,1));

                bcs.Set(-100, "ux uy wz", 0.0, 0.0, 0.0);
                bcs.Set(  -1, "qn",      -1.0);
                break;
            }
            case 4:
            {
                mesh.SetSize (3, 2);
                mesh.SetVert (0, -100,  0.0, 0.0);
                mesh.SetVert (1,    0, L/2., 0.0);
                mesh.SetVert (2, -200,    L, 0.0);
                mesh.SetCell (0,   -1, Array<int>(0,1));
                mesh.SetCell (1,   -1, Array<int>(1,2));

                bcs.Set(-100, "ux uy wz", 0.0, 0.0, 0.0);
                bcs.Set(  -1, "qn",      -1.0);
                break;
            }
            case 5:
            {
                mesh.SetSize (6, 5);
                mesh.SetVert (0, -100,  0.0, 2.5);
                mesh.SetVert (1, -100,  2.5, 0.0);
                mesh.SetVert (2, -100, 10.0, 5.0);
                mesh.SetVert (3,    0,  7.5, 2.5);
                mesh.SetVert (4,    0,  2.5, 2.5);
                mesh.SetVert (5, -200,  5.0, 2.5);
                mesh.SetCell (0,   -3, Array<int>(0,4));
                mesh.SetCell (1,   -1, Array<int>(4,1));
                mesh.SetCell (2,   -1, Array<int>(3,2));
                mesh.SetCell (3,   -2, Array<int>(3,5));
                mesh.SetCell (4,   -2, Array<int>(5,4));
                //mesh.Refine  ();
                //mesh.Refine  ();

                dead_load = false;
                bcs.Set (-100, "ux uy", 0.0, 0.0);
                if (dead_load)
                {
                    bcs.Set (-1, "gravity", gra);
                    bcs.Set (-2, "gravity", gra);
                    bcs.Set (-3, "gravity", gra);
                }
                else
                {
                    //bcs.Set (-200, "mz",   12.5);
                    bcs.Set (  -2, "qn",   -9.2);
                }

                prps.Set (-2, "prob fra  rho E A Izz", PROB("Beam"), 1.,  rho, E, A, Izz);
                prps.Set (-3, "prob fra  rho E A Izz", PROB("Beam"), 1.,  rho, E, A, Izz);
                //prps.Set (-3, "prob fra E A",         PROB("Rod"),  1.0, 1.0, 1.0);

                sf = 0.005;
                break;
            }
            case 6:
            {
                mesh.SetSize  (6, 5);
                mesh.SetVert  (0, -100, 0.0, 0.0);
                mesh.SetVert  (1,    0, 0.0,   L);
                mesh.SetVert  (2,    0,   L, 2*L);
                mesh.SetVert  (3,    0, 2*L, 2*L);
                mesh.SetVert  (4,    0, 3*L,   L);
                mesh.SetVert  (5, -100, 3*L, 0.0);
                mesh.SetCell  (0,   -1, Array<int>(0,1));
                mesh.SetCell  (1,   -1, Array<int>(1,2));
                mesh.SetCell  (2,   -2, Array<int>(2,3));
                mesh.SetCell  (3,   -1, Array<int>(3,4));
                mesh.SetCell  (4,   -1, Array<int>(4,5));
                //mesh.SetCell  (0,   -1, Array<int>(1,0));
                //mesh.SetCell  (1,   -1, Array<int>(2,1));
                //mesh.SetCell  (2,   -2, Array<int>(3,2));
                //mesh.SetCell  (3,   -1, Array<int>(3,4));
                //mesh.SetCell  (4,   -1, Array<int>(5,4));
                
                bcs.Set (-100, "ux uy wz", 0.0, 0.0, 0.0);
                bcs.Set (  -2, "qn",      -1.0);

                prps.Set (-2, "prob fra  rho E A Izz", PROB("Beam"), 1.,  rho, E, A, Izz);

                sf = 0.5;
                break;
            }
            case 7:
            {
                prps.Set (-1, "prob fra  rho E A Izz", PROB("Beam"), 1.,  1.,1.,1.,1.);

                mesh.SetSize (3, 2);
                mesh.SetVert (0, -100,  0.0, 0.0);
                mesh.SetVert (1, -300, L/2., 0.0);
                mesh.SetVert (2, -200,    L, 0.0);
                mesh.SetCell (0,   -1, Array<int>(0,1));
                mesh.SetCell (1,   -1, Array<int>(1,2));

                mesh.AddPin (-300);

                bcs.Set(-100, "ux uy wz", 0.0, 0.0, 0.0);
                bcs.Set(-200, "ux uy",    0.0, 0.0);
                bcs.Set(  -1, "qn",      -10.0);
                break;
            }
            default: throw new Fatal("main: Test = %d is not available",tst);
        }
        //mesh.WriteMPY ("beam01_mesh");

        // domain
        FEM::Beam::DrwAxial = false;
        FEM::Beam::DrwShear = false;
        FEM::Domain dom(mesh, prps, Dict(), Dict());

        // solver
        FEM::Solver sol(dom);
        sol.SetScheme ("FE");

        // run
        dom.SetBCs (bcs);
        //cout << dom << endl;
        sol.Solve  (/*NDiv*/1);
        //dom.PrintResults ("%15.6e",false);

        // output
        //dom.WriteMPY ("beam01_res", sf);
        //dom.WriteVTU ("beam01_res");

        // check
        cout << endl;
        double N,V,M;
        double errors = 0.0;
        double err_u  = 0.0;
        double tol    = 1.0e-15;
        double tol_u  = 1.0e-15;
        switch (tst)
        {
            case 1:
            {
                static_cast<FEM::Beam*>(dom.Eles[0])->CalcRes (0.5, N,V,M);
                double Mcor;
                if (dead_load) Mcor = rho*gra*A/8.;
                else           Mcor = 1.0/8.0;
                cout << "M(max) = " << M << "  => " << Mcor << endl;
                errors += fabs(M-Mcor);
                break;
            }
            case 2:
            {
                static_cast<FEM::Beam*>(dom.Eles[0])->CalcRes (0.5, N,V,M);
                double Mcor;
                if (dead_load) Mcor = rho*gra*A/8.;
                else           Mcor = 1.0/8.0;
                cout << "M(max) = " << M << "  => " << Mcor << endl;
                errors += fabs(M-Mcor);
                break;
            }
            case 3:
            {
                static_cast<FEM::Beam*>(dom.Eles[0])->CalcRes (0.0, N,V,M);
                double Mcor = -1.0/2.0;
                cout << "M(max) = " << M << "  => " << Mcor << endl;
                errors += fabs(M-Mcor);
                break;
            }
            case 4:
            {
                break;
            }
            case 5:
            {
                FEM::Beam const * e0 = static_cast<FEM::Beam const*>(dom.Eles[0]);
                FEM::Beam const * e1 = static_cast<FEM::Beam const*>(dom.Eles[1]);
                FEM::Beam const * e2 = static_cast<FEM::Beam const*>(dom.Eles[2]);
                FEM::Beam const * e3 = static_cast<FEM::Beam const*>(dom.Eles[3]);
                FEM::Beam const * e4 = static_cast<FEM::Beam const*>(dom.Eles[4]);

                if (dead_load)
                {
                    e0->CalcRes(0.0,N,V,M);  errors += fabs(M - (  0.0000));  printf("M: %g  =>  %g\n",M,  0.0000);
                    e0->CalcRes(1.0,N,V,M);  errors += fabs(M - (-10.3336));  printf("M: %g  =>  %g\n",M,-10.3336);
                    e1->CalcRes(0.0,N,V,M);  errors += fabs(M - (  4.6232));  printf("M: %g  =>  %g\n",M,  4.6232);
                    e1->CalcRes(1.0,N,V,M);  errors += fabs(M - (  0.0000));  printf("M: %g  =>  %g\n",M,  0.0000);
                    e2->CalcRes(0.0,N,V,M);  errors += fabs(M - (-11.8340));  printf("M: %g  =>  %g\n",M,-11.8340);
                    e2->CalcRes(1.0,N,V,M);  errors += fabs(M - (  0.0000));  printf("M: %g  =>  %g\n",M,  0.0000);
                    e3->CalcRes(0.0,N,V,M);  errors += fabs(M - (-11.8340));  printf("M: %g  =>  %g\n",M,-11.8340);
                    e3->CalcRes(1.0,N,V,M);  errors += fabs(M - (  9.5987));  printf("M: %g  =>  %g\n",M,  9.5987);
                    e4->CalcRes(0.0,N,V,M);  errors += fabs(M - (  9.5968));  printf("M: %g  =>  %g\n",M,  9.5968);
                    e4->CalcRes(1.0,N,V,M);  errors += fabs(M - (-14.9567));  printf("M: %g  =>  %g\n",M,-14.9567);
                    tol = 1.0e-1;

                    double ux = dom.Nods[4]->U("ux")*1000.0;
                    double uy = dom.Nods[4]->U("uy")*1000.0;
                    double wz = dom.Nods[4]->U("wz");
                    printf("ux, uy, uz: %g => %g,  %g => %g,  %g => %g  [mm, mm,  ]\n", ux,8.193e-3, uy,-9.111e-3, wz,4.834e-5);
                    err_u += fabs(ux-( 8.193e-3));
                    err_u += fabs(uy-(-9.111e-3));
                    err_u += fabs(wz-( 4.834e-2));
                    ux = dom.Nods[3]->U("ux")*1000.0;
                    uy = dom.Nods[3]->U("uy")*1000.0;
                    wz = dom.Nods[3]->U("wz");
                    printf("ux, uy, uz: %g => %g,  %g => %g,  %g => %g  [mm, mm,  ]\n", ux,2.377e-2, uy,-4.432e-2, wz,-6.469e-5);
                    err_u += fabs(ux-( 2.377e-2));
                    err_u += fabs(uy-(-4.432e-2));
                    err_u += fabs(wz-(-6.469e-5));
                    tol_u = 1.0e-1;
                }
                else
                {
                    e0->CalcRes(0.0,N,V,M);  errors += fabs(M - (  0.0000));  printf("M: %g  =>  %g\n",M,  0.0000);
                    e0->CalcRes(1.0,N,V,M);  errors += fabs(M - ( -9.3832));  printf("M: %g  =>  %g\n",M, -9.3832);
                    e1->CalcRes(0.0,N,V,M);  errors += fabs(M - (  9.3724));  printf("M: %g  =>  %g\n",M,  9.3724);
                    e1->CalcRes(1.0,N,V,M);  errors += fabs(M - (  0.0000));  printf("M: %g  =>  %g\n",M,  0.0000);
                    e2->CalcRes(0.0,N,V,M);  errors += fabs(M - (-10.8887));  printf("M: %g  =>  %g\n",M,-10.8887);
                    e2->CalcRes(1.0,N,V,M);  errors += fabs(M - (  0.0000));  printf("M: %g  =>  %g\n",M,  0.0000);
                    e3->CalcRes(0.0,N,V,M);  errors += fabs(M - (-10.8887));  printf("M: %g  =>  %g\n",M,-10.8887);
                    e3->CalcRes(1.0,N,V,M);  errors += fabs(M - ( 13.9279));  printf("M: %g  =>  %g\n",M, 13.9279);
                    e4->CalcRes(0.0,N,V,M);  errors += fabs(M - ( 13.9279));  printf("M: %g  =>  %g\n",M, 13.9279);
                    e4->CalcRes(1.0,N,V,M);  errors += fabs(M - (-18.7555));  printf("M: %g  =>  %g\n",M,-18.7555);
                    tol = 1.0e-3;

                    double ux = dom.Nods[4]->U("ux")*1000.0;
                    double uy = dom.Nods[4]->U("uy")*1000.0;
                    double wz = dom.Nods[4]->U("wz");
                    printf("ux, uy, uz: %g => %g,  %g => %g,  %g => %g  [mm, mm,  ]\n", ux,6.476e-3, uy,-6.212e-3, wz,9.394e-5);
                    err_u += fabs(ux-( 6.476e-3));
                    err_u += fabs(uy-(-6.212e-3));
                    err_u += fabs(wz-( 9.394e-5));
                    ux = dom.Nods[3]->U("ux")*1000.0;
                    uy = dom.Nods[3]->U("uy")*1000.0;
                    wz = dom.Nods[3]->U("wz");
                    printf("ux, uy, uz: %g => %g,  %g => %g,  %g => %g  [mm, mm,  ]\n", ux,1.778e-2, uy,-3.243e-2, wz,-1.601e-4);
                    err_u += fabs(ux-( 1.778e-2));
                    err_u += fabs(uy-(-3.243e-2));
                    err_u += fabs(wz-(-1.601e-4));
                    tol_u = 1.0e-3;
                }

                break;
            }
            case 6:
            {
                break;
            }
            case 7:
            {
                FEM::Beam const * e0 = static_cast<FEM::Beam const*>(dom.Eles[0]);
                FEM::Beam const * e1 = static_cast<FEM::Beam const*>(dom.Eles[1]);
                e0->CalcRes(0.,  N,V,M);   errors += fabs(M-(-2.5));
                e0->CalcRes(1.,  N,V,M);   errors += fabs(M-0.0);
                e1->CalcRes(0.5, N,V,M);   errors += fabs(M-(10.0*0.5*0.5/8.));
                err_u += fabs(dom.Nods[1]->U("uy") - (-0.182292));
                err_u += fabs(dom.Nods[1]->U("wz") - (-25./48.));
                err_u += fabs(dom.Nods[3]->U("uy") - (-0.182292));
                err_u += fabs(dom.Nods[3]->U("wz") - ( 0.31250));
                tol_u = 1.0e-6;
                break;
            }
        }

        printf("Errors (all)                = %s%.8e%s\n", (errors>tol  ?TERM_RED:TERM_GREEN), errors, TERM_RST);
        printf("Errors (displacements, all) = %s%.8e%s\n", (err_u >tol_u?TERM_RED:TERM_GREEN), err_u,  TERM_RST);
        if (errors>tol)  with_error = true;
        if (err_u>tol_u) with_error = true;
    }

    // end
    return with_error;
}
MECHSYS_CATCH
