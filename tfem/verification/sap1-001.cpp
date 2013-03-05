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
#include <mechsys/fem/fem.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;
using FEM::Domain;
using FEM::Solver;

double lcase0_M_on_elem1 (double r)
{
    double q = 0.15;
    double x = (1.0-r)*10.0;
    return 1.2*q*x - 0.4*q*x*x;
}

double lcase3_M_on_elem0 (double r)
{
    double x = r*10.0;
    return 12.0*x-x*x;
}

double lcase3_M_on_elem1 (double r)
{
    double x = (1.0-r)*10.0;
    return 12.0*x-x*x;
}

double lcase4_M_on_elem0 (double r)
{
    double x = r*10.0;
    return 5.0752*x - 0.9984*x*x/3.0 - (16.0-x)*0.0624*x*x/6.0;
}

double lcase4_M_on_elem1 (double r)
{
    double x = (1.0-r)*10.0;
    return 1.7472*x - 0.6*0.0624*x*x*x/6.0;
}

int main(int argc, char **argv) try
{
    // mesh
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize  (4, 3);
    mesh.SetVert  (0, -100,  4.0,  0.0);
    mesh.SetVert  (1, -400,  4.0, 10.0);
    mesh.SetVert  (2, -200, 12.0, 16.0);
    mesh.SetVert  (3, -300,  0.0, 10.0);
    mesh.SetCell  (0,   -1, Array<int>(0,1));
    mesh.SetCell  (1,   -2, Array<int>(1,2));
    mesh.SetCell  (2,   -3, Array<int>(3,1));
    mesh.WriteMPY ("sap1-001_mesh");

    // elements properties
    double gra = 32.18503937; // ft/s2
    double gam = 0.15;        // kip/ft3
    double rho = gam/gra;     //
    double E   = 518400.0;    // kip/ft2
    double A   = 1.0;         // ft2
    double Izz = 1.0/12.0;    // ft4
    Dict prps;
    prps.Set (-1, "prob fra  rho E A Izz", PROB("Beam"), 1.,  rho, E, A, Izz);
    prps.Set (-2, "prob fra  rho E A Izz", PROB("Beam"), 1.,  rho, E, A, Izz);
    prps.Set (-3, "prob fra  rho E A Izz", PROB("Beam"), 1.,  rho, E, A, Izz);

    // load cases
    Array<Dict> BC(7);
    for (size_t i=0; i<7; ++i) // default restraints
    {
        BC[i].Set (-100, "ux uy", 0.0, 0.0);
        BC[i].Set (-200, "ux",    0.0);
    }

    // 0) self weight
    BC[0].Set (-1, "gravity", gra);
    BC[0].Set (-2, "gravity", gra);
    BC[0].Set (-3, "gravity", gra);

    // 1) global uniform distributed load on frame element 2, plus concentrated load on node 3
    BC[1].Set (-300, "fy", -10.0);
    BC[1].Set (-3,   "qn",  -1.8);

    // 2) global joint force and moment at node 1
    BC[2].Set (-400, "fy", -17.2);
    BC[2].Set (-400, "mz",  54.4);

    // 3) uniformly distributed load on elements 0 and 1
    BC[3].Set (-1, "qn", -2.0);
    BC[3].Set (-2, "qn", -2.0);

    // 4) trapezoidal load on elements 0 and 1
    BC[4].Set (-1, "qnl_qnr_-0.9984_-0.3744", 1.);
    BC[4].Set (-2, "qnl_qnr_-0.3744_0.0", 1.);

    // solve
    double errors = 0.0;
    double err_u  = 0.0;
    double tol    = 1.0e-10;
    double tol_u  = 1.0e-4;
    double V,N,M;
    for (size_t lcase=0; lcase<5; ++lcase)
    {
        // domain
        FEM::Domain::WithInfo = false;
        FEM::Domain dom(mesh, prps, Dict(), Dict());

        // solve
        SDPair flags;
        flags.Set ("fe", 1.);
        FEM::STDSolver sol(dom, flags);
        dom.SetBCs (BC[lcase]);
        sol.Solve  (1);
        printf("\n%s===================================== Load case # %zd =====================================%s\n",TERM_YELLOW_BLUE,lcase,TERM_RST);

        // output
        String buf;
        buf.Printf("sap1-001_lcase%d",lcase);
        //dom.PrintResults ("%10.5f", /*with_elems*/false);
        FEM::MPyPrms mpy_prms;
        dom.WriteMPY     (buf.CStr(), mpy_prms);
        dom.WriteVTU     (buf.CStr());
        printf("\nNorm(R)                = %s%8e%s\n",(sol.ResidOK()?TERM_GREEN:TERM_RED),sol.NormR,TERM_RST);

        // elements
        FEM::Beam const * e0 = static_cast<FEM::Beam const *>(dom.Eles[0]);
        FEM::Beam const * e1 = static_cast<FEM::Beam const *>(dom.Eles[1]);
        FEM::Beam const * e2 = static_cast<FEM::Beam const *>(dom.Eles[2]);

        // check
        if (lcase==0)
        {
            double q = gam*A;
            e0->CalcRes(1.0,V,N,M);   errors += fabs(M - (-20.0*q)); // M at right of beam #0 == -20q
            for (size_t i=0; i<11; ++i)
            {
                double r = i/10.0;
                e1->CalcRes (r, V,N,M);
                //printf("M along E1, at r=%.1f: Num = %g,  An = %g\n",r,M,lcase0_M_on_elem1(r));
                errors += fabs(M - lcase0_M_on_elem1(r));
            }
            e2->CalcRes(1.0,V,N,M);   errors += fabs(M - (-q*4.0*4.0/2.0)); // M at right of beam #2 == -qL^2/2
            Table disp;
            disp.Set("         ux            uy            wz", 4,
                       0.0000E+00,   0.0000E+00,   2.8001E-04,
                      -1.9712E-02,  -6.5972E-04,  -6.7216E-05,
                       0.0000E+00,  -2.7208E-02,  -3.2185E-04,
                      -1.9712E-02,   1.2333E-03,  -3.0179E-05);
            for (size_t i=0; i<4; ++i)
            {
                err_u += fabs(dom.Nods[i]->U("ux")*12.0 - disp("ux",i));
                err_u += fabs(dom.Nods[i]->U("uy")*12.0 - disp("uy",i));
                err_u += fabs(dom.Nods[i]->U("wz")      - disp("wz",i));
                //printf("Displacements: %g/%g,  %g/%g,  %g/%g\n", dom.Nods[i]->U[dom.Nods[i]->UMap("ux")]*12.0, disp("ux",i), dom.Nods[i]->U[dom.Nods[i]->UMap("uy")]*12.0, disp("uy",i), dom.Nods[i]->U[dom.Nods[i]->UMap("wz")], disp("wz",i));
            }
            printf("Errors                 = %s%.8e%s\n", (errors>tol  ?TERM_RED:TERM_GREEN), errors, TERM_RST);
            printf("Errors (displacements) = %s%.8e%s\n", (err_u >tol_u?TERM_RED:TERM_GREEN), err_u,  TERM_RST);
        }
        else if (lcase==1)
        {
            e0->CalcRes(1.0,V,N,M);  errors += fabs(M - 34.0);
            e1->CalcRes(0.0,V,N,M);  errors += fabs(M + 20.4);
            e2->CalcRes(1.0,V,N,M);  errors += fabs(M + (10.0*4.0+1.8*4.0*4.0/2.0));
            printf("Errors                 = %s%.8e%s\n", (errors>tol  ?TERM_RED:TERM_GREEN), errors, TERM_RST);
        }
        else if (lcase==2)
        {
            e0->CalcRes(1.0,V,N,M);  errors += fabs(M - 34.0);
            e1->CalcRes(0.0,V,N,M);  errors += fabs(M + 20.4);
            printf("Errors                 = %s%.8e%s\n", (errors>tol  ?TERM_RED:TERM_GREEN), errors, TERM_RST);
        }
        else if (lcase==3)
        {
            for (size_t i=0; i<11; ++i)
            {
                double r = i/10.0;
                e0->CalcRes (r, V,N,M);
                //printf("M along E0, at r=%.1f: Num = %g,  An = %g\n",r,M,lcase3_M_on_elem0(r));
                errors += fabs(M - lcase3_M_on_elem0(r));
            }
            for (size_t i=0; i<11; ++i)
            {
                double r = i/10.0;
                e1->CalcRes (r, V,N,M);
                //printf("M along E1, at r=%.1f: Num = %g,  An = %g\n",r,M,lcase3_M_on_elem1(r));
                errors += fabs(M - lcase3_M_on_elem1(r));
            }
            printf("Errors                 = %s%.8e%s\n", (errors>tol  ?TERM_RED:TERM_GREEN), errors, TERM_RST);
        }
        else if (lcase==4)
        {
            for (size_t i=0; i<11; ++i)
            {
                double r = i/10.0;
                e0->CalcRes (r, V,N,M);
                //printf("M along E0, at r=%.1f: Num = %g,  An = %g\n",r,M,lcase4_M_on_elem0(r));
                errors += fabs(M - lcase4_M_on_elem0(r));
            }
            for (size_t i=0; i<11; ++i)
            {
                double r = i/10.0;
                e1->CalcRes (r, V,N,M);
                //printf("M along E1, at r=%.1f: Num = %g,  An = %g\n",r,M,lcase4_M_on_elem1(r));
                errors += fabs(M - lcase4_M_on_elem1(r));
            }
            Table disp;
            disp.Set("         ux            uy            wz", 4,
                       0.0000E+00,   0.0000E+00,  -2.0806E-03,
                       1.1583E-01,  -3.4667E-04,   5.4343E-04,
                       0.0000E+00,   1.5319E-01,   2.2045E-03,
                       1.1583E-01,  -2.6431E-02,   5.4343E-04);
            for (size_t i=0; i<4; ++i)
            {
                err_u += fabs(dom.Nods[i]->U("ux")*12.0 - disp("ux",i));
                err_u += fabs(dom.Nods[i]->U("uy")*12.0 - disp("uy",i));
                err_u += fabs(dom.Nods[i]->U("wz")      - disp("wz",i));
                //printf("Displacements: %g/%g,  %g/%g,  %g/%g\n", dom.Nods[i]->U[dom.Nods[i]->UMap("ux")]*12.0, disp("ux",i), dom.Nods[i]->U[dom.Nods[i]->UMap("uy")]*12.0, disp("uy",i), dom.Nods[i]->U[dom.Nods[i]->UMap("wz")], disp("wz",i));
            }
            printf("Errors                 = %s%.8e%s\n", (errors>tol  ?TERM_RED:TERM_GREEN), errors, TERM_RST);
            printf("Errors (displacements) = %s%.8e%s\n", (err_u >tol_u?TERM_RED:TERM_GREEN), err_u,  TERM_RST);
        }
    }

    // end
    printf("\n%s========================================= End ===========================================%s\n\n",TERM_YELLOW_BLUE,TERM_RST);
    printf("Errors (all)                = %s%.8e%s\n", (errors>tol  ?TERM_RED:TERM_GREEN), errors, TERM_RST);
    printf("Errors (displacements, all) = %s%.8e%s\n", (err_u >tol_u?TERM_RED:TERM_GREEN), err_u,  TERM_RST);
    return (errors>tol);
    return (err_u >tol_u);
}
MECHSYS_CATCH
