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
#include "mesh/structured.h"
#include "fem/elems/quad8.h"
#include "fem/equilibelem.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "models/linelastic.h"
#include "util/maps.h"
#include "util/util.h"
#include "util/fatal.h"
#include "draw.h"

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

void CentroidSolution (double Pressure, double radius, double Radius, FEM::Element const & E, Table & Sol, bool Silent=true)
{
    // results
    Vec_t  ct;
    SDPair res;
    E.Centroid (ct);
    E.GetState (res);
    double x = ct(0);
    double y = ct(1);

    // solution in polar coordinates
    double p      = Pressure;
    double ri     = radius;
    double ro     = Radius;
    double r      = sqrt(x*x+y*y);
    double coef   = p*ri*ri/(ro*ro-ri*ri);
    double sig_r  = coef*(1.0-ro*ro/(r*r));
    double sig_t  = coef*(1.0+ro*ro/(r*r));
    double sig_rt = 0.0;
    
    // rotation to x-y coordinates
    double c      = x/r;
    double s      = y/r;
    double cc     = c*c;
    double ss     = s*s;
    double cs     = c*s;
    double sig_x  = cc*sig_r + ss*sig_t -  2.0*cs*sig_rt;
    double sig_y  = ss*sig_r + cc*sig_t +  2.0*cs*sig_rt;
    double sig_xy = cs*sig_r - cs*sig_t + (cc-ss)*sig_rt;

    if (!Silent)
    {
        // rotation to r-t coordinates
        double sx  = res("sx");
        double sy  = res("sy");
        double sxy = res("sxy");
        double sr  =  cc*sx + ss*sy +  2.0*cs*sxy;
        double st  =  ss*sx + cc*sy -  2.0*cs*sxy;
        double srt = -cs*sx + cs*sy + (cc-ss)*sxy;

        // output
        cout << Util::_4 << E.Cell.ID;
        cout << Util::_6_3  << x      << Util::_6_3  << y   << Util::_6_3 << r << "  ";
        cout << Util::_10_6 << sig_r  << Util::_10_6 << sr  << "  ";
        cout << Util::_10_6 << sig_t  << Util::_10_6 << st  << "  ";
        cout << Util::_10_6 << sig_rt << Util::_10_6 << srt << endl;
    }

    // error
    Sol.SetVal ("sx",  E.Cell.ID, sig_x);
    Sol.SetVal ("sy",  E.Cell.ID, sig_y);
    Sol.SetVal ("sxy", E.Cell.ID, sig_xy);
}

void ElemSolution (double Pressure, double radius, double Radius, FEM::Element const & E, Table & Sol, bool Silent=true)
{
    Array<SDPair> res(E.GE->NIP);
    E.GetState (res);

    for (size_t i=0; i<E.GE->NIP; ++i)
    {
        // coordinates of IP
        Vec_t X;
        E.CoordsOfIP (i, X);
        double x = X(0);
        double y = X(1);

        // solution in polar coordinates
        double p      = Pressure;
        double ri     = radius;
        double ro     = Radius;
        double r      = sqrt(x*x+y*y);
        double coef   = p*ri*ri/(ro*ro-ri*ri);
        double sig_r  = coef*(1.0-ro*ro/(r*r));
        double sig_t  = coef*(1.0+ro*ro/(r*r));
        double sig_rt = 0.0;
        
        // rotation to x-y coordinates
        double c      = x/r;
        double s      = y/r;
        double cc     = c*c;
        double ss     = s*s;
        double cs     = c*s;
        double sig_x  = cc*sig_r + ss*sig_t -  2.0*cs*sig_rt;
        double sig_y  = ss*sig_r + cc*sig_t +  2.0*cs*sig_rt;
        double sig_xy = cs*sig_r - cs*sig_t + (cc-ss)*sig_rt;

        if (!Silent)
        {
            // rotation to r-t coordinates
            double sx  = res[i]("sx");
            double sy  = res[i]("sy");
            double sxy = res[i]("sxy");
            double sr  =  cc*sx + ss*sy +  2.0*cs*sxy;
            double st  =  ss*sx + cc*sy -  2.0*cs*sxy;
            double srt = -cs*sx + cs*sy + (cc-ss)*sxy;

            // output
            cout << Util::_4 << E.Cell.ID;
            cout << Util::_6_3  << x      << Util::_6_3  << y   << Util::_6_3 << r << "  ";
            cout << Util::_10_6 << sig_r  << Util::_10_6 << sr  << "  ";
            cout << Util::_10_6 << sig_t  << Util::_10_6 << st  << "  ";
            cout << Util::_10_6 << sig_rt << Util::_10_6 << srt << endl;
        }

        // error
        Sol.SetVal ("sx",  i + E.Cell.ID*E.GE->NIP, sig_x);
        Sol.SetVal ("sy",  i + E.Cell.ID*E.GE->NIP, sig_y);
        Sol.SetVal ("sxy", i + E.Cell.ID*E.GE->NIP, sig_xy);
    }
}

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    double r =  5.0; // inner radius
    double R = 15.0; // outer radius
    double a = r*cos(Util::PI/4.0);
    double A = R*cos(Util::PI/4.0);
    double b = r*cos(Util::PI/8.0);
    double B = R*cos(Util::PI/8.0);
    double c = r*sin(Util::PI/8.0);
    double C = R*sin(Util::PI/8.0);
    double m = (R+r)/2.0;
    double M = (R+r)*sin(Util::PI/4.0)/2.0;

    Array<Mesh::Block> blks(2);
    blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/8,
                  0.0,    r, 0.0, 
                  0.0,    R, 0.0, 
                  0.0,    A,   A, 
                  0.0,    a,   a, 
                  0.0,    m, 0.0,
                  0.0,    B,   C, 
                  0.0,    M,   M,
                  0.0,    b,   c, 
                 -10.0,0.0,0.0,-30.0);
    blks[1].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/8,
                  0.0,    a,   a, 
                  0.0,    A,   A, 
                  0.0,  0.0,   R, 
                  0.0,  0.0,   r, 
                  0.0,    M,   M,
                  0.0,    C,   B, 
                  0.0,  0.0,   m,
                  0.0,    c,   b, 
                  0.0,0.0,-20.0,-30.0);
    blks[0].SetNx (7, /*Ax*/1.0, /*NonLin*/false);
    blks[0].SetNy (5);
    blks[1].SetNx (7, /*Ax*/1.0, /*NonLin*/false);
    blks[1].SetNy (5);
    //blks[0].SetNx (1, /*Ax*/1.0, /*NonLin*/false);
    //blks[0].SetNy (1);
    //blks[1].SetNx (1, /*Ax*/1.0, /*NonLin*/false);
    //blks[1].SetNy (1);
    Mesh::Structured mesh(/*NDim*/2);
    mesh.Generate (blks,/*O2*/true);
    //mesh.WriteVTU ("ex79");
    mesh.WriteMPY  ("ex79",/*OnlyMesh*/true);

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom active  gam psa  nip", PROB("Equilib"), GEOM("Quad8"), TRUE,  0.0, TRUE, 9.0);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu psa", MODEL("LinElastic"), 10000.0, 0.25, TRUE);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);

    // solver
    FEM::Solver sol(dom);

    // stage # 1 -----------------------------------------------------------
    double p = 20; // pressure
    Dict bcs;
    bcs.Set(-10, "uy", 0.0)
       .Set(-20, "ux", 0.0)
       .Set(-30, "qn", -p);
    dom.SetBCs (bcs);
    sol.Solve  ();

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    //cout << dom;
    //dom.PrintResults (cout, Util::_12_6);

    // draw elements with IPs
    std::ofstream of("ex79_elements.draw", std::ios::out);
    MPL::Header   (of);
    for (size_t i=0; i<dom.Eles.Size(); ++i) dom.Eles[i]->Draw (of);
    MPL::AddPatch (of);
    MPL::Show     (of);
    of.close      ();

    //////////////////////////////////////////////////////////////////////////////////////// Check /////

    bool silent = true;
    if (!silent)
    {
        cout << "\n[1;37m--- Solution -----------------------------------------------------------------\n";
        cout << "\n" << Util::_4 << "Elem";
        cout << Util::_6_3  << "x"      << Util::_6_3  << "y"   << Util::_6_3 << "r" << "  ";
        cout << Util::_10_6 << "(c)sr"  << Util::_10_6 << "sr"  << "  ";
        cout << Util::_10_6 << "(c)st"  << Util::_10_6 << "st"  << "  ";
        cout << Util::_10_6 << "(c)srt" << Util::_10_6 << "srt" << "[0m\n";
    }

    // solution
    Table nod_sol, ele_sol, elem_sol;
    ele_sol .SetZero("sx sy sxy", dom.Eles.Size());
    elem_sol.SetZero("sx sy sxy", dom.Eles.Size()*dom.Eles[0]->GE->NIP);
    for (size_t i=0; i<dom.Eles.Size(); ++i)
    {
        CentroidSolution (p, r, R, (*dom.Eles[i]), ele_sol,  silent);
        ElemSolution     (p, r, R, (*dom.Eles[i]), elem_sol, silent);
    }

    // error tolerance
    SDPair nod_tol, ele_tol;
    ele_tol.Set("sx sy sxy", 2.0e-1,2.0e-1,1.0e-1);

    // return error flag
           dom.CheckError (cout, elem_sol, ele_tol);
    return dom.CheckError (cout, nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH
