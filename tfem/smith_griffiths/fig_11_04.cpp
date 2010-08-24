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

/*  Smith & Griffiths (2004): Figure 11.4 p475 
 *  ==========================================  */

// STL
#include <iostream>
#include <cmath>    // for cos

// MechSys
#include <mechsys/mesh/mesh.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/fem/elems/quad8.h>
#include <mechsys/fem/equilibelem.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solver.h>
#include <mechsys/models/linelastic.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/util.h>
#include <mechsys/util/fatal.h>
#include <mechsys/draw.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

void OutFun (FEM::Solver const & Sol, void * Dat)
{
    size_t idx_nod = 17;
    int    eq_ux   = Sol.Dom.Nods[idx_nod]->Eq("ux");
    int    eq_uy   = Sol.Dom.Nods[idx_nod]->Eq("uy");
    double ux      = Sol.U(eq_ux);
    double uy      = Sol.U(eq_uy);
    double vx      = Sol.V(eq_ux);
    double vy      = Sol.V(eq_uy);
    double ax      = Sol.A(eq_ux);
    double ay      = Sol.A(eq_uy);
    if (Sol.IdxOut==0)
    {
        std::ofstream of("fig_11_04.out", std::ios::out);
        String str;
        str.Printf ("%15s  %15s %15s  %15s %15s  %15s %15s\n", "Time","ux","uy","vx","vy","ax","ay");
        of << str;
        str.Printf ("%15.8e  %15.8e %15.8e  %15.8e %15.8e  %15.8e %15.8e\n", Sol.Time,ux,uy,vx,vy,ax,ay);
        of << str;
        //of << "Time    ux   uy     vx    vy     ax    ay\n";
        //of << Sol.Time << " " << ux << " " << uy << " " << vx << " " << vy << " " << ax << " " << ay << endl;
        of.close();
    }
    else
    {
        std::ofstream of("fig_11_04.out", std::ios::app);
        String str;
        str.Printf ("%15.8e  %15.8e %15.8e  %15.8e %15.8e  %15.8e %15.8e\n", Sol.Time,ux,uy,vx,vy,ax,ay);
        of << str;
        //of << Sol.Time << " " << ux << " " << uy << " " << vx << " " << vy << " " << ax << " " << ay << endl;
        of.close();
    }
}

double Multiplier (double t)
{
    double omega = 0.3;
    return cos(omega*t);
}

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
    
    double H = 1.0;  // height
    double h = 0.5;  // half height
    double L = 4.0;  // length
    double l = L/3.; // element length
    double m = l/2.; // half length (for mid nodes)
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize    (18/*verts*/, 3/*cells*/);
    mesh.SetVert    ( 0,    0,  0.0,    H, 0);
    mesh.SetVert    ( 1,    0,  0.0,    h, 0);
    mesh.SetVert    ( 2,    0,  0.0,  0.0, 0);
    mesh.SetVert    ( 3,    0,    m,    H, 0);
    mesh.SetVert    ( 4,    0,    m,  0.0, 0);
    mesh.SetVert    ( 5,    0,    l,    H, 0);
    mesh.SetVert    ( 6,    0,    l,    h, 0);
    mesh.SetVert    ( 7,    0,    l,  0.0, 0);
    mesh.SetVert    ( 8,    0,  l+m,    H, 0);
    mesh.SetVert    ( 9,    0,  l+m,  0.0, 0);
    mesh.SetVert    (10,    0,  l+l,    H, 0);
    mesh.SetVert    (11,    0,  l+l,    h, 0);
    mesh.SetVert    (12,    0,  l+l,  0.0, 0);
    mesh.SetVert    (13,    0,  L-m,    H, 0);
    mesh.SetVert    (14,    0,  L-m,  0.0, 0);
    mesh.SetVert    (15,    0,    L,    H, 0);
    mesh.SetVert    (16,    0,    L,    h, 0);
    mesh.SetVert    (17, -100,    L,  0.0, 0);
    mesh.SetCell    ( 0,   -1, Array<int>( 2, 7, 5, 0,  4, 6, 3, 1));
    mesh.SetCell    ( 1,   -1, Array<int>( 7,12,10, 5,  9,11, 8, 6));
    mesh.SetCell    ( 2,   -1, Array<int>(12,17,15,10, 14,16,13,11));
    mesh.SetBryTag  ( 0, 3, -10);
    mesh.WriteMPY   ("fig_11_04");
    
    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom active psa rho", PROB("Equilib"), GEOM("Quad8"), 1.0, 1.0, 1.0);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu psa", MODEL("LinElastic"), 1.0, 0.3, 1.0);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);
    dom.MFuncs[-100] = &Multiplier; // set database of callbacks

    // solver
    FEM::Solver sol(dom, &OutFun);
    //sol.DScheme = FEM::Solver::SS22_t;
    sol.DScheme = FEM::Solver::GN22_t;
    sol.DampAm  = 0.005;
    sol.DampAk  = 0.272;
    sol.DampTy  = FEM::Solver::Rayleigh_t;
    //sol.CteTg   = true;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set( -10, "ux uy",    0.0, 0.0);
    bcs.Set(-100, "fy mfunc", 1.0, 0.0);
    dom.SetBCs (bcs);
    //cout << dom << endl;
    sol.DynSolve (/*tf*/100, /*dt*/1.0, /*dtOut*/1.0, "fig_11_04");

    return 0.0;
}
MECHSYS_CATCH
