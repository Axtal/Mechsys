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

/*  Zienkiewicz & Shiomi (1984) Example 1 (Fig 1)  p85             *
 *  Int J Num Methods and Analytical in Geomechanics Vol 8, 71-96  *
 *  ============================================================== */

// STL
#include <iostream>
#include <fstream>

// MechSys
#include <mechsys/mesh/structured.h>
#include <mechsys/fem/elems/quad4.h>
#include <mechsys/fem/elems/quad8.h>
#include <mechsys/fem/hydromechelem.h>
#include <mechsys/models/linelastic.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solver.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;
using Util::_4;
using Util::_6_3;
using Util::_8s;
using Util::PI;

void DbgFun (FEM::Solver const & Sol, void * Dat)
{
    String fn;
    fn.Printf ("zienk_shiomi_01_%06.3f.res", Sol.Time);
    std::ofstream of(fn.CStr(),std::ios::out);
    of<< _8s << "x" << _8s << "y" << _8s << "pw" << "\n";
    for (size_t i=0; i<Sol.Dom.Nods.Size(); ++i)
    {
        double x = Sol.Dom.Nods[i]->Vert.C(0);
        double y = Sol.Dom.Nods[i]->Vert.C(1);
        if (fabs(x)<0.001)
        {
            size_t ipw = Sol.Dom.Nods[i]->UMap("pw");
            long   eq  = Sol.Dom.Nods[i]->EQ[ipw];
            double pw  = Sol.U[eq];
            of << _8s << x << _8s << y << _8s << pw << "\n";
        }
    }
    of.close();
}

double Multiplier (double t)
{
    double tsw = 0.1;
    if (t<tsw) return sin(PI*t/(2.0*tsw));
    else       return 1.0;
}

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
    
    bool o2 = true;
    Array<Mesh::Block> blks(1);
    blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                 -100., 0.0,  0.0,
                    0., 1.0,  0.0,
                    0., 1.0, 30.0,
                 -200., 0.0, 30.0,  -10., -20., -30., -20.);
    blks[0].SetNx (1);
    blks[0].SetNy (10);
    Mesh::Structured mesh(/*NDim*/2);
    mesh.Generate (blks,o2);
    mesh.WriteMPY ("zienk_shiomi_01");

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // constants
    double qn   =  -1.0;    // kPa (loading)
    double E    =  30.0e+3; // kPa
    double nu   =   0.2;
    double Kf   = 100.0e+3; // kPa
    double Ks   =  10.0e+3; // kPa
    double alp  =   1.0;
    double n    =   0.3;
    double rhoS =   2.0;    // Mg/m^3
    double rhoF =   1.0;    // Mg/m^3
    double k    =   1.0e-2/9.81; // m/s
    double Q    = 1.0/(n/Kf+(alp-n)/Ks);
    cout << "Q = " << Q << endl;

    // elements properties
    Dict prps;
    if (o2) prps.Set(-1, "prob geom psa  alp n rhoS rhoF k Q", PROB("HMEquilib"), GEOM("Quad8"), TRUE, alp, n, rhoS, rhoF, k, Q);
    else    prps.Set(-1, "prob geom psa  alp n rhoS rhoF k Q", PROB("HMEquilib"), GEOM("Quad4"), TRUE, alp, n, rhoS, rhoF, k, Q);

    // models
    Dict mdls;
    mdls.Set(-1, "name psa  E nu", MODEL("LinElastic"), TRUE, E, nu);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);
    if (o2) dom.SetOutNods ("zienk_shiomi_01", Array<int>(-100, -200));
    else    dom.SetOutNods ("zienk_shiomi_01", Array<int>(-100, -200, 2,3, 10,11, 18,19));
            dom.SetOutEles ("zienk_shiomi_01", Array<int>(0, 5, 9));
    dom.MFuncs[-30] = &Multiplier;

    // solver
    FEM::Solver sol(dom, &DbgFun);
    sol.DampTy = FEM::Solver::HMCoup_t;
    sol.DScheme = FEM::Solver::GNHMCoup_t;
    //sol.DynTh1 = 0.6;
    //sol.DynTh2 = 0.605;

    // initial state -------------------------------------------------------
    //double gamw = 10.0; // kN/m^3
    //for (size_t i=0; i<dom.Nods.Size(); ++i)
    //{
        //size_t ipw = dom.Nods[i]->UMap("pw");
        //double hp  = 30.0 - dom.Nods[i]->Vert.C(1);
        //dom.Nods[i]->U[ipw] = gamw*hp;
    //}

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    //bcs.Set(-10, "ux uy flux",  0.0, 0.0, 0.0);
    //bcs.Set(-20, "ux flux",     0.0, 0.0);
    bcs.Set(-10, "ux uy Ux Uy flux",  0.0, 0.0, 0.0, 0.0, 0.0);
    bcs.Set(-20, "ux Ux flux",        0.0, 0.0, 0.0);
    bcs.Set(-30, "qn mfunc pw",        qn, 0.0, 0.0);
    dom.SetBCs (bcs);
    cout << dom << endl;
    sol.DynSolve (2.0, 0.01, 0.05);//, "zienk_shiomi_01");

    return 0;
}
MECHSYS_CATCH
