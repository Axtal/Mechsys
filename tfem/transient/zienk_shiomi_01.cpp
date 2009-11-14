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

// MechSys
#include "mesh/structured.h"
#include "fem/elems/quad8.h"
#include "fem/hydromechelem.h"
#include "models/linelastic.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "util/maps.h"
#include "util/fatal.h"

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;
using Util::_4;
using Util::_6_3;
using Util::_8s;
using Util::PI;

double Multiplier (double t)
{
    double tsw = 0.1;
    if (t<tsw) return sin(PI*t/(2.0*tsw));
    else       return 1.0;
}

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
    
    Array<Mesh::Block> blks(1);
    blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                 -100., 0.0,  0.0,
                    0., 1.0,  0.0,
                    0., 1.0, 30.0,
                 -200., 0.0, 30.0,  -10., -20., -30., -20.);
    blks[0].SetNx (1);
    blks[0].SetNy (10);
    Mesh::Structured mesh(/*NDim*/2);
    mesh.Generate (blks,/*O2*/true);
    mesh.WriteMPY ("zienk_shiomi_01");

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    double Kf   = 100.0;
    double Ks   =  10.0;
    double alp  =   1.0;
    double n    =   0.3;
    double rhoS =   2.0;
    double rhoF =   1.0;
    double k    =   1.0e-2;
    double Q    = 1.0/(n/Kf+(alp-n)/Ks);
    Dict prps;
    prps.Set(-1, "prob geom psa  alp n rhoS rhoF k Q", PROB("HMEquilib"), GEOM("Quad8"), TRUE, alp, n, rhoS, rhoF, k, Q);

    // models
    double E  = 30.0;
    double nu =  0.2;
    Dict mdls;
    mdls.Set(-1, "name psa  E nu", MODEL("LinElastic"), TRUE, E, nu);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);
    dom.SetOutNods ("zienk_shiomi_01", Array<int>(-100, -200));
    dom.SetOutEles ("zienk_shiomi_01", Array<int>(0, /*JustOne*/true));
    dom.MFuncs[-30] = &Multiplier;

    // solver
    FEM::Solver sol(dom);
    sol.DampTy = FEM::Solver::HMCoup_t;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set(-10, "ux uy Ux Uy",  0.0, 0.0, 0.0, 0.0);
    bcs.Set(-20, "ux Ux",        0.0, 0.0);
    bcs.Set(-30, "qn mfunc",    -1.0, 0.0);
    bcs.Set(-30, "pw",           0.0);
    dom.SetBCs (bcs);
    sol.DynSolve (2.0, 0.01, 0.05, "zienk_shiomi_01");

    return 0;
}
MECHSYS_CATCH
