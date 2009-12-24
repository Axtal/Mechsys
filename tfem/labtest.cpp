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
#include <mechsys/mesh/structured.h>
#include <mechsys/fem/elems/hex8.h>
#include <mechsys/fem/elems/hex20.h>
#include <mechsys/fem/equilibelem.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solver.h>
#include <mechsys/models/linelastic.h>
#include <mechsys/models/elastoplastic.h>
#include <mechsys/models/camclay.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;
using Util::_6_3;
using Util::_8s;

struct DbgDat
{
    size_t node_id;
    size_t dof_idx;
    std::ofstream of;
    ~DbgDat () { of.close (); }
     DbgDat () : node_id(0), dof_idx(0)
    {
         of.open ("labtest_uf.res",std::ios::out);
         of<<_6_3<<"Time"<<_8s<<"u"<< _8s<<"f_int"<<_8s<<"f_ext\n";
    }
};

void DbgFun (FEM::Solver const & Sol, void * Dat)
{
    DbgDat * dat = static_cast<DbgDat*>(Dat);
    long eq = Sol.Dom.Nods[dat->node_id]->EQ[dat->dof_idx];
    dat->of << _6_3 << Sol.Time << _8s << Sol.U(eq) << _8s << Sol.F_int(eq) << _8s << Sol.F(eq) << endl;
}

int main(int argc, char **argv) try
{
    // index of model
    size_t idx_mdl = 0;
    if (argc>1) idx_mdl = atoi(argv[1]);

    // path
    double qx =  0.0; // delta pressure on the x-side (keep it constant)
    double qy =  0.0; // delta pressure on the y-side (keep it constant)
    double uz = -0.07;

    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    bool   o2 = false;
    size_t nd = 1;
    Mesh::Structured mesh(/*NDim*/3);
    mesh.GenBox (/*O2*/o2,/*Nx*/nd,/*Ny*/nd,/*Nz*/nd);

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    if (o2) prps.Set(-1, "prob geom", PROB("Equilib"), GEOM("Hex20"));
    else    prps.Set(-1, "prob geom", PROB("Equilib"), GEOM("Hex8"));

    // models
    Dict mdls;
    switch (idx_mdl)
    {
        case 0: // Cam clay
        {
            cout << "\n[1;33m====================================== Cam clay ====================================[0m\n";
            mdls.Set(-1, "name  lam kap nu phi", MODEL("CamClay"), 0.01, 0.001, 0.3, M2Phi(1.0,"cam"));
            break;
        }
        case 1: // elasto-perfect-plastic with von Mises FC
        {
            cout << "\n[1;33m======================= Elasto-perfect-plastic with von Mises ======================[0m\n";
            mdls.Set(-1, "name E nu fc sY", MODEL("ElastoPlastic"), 4000.0, 0.3, FAILCRIT("VM"), 150.0);
            break;
        }
        case 2: // linear elastic
        {
            cout << "\n[1;33m==================================== Linear elastic ================================[0m\n";
            mdls.Set(-1, "name E nu", MODEL("LinElastic"), 4000.0, 0.3);
            break;
        }
        default: throw new Fatal("Index of model == %d is invalid. Valid values are 0,1,2",idx_mdl);
    }

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz  v0", -100.0,-100.0,-100.0, 2.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);
    dom.SetOutNods ("labtest", Array<int>((nd==3 ?  7 : 7), /*JustOne*/true));
    dom.SetOutEles ("labtest", Array<int>((nd==3 ? 13 : 0), /*JustOne*/true));

    // debug data
    DbgDat dat;
    dat.node_id = 7;
    dat.dof_idx = dom.Nods[dat.node_id]->FMap("fz");
 
    // solver
    FEM::Solver sol(dom, &DbgFun, &dat);
    sol.TolR   = 1.0e-3;
    //sol.dTini = 0.01;
    //sol.mMax = 2.0;
    //sol.mMin = 0.01;

    ////////////////////////////////////////////////////////////////////////////////////////// Run /////
    
    Dict bcs;
    bcs.Set(-10, "ux", 0.0);
    bcs.Set(-30, "uy", 0.0);
    bcs.Set(-50, "uz", 0.0);
    bcs.Set(-20, "qn", qx);
    bcs.Set(-40, "qn", qy);
    bcs.Set(-60, "uz", uz);
    dom.SetBCs (bcs);
    sol.Solve  (/*NDiv*/20);

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.WriteVTU ("labtest");

    return 0;
}
MECHSYS_CATCH
