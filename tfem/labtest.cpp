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
#include <mechsys/models/unconv03.h>
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
    bool   NR      = false;
    bool   pstrain = false;
    if (argc>1) idx_mdl = atoi(argv[1]);
    if (argc>2) NR      = atoi(argv[2]);
    if (argc>3) pstrain = atoi(argv[3]);

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

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz  v0", -100.0,-100.0,-100.0, 2.0);

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
            double cu = 120.0/sqrt(3.0);
            cout << "\n[1;33m======================= Elasto-perfect-plastic with von Mises ======================[0m\n";
            mdls.Set(-1, "name E nu fc cu", MODEL("ElastoPlastic"), 4000.0, 0.3, FAILCRIT("VM"), cu);
            break;
        }
        case 2: // elasto-perfect-plastic with Mohr-Coulomb FC
        {
            cout << "\n[1;33m======================= Elasto-perfect-plastic with Mohr-Coulomb ===================[0m\n";
            mdls.Set(-1, "name E nu fc c phi", MODEL("ElastoPlastic"), 4000.0, 0.3, FAILCRIT("MC"), 0.0, M2Phi(1.0,"cam"));
            break;
        }
        case 3: // linear elastic
        {
            cout << "\n[1;33m==================================== Linear elastic ================================[0m\n";
            mdls.Set(-1, "name E nu", MODEL("LinElastic"), 4000.0, 0.3);
            break;
        }
        case 4:
        {
            cout << "\n[1;33m==================================== Unconv 03 =====================================[0m\n";
            double l0    = 0.001;
            double l1    = 0.005;
            double l3    = 0.008;
            double betb  = 100.0;
            double betbb = 100.0;
            double v0    = 2.0;
            double xR10  = log(150.0*sqrt(3.0));
            double xR30  = xR10+0.1;
            double phi   = M2Phi(1.0,"cam");
            double nu    = 0.3;
            mdls.Set (-1, "name l0 l1 l3 betb betbb phi nu", MODEL("Unconv03"), l0, l1, l3, betb, betbb, phi, nu);
            inis.Set (-1, "sx sy sz v0 xR10 xR30", -100.0,-100.0,-100.0, v0, xR10, xR30);
            break;
        }
        default: throw new Fatal("Index of model == %d is invalid. Valid values are 0,1,2",idx_mdl);
    }

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
    if (NR) sol.SetScheme("NR");
    sol.TolR   = 1.0e-3;
    //sol.dTini = 0.01;
    //sol.mMax = 2.0;
    //sol.mMin = 0.01;

    ////////////////////////////////////////////////////////////////////////////////////////// Run /////
    
    Dict bcs;
    bcs.Set(-10, "ux",  0.0);
    bcs.Set(-30, "uy",  0.0);
    bcs.Set(-50, "uz",  0.0);
    if (pstrain)
    {
        bcs.Set(-20, "ux",  0.0); // delta pressure on the x-side (keep it constant)
        //bcs.Set(-40, "qn",  0.0); // delta pressure on the y-side (keep it constant)
        bcs.Set(-60, "uz", -0.07);
        dom.SetBCs (bcs);
        sol.Solve  (/*NDiv*/20, NULL,NULL,true);
    }
    else
    {
        bool   disp = false;
        double dqn  = 150.0;
        if (disp)
        {
            bcs.Set(-20, "qn",  0.0); // delta pressure on the x-side (keep it constant)
            bcs.Set(-40, "qn",  0.0); // delta pressure on the y-side (keep it constant)
            bcs.Set(-60, "uz", -0.07);
        }
        else
        {
            bcs.Set(-20, "qn",  0.0); // delta pressure on the x-side (keep it constant)
            bcs.Set(-40, "qn",  0.0); // delta pressure on the y-side (keep it constant)
            bcs.Set(-60, "qn", -dqn);
        }
        dom.SetBCs (bcs);
        sol.Solve  (/*NDiv*/20, NULL,NULL,true);

        // TODO: crazy behaviour
        //bcs.Set(-60, "uz", -uz);
        //dom.SetBCs (bcs);
        //sol.Solve  (/*NDiv*/20);

        bcs.clear();
        bcs.Set(-10, "ux", 0.0);
        bcs.Set(-30, "uy", 0.0);
        bcs.Set(-50, "uz", 0.0);
        if (!disp)
        {
            bcs.Set(-20, "qn", 0.0);
            bcs.Set(-40, "qn", 0.0);
            bcs.Set(-60, "qn", dqn);
            dom.SetBCs (bcs);
            sol.Solve  (/*NDiv*/20);
        }

        bcs.clear();
        bcs.Set(-10, "ux", 0.0);
        bcs.Set(-30, "uy", 0.0);
        bcs.Set(-50, "uz", 0.0);
        if (!disp)
        {
            bcs.Set(-20, "qn", 0.0);
            bcs.Set(-40, "qn", 0.0);
            bcs.Set(-60, "qn", -dqn);
            dom.SetBCs (bcs);
            sol.SSOut = true;
            sol.Solve  (/*NDiv*/20);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.WriteVTU ("labtest");

    return 0;
}
MECHSYS_CATCH
