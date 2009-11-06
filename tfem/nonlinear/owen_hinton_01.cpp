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

/*  Owen & Hinton (1980): Example 3.11.3, p78  *
 *  Finite Elements in Plasticity              *
 *  =========================================  */

// STL
#include <iostream>

// MechSys
#include "mesh/mesh.h"
#include "fem/rod.h"
#include "fem/nlrod.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "util/maps.h"
#include "util/fatal.h"

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;
using Util::_6_3;
using Util::_8s;

struct DbgDat
{
    long eqx; // equation number for output
    std::ofstream of;
    ~DbgDat () { of.close (); }
     DbgDat ()
    {
         of.open ("owen_hinton_01.res",std::ios::out);
         of<<_6_3<<"Time"<<_8s<<"u"<< _8s<<"fint"<<_8s<<"fext\n";
         of<<_6_3<<  0   <<_8s<< 0 << _8s<<  0   <<_8s<<  0 << endl;
    }
};

void DbgFun (FEM::Solver const & Sol, void * Dat)
{
    DbgDat * dat = static_cast<DbgDat*>(Dat);
    dat->of << _6_3 << Sol.Time << _8s << Sol.U(dat->eqx) << _8s << Sol.F_int(dat->eqx) << _8s << Sol.F(dat->eqx) << endl;
}

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    double l = 5.;
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize  (3/*nodes*/, 2/*cells*/);
    mesh.SetVert  (0, -100,  0.0,  2.*l);
    mesh.SetVert  (1, -200,  0.0,     l);
    mesh.SetVert  (2, -300,  0.0,   0.0);
    mesh.SetCell  (0,   -1, Array<int>(0, 1));
    mesh.SetCell  (1,   -1, Array<int>(1, 2));
    mesh.WriteMPY ("owen_hinton_01");

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob E0 alp A fra", PROB("NLRod"), 200.0, 5.0, 1.0, TRUE);
    //prps.Set(-1, "prob E A fra", PROB("Rod"), 200.0, 1.0, TRUE);

    // domain
    FEM::Domain dom(mesh, prps, /*mdls*/Dict(), /*inis*/Dict());

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set(-100, "ux uy", 0.0,0.0);
    bcs.Set(-200, "ux",    0.0);
    bcs.Set(-300, "ux fy", 0.0, -10.0);
    dom.SetBCs (bcs);

    // weights
    Array<double> weights(2);
    weights = 0.8, 0.2;

    // debug data
    DbgDat dat;
    dat.eqx = 5; // eq for output

    // solver
    FEM::Solver sol(dom, &DbgFun, &dat);
    sol.Scheme = FEM::Solver::NR_t;
    sol.TolR   = 1.0e-5;

    // solve
    sol.Solve (weights.Size(), &weights);

    return 0;
}
MECHSYS_CATCH
