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
using Util::_6_3;
using Util::_8s;

struct DbgDat
{
    long eqx; // equation number for output
    std::ofstream of;
    ~DbgDat () { of.close (); }
     DbgDat ()
    {
         of.open ("ep_test.res",std::ios::out);
         of<<_6_3<<"Time"<<_8s<<"u"<< _8s<<"fint"<<_8s<<"fext\n";
         of<<_6_3<<  0   <<_8s<< 0 << _8s<<  0   <<_8s<<  0 << endl;
    }
};

void DbgFun (FEM::Solver const & Sol, void * Dat)
{
    DbgDat * dat = static_cast<DbgDat*>(Dat);
    dat->of << _6_3 << Sol.Dom.Time << _8s << Sol.U(dat->eqx) << _8s << Sol.F_int(dat->eqx) << _8s << Sol.F(dat->eqx) << endl;
}

int main(int argc, char **argv) try
{
    // mesh
    Mesh::Structured mesh(2);
    mesh.SetSize   (4, 1);
    mesh.SetVert   (0, -100, 0.0, 0.0);
    mesh.SetVert   (1, -200, 1.0, 0.0);
    mesh.SetVert   (2, -300, 1.0, 1.0);
    mesh.SetVert   (3, -400, 0.0, 1.0);
    mesh.SetCell   (0, -1, Array<int>(0,1,2,3));
    mesh.SetBryTag (0, 2, -10);
    mesh.WriteMPY  ("ep_test");

    // parameters
    double E  = 10.0e+6;
    double nu = 0.33;
    double sY = 10.0e+3;
    double Hp = 0.0;

    // props, domain, and solver
    Array<int> out_verts(1,2,3);
    Dict prps, mdls;
    prps.Set(-1, "prob geom psa", PROB("Equilib"), GEOM("Quad4"), 1.0);
    mdls.Set(-1, "name E nu fc sY Hp psa", MODEL("ElastoPlastic"), E, nu, FAILCRIT("VM"), sY, Hp, 1.0);
    //mdls.Set(-1, "name E nu axs", MODEL("LinElastic"), E, nu, 1.0);
    FEM::Domain dom(mesh, prps, mdls, /*inis*/Dict(), "ep_test", &out_verts);

    // debug data
    DbgDat dat;
    dat.eqx = 5; // eq for output

    // solver
    FEM::Solver sol(dom, NULL, NULL, &DbgFun, &dat);
    sol.Scheme = FEM::Solver::NR_t;

    // solve
    Dict bcs;
    bcs.Set      (-100, "ux uy", 0.0, 0.0);
    bcs.Set      (-200, "uy",    0.0);
    bcs.Set      (-400, "ux",    0.0);
    bcs.Set      (-10,  "qn", -11500.0);
    dom.SetBCs   (bcs);
    sol.Solve    (10);
    dom.WriteVTU ("ep_test");

    // end
	return 0;
}
MECHSYS_CATCH
