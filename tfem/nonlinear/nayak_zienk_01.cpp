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
         of.open ("nayak_zienk_01.res",std::ios::out);
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
    //if (argc>1) is2d = atoi(argv[1]);

    double r  = 3.0;
    double R  = 6.0;
    double th = 4.0*Util::PI/180.0;

    // mesh
    Mesh::Structured mesh(2);
    mesh.GenSector (/*Nr*/4, /*Nth*/2, r, R, th);
    mesh.WriteMPY  ("nayak_zienk_01");

    // parameters
    double E  = 10.0e+6;
    double nu = 0.33;
    double sY = 10.0e+3;
    double Hp = 0.0;

    // props, domain, and solver
    Dict prps, mdls;
    prps.Set(-1, "prob geom axs", PROB("Equilib"), GEOM("Quad8"), 1.0);
    mdls.Set(-1, "name E nu fc sY Hp axs", MODEL("ElastoPlastic"), E, nu, FAILCRIT("VM"), sY, Hp, 1.0);
    //mdls.Set(-1, "name E nu axs", MODEL("LinElastic"), E, nu, 1.0);
    FEM::Domain dom(mesh, prps, mdls, /*inis*/Dict());

    // debug data
    DbgDat dat;
    dat.eqx = 4; // eq for output

    // solver
    FEM::Solver sol(dom, NULL, NULL, &DbgFun, &dat);
    //sol.Scheme = FEM::Solver::FE_t;
    //sol.Scheme = FEM::Solver::NR_t;
    sol.TolR  = 1.0e-2;
    sol.STOL  = 1.0e-7;
    sol.dTini = 0.1;
    sol.SSOut = false;

    // solve
    dom.SetOutNods ("nayak_zienk_01", Array<int>(0,1,2));
    dom.SetOutEles ("nayak_zienk_01", Array<int>(0,1,2));
    Dict bcs;
    bcs.Set      (-100, "inclsupport alpha", 1.0, th);
    bcs.Set      (-30,  "uy", 0.0);
    bcs.Set      (-10,  "qn", -13900.0);
    dom.SetBCs   (bcs);
    sol.Solve    (100);
    dom.WriteVTU ("nayak_zienk_01");

    // end
	return 0;
}
MECHSYS_CATCH
