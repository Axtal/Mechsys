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
    int  idx_out;
    std::ofstream of;
    ~DbgDat () { of.close (); }
     DbgDat () : idx_out(0)
    {
         of.open ("nayak_zienk_01.res",std::ios::out);
         of<<_6_3<<"Time"<<_8s<<"u"<< _8s<<"fint"<<_8s<<"fext\n";
         of<<_6_3<<  0   <<_8s<< 0 << _8s<<  0   <<_8s<<  0 << endl;
    }
};

void DbgFun (FEM::Solver const & Sol, void * Dat)
{
    DbgDat * dat = static_cast<DbgDat*>(Dat);
    dat->of << _6_3 << Sol.Dom.Time << _8s << Sol.U(dat->eqx) << _8s << Sol.F_int(dat->eqx) << _8s << Sol.F(dat->eqx) << endl;

    /*
    String fn;
    fn.Printf ("nayak_zienk_01_%04d", dat->idx_out);
    Sol.Dom.WriteMPY (fn.CStr());
    dat->idx_out++;
    */
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
    //mesh.WriteMPY  ("nayak_zienk_01");

    // parameters
    double E  = 10.0e+6;
    double nu = 0.33;
    double sY = 10.0e+3;
    double Hp = 0.0;

    // props, domain, and solver
    Array<int> out_verts(0,1,2,17);
    Dict prps, mdls;
    prps.Set(-1, "prob geom axs", PROB("Equilib"), GEOM("Quad8"), 1.0);
    mdls.Set(-1, "name E nu fc sY Hp axs", MODEL("ElastoPlastic"), E, nu, FAILCRIT("VM"), sY, Hp, 1.0);
    //mdls.Set(-1, "name E nu axs", MODEL("LinElastic"), E, nu, 1.0);
    FEM::Domain dom(mesh, prps, mdls, /*inis*/Dict(), "nayak_zienk_01", &out_verts);

    // debug data
    DbgDat dat;
    dat.eqx = 4; // eq for output

    // solver
    FEM::Solver sol(dom, NULL, NULL, &DbgFun, &dat);
    //sol.SetScheme ("NR");
    sol.SSOut = true;
    //sol.SetIncsW (40, /*NonLinWei*/true);

    // solve
    Dict bcs;
    bcs.Set      (-100, "incsup alpha", 1.0, 4.0);
    bcs.Set      (-30,  "uy", 0.0);
    bcs.Set      (-10,  "qn", -13920.0);
    dom.SetBCs   (bcs);
    cout << dom << endl;
    sol.Solve    (40);
    //sol.Solve    (40, "nayak_zienk_01");
    dom.WriteVTU ("nayak_zienk_01");

    /* It seems that as the plastification zone advances, the two matrices Ke and Kep
     * become different. Each time a new plastification happens in a subincrement, the
     * new K is different, then the subincrement must be subdivided.
     *
     * After a first plastification occurs, the two matrices FE and ME will be similar
     * and then the update is faster. If plastification happens inside the subincrement
     * it will be slow.
     *
     * This my be the reason why ME uses many NSS intercalated with a few. Ex.: nss=1
     * nss=1 nss=1 nss=48 nss=1 nss=3 nss=1 nss=29 ...
     *
     * Thus, it seems that the ME is capturing the moment of plastification */

    // end
	return 0;
}
MECHSYS_CATCH
