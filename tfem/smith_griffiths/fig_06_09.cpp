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
#include <cmath>

// MechSys
#include <mechsys/fem/fem.h>

using std::cout;
using std::endl;
using Util::_6_3;
using Util::_8s;
using FEM::PROB;
using FEM::GEOM;
const double TRUE = 1.0;

double DelP = 520.0;

struct OutDat
{
    long eqx; // equation number for output
    std::ofstream of;
    double prev_time;
    ~OutDat () { of.close (); }
     OutDat (char const * FileKey)
    {
        char buf[64];
        sprintf (buf,"%s.res",FileKey);
        of.open (buf,std::ios::out);
        of<<_6_3<<"Time"<<_8s<<"u"<< _8s<<"fint"<<_8s<<"fext\n";
    }
};

void OutFun (FEM::Solver const & Sol, void * Dat)
{
    OutDat * dat = static_cast<OutDat*>(Dat);
    dat->of << _6_3 << Sol.Dom.Time << _8s << Sol.U(dat->eqx) << _8s << Sol.F_int(dat->eqx) << _8s << Sol.F(dat->eqx) << endl;
}

int main(int argc, char **argv) try
{
    // input
    char const * scheme = "ME";
    if (argc>1) scheme = argv[1];
    char filekey[128];
    sprintf (filekey,"fig_06_09_%s",scheme);

    // mesh
    Array<double> X(0.0, 1.0, 2.0, 3.0, 4.0, 5.5, 7.0, 9.0, 12.0);
    Array<double> Y(0.0, -1.25, -2.5, -3.75, -5.0);
    Mesh::Generic mesh(/*NDim*/2);
    mesh.GenGroundSG (X, Y, /*length_of_footing*/2.0);
    mesh.WriteMPY    (filekey);

    // parameters
    double E  = 1.0e+5;
    double nu = 0.3;
    double cu = 100.0;
    double Hp = 0.0;
    
    // props, models, and domain
    Dict prps, mdls;
    prps.Set(-1, "prob geom psa nip", PROB("Equilib"), GEOM("Quad8"), TRUE, 4.0);
    mdls.Set(-1, "name E nu VM cu Hp psa", MODEL("ElastoPlastic"), E, nu, TRUE, cu, Hp, TRUE);
    //mdls.Set(-1, "name E nu axs", MODEL("LinElastic"), E, nu, 1.0);
    FEM::Domain dom(mesh, prps, mdls, /*inis*/Dict());

    // debug data
    OutDat dat(filekey);
    dat.eqx = 1; // eq for output

    // solver
    FEM::Solver sol(dom, &OutFun, &dat);
    sol.SetScheme (scheme);
    if (strcmp(scheme,"FE")==0) sol.nSS = 1000;
    //sol.STOL = 1.0e+9;

    // weights
    sol.IncsW.Resize (10);
    sol.IncsW = 200., 100., 50., 50., 50., 30., 20., 10., 5., 5.;
    for (size_t i=0; i<sol.IncsW.Size(); ++i) sol.IncsW[i] /= DelP;

    // solve
    Dict bcs;
    bcs.Set      (-10, "ux",    0.0);
    bcs.Set      (-20, "ux",  0.0);
    bcs.Set      (-30, "ux uy", 0.0, 0.0);
    bcs.Set      (-40, "qn",    -DelP);
    dom.SetBCs   (bcs);
    sol.Solve    (sol.IncsW.Size());
    dom.WriteVTU (filekey);
}
MECHSYS_CATCH
