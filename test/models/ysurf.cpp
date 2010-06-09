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
#include <mechsys/models/elastoplastic.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/fatal.h>
#include <mechsys/vtk/meshgrid.h>
#include <mechsys/vtk/sphere.h>
#include <mechsys/vtk/axes.h>
#include <mechsys/vtk/isosurf.h>

using std::cout;
using std::endl;
using Util::PI;

void Func (Vec3_t X, double & F, Vec3_t & V, void * UserData)
{
    ElastoPlastic const * mdl = static_cast<ElastoPlastic*>(UserData);

    EquilibState sta(/*NDim*/3);
    sta.Sig = X(0), X(1), X(2), 0.0, 0.0, 0.0;
    F = mdl->YieldFunc (&sta);

    cout << F << endl;

    V = 1.0, 0.0, 0.0;
}

int main(int argc, char **argv) try
{
    // input
    int np = 10;
    if (argc>1) np = atoi(argv[1]);

    // model
    SDPair prms;
    prms.Set("E nu fc c phi", 1000.0, 0.3, FAILCRIT("MC"), 1.0, 30.0);
    ElastoPlastic mdl(/*NDim*/3, prms);

    // grid
    MeshGrid oct(0.0,  10.0,   np,  // p
                 0.0,  10.0,   np,  // q
                 -PI,    PI, 2*np); // th
    MeshGrid sig(0,0,np, 0,0,np, 0,0,2*np);
    for (int i=0; i<oct.Length(); ++i)
    {
        Vec3_t L;
        pqth2L (oct.X[i], oct.Y[i], oct.Z[i], L, "cam");
        sig.X[i] = L(0);
        sig.Y[i] = L(1);
        sig.Z[i] = L(2);
    }

    // window and axes
    VTKWin win;
    win.SetViewDefault (/*reverse*/true);
    Axes ax(/*scale*/20, /*hydroline*/true, /*reverse*/true);
    ax.AddActorsTo (win);

    // isosurf
    IsoSurf iso(sig, &Func, &mdl);
    iso.ShowPoints  = true;
    iso.ShowVectors = false;
    iso.AddActorsTo (win);
    iso.WriteFile  ("ysurf.vtk");

    // end
    win.Show ();
    return 0;
}
MECHSYS_CATCH
