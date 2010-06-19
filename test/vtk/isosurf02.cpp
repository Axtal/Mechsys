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

// Std Lib
#include <iostream>
#include <cmath>

// MechSys
#include <mechsys/vtk/sgrid.h>
#include <mechsys/vtk/isosurf.h>
#include <mechsys/vtk/win.h>
#include <mechsys/vtk/axes.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;

struct Data
{
    Vec3_t c; // center of sphere
    double r; // radius of sphere
};

void Func (Vec3_t X, double & F, Vec3_t & V, void * UserData)
{
    Data const * d = static_cast<Data*>(UserData);

    F =   pow(X(0) - d->c(0), 2.0)
        + pow(X(1) - d->c(1), 2.0)
        + pow(X(2) - d->c(2), 2.0)
        - pow(d->r,           2.0);

    V = 2.0*(X(0) - d->c(0)),
        2.0*(X(1) - d->c(1)),
        2.0*(X(2) - d->c(2));
}

int main(int argc, char **argv) try
{
    // number:  nx ny nz
    Array<int> N(11, 11, 21);
    if (argc>1) N[0] = atoi(argv[1]);
    if (argc>2) N[1] = atoi(argv[2]);
    if (argc>3) N[2] = atoi(argv[3]);
    if (N[0]<2) throw new Fatal("nx must be greater than 1");
    if (N[1]<2) throw new Fatal("ny must be greater than 1");
    if (N[2]<2) throw new Fatal("nz must be greater than 1");

    // limits
    Array<double> L(6);
    L = 0,1, 0,1, 0,1;

    // data
    Data dat;
    dat.c = 0.5, 0.5, 0.5;
    dat.r = 0.25;

    // sgrid
    VTK::SGrid gri(N,L, &Func, &dat);
    gri.FilterV (0.0, /*Tol*/1.0e-2);

    // iso surface
    VTK::IsoSurf iso(gri);
    iso.ShowVectors = true;
    iso.SetVecScale (0.1);

    // axes
    VTK::Axes axe(1.2, /*HydroLine*/true);

    // window
    VTK::Win win;
    iso.AddTo (win);
    axe.AddTo (win);
    win.Show  ();

    // end
    return 0;
}
MECHSYS_CATCH
