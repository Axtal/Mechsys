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
#include <mechsys/vtk/isosurf.h>
#include <mechsys/vtk/vtkwin.h>
#include <mechsys/vtk/axes.h>
#include <mechsys/vtk/meshgrid.h>
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
    // input
    size_t ndiv = 20;
    if (argc>1) ndiv = atoi(argv[1]);

    // meshgrid
    MeshGrid mg(0.0,1.0,ndiv, 0.0,1.0,ndiv, 0.0,1.0,ndiv);

    // data
    Data dat;
    dat.c = 0.5, 0.5, 0.5;
    dat.r = 0.25;

    // iso surface
    IsoSurf iso(mg, &Func, &dat);
    iso.ShowVectors = true;
    iso.SetVecScale (0.1);
    iso.SetVecF     (0.0, /*Tol*/1.0e-2);

    // axes
    Axes ax(1.2, /*HydroLine*/true);

    // window
    VTKWin win;
    iso.AddActorsTo (win);
    ax .AddActorsTo (win);
    win.Show();

    // end
    return 0;
}
MECHSYS_CATCH
