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
#include <mechsys/vtk/axes.h>
#include <mechsys/vtk/arrows.h>
#include <mechsys/vtk/sgrid.h>
#include <mechsys/vtk/isosurf.h>

using std::cout;
using std::endl;
using Util::PI;
using Util::SQ3;

void Func (Vec3_t X, double & F, Vec3_t & V, void * UserData)
{
    F = 0;
    if (X(2)>0.1) V = 0.0, 0.0, -0.5;
    else          V = 0.0, 0.0, 0.5;
}

int main(int argc, char **argv) try
{
    // number:  nx ny nz
    Array<int> N(3, 3, 2);
    if (argc>1) N[0] = atoi(argv[1]);
    if (argc>2) N[1] = atoi(argv[2]);
    if (argc>3) N[2] = atoi(argv[3]);
    if (N[0]<2) throw new Fatal("nx must be greater than 1");
    if (N[1]<2) throw new Fatal("ny must be greater than 1");
    if (N[2]<2) throw new Fatal("nz must be greater than 1");

    // limits
    Array<double> L(0,1, 0,1, 0,1);

    // grid
    VTK::SGrid gri(N, L, &Func);
    gri.ShowPoints ();

    // arrows
    VTK::Arrows arr(gri);

    // window and axes
    VTK::Win  win;
    VTK::Axes axe(/*scale*/1.1);
    axe.AddTo (win);
    gri.AddTo (win);
    arr.AddTo (win);
    win.Show  ();

    // end
    return 0;
}
MECHSYS_CATCH
