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

// MechSys
#include <mechsys/vtk/axes.h>
#include <mechsys/vtk/sphere.h>
#include <mechsys/vtk/cube.h>
#include <mechsys/vtk/plane.h>
#include <mechsys/vtk/vtkwin.h>
#include <mechsys/util/fatal.h>

int main(int argc, char **argv) try
{
    Axes   ax(/*scale*/2, /*hydroline*/true);
    Sphere sp(/*x*/  Vec3_t(1,1,1), /*R*/0.5);
    Cube   cu(/*Cen*/Vec3_t(0.5,1.0,0.5), /*lx*/1.0, /*ly*/2.0, /*lz*/0.5);
    Plane  pl(/*Ori*/Vec3_t(0,0,0), /*P1*/Vec3_t(2,0,0), /*P2*/Vec3_t(0,2,0), /*normal*/Vec3_t(0,0,1));
    Plane  p2(/*Cen*/Vec3_t(1,1,1), /*n*/Vec3_t(1,1,1));
    p2.SetColor ("green");

    VTKWin win;
    ax.AddActorsTo (win);
    sp.AddActorsTo (win);
    cu.AddActorsTo (win);
    pl.AddActorsTo (win);
    p2.AddActorsTo (win);
    win.Show();
    win.WritePNG("objects.png");
    cout << "file <[1;34mobjects.png[0m> written" << endl;

    return 0;
}
MECHSYS_CATCH
