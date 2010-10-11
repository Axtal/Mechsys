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
#include <mechsys/vtk/arrow.h>
#include <mechsys/vtk/axes.h>
#include <mechsys/vtk/cube.h>
#include <mechsys/vtk/cylinder.h>
#include <mechsys/vtk/disk.h>
#include <mechsys/vtk/plane.h>
#include <mechsys/vtk/sphere.h>
#include <mechsys/vtk/spheres.h>
#include <mechsys/vtk/text.h>
#include <mechsys/vtk/win.h>
#include <mechsys/util/fatal.h>

int main(int argc, char **argv) try
{
    VTK::Axes     ax(/*scale*/1, /*hydroline*/true);
    VTK::Arrow    ar(/*x*/Vec3_t(0.0,0.0,0.0), /*v*/Vec3_t(0.5,0.5,0.5));
    VTK::Sphere   sp(/*x*/Vec3_t(0.75,0.75,0.75), /*R*/0.25);
    VTK::Cube     cu(/*Cen*/Vec3_t(0.0,0.0,0.0), /*lx*/1.0, /*ly*/1.0, /*lz*/1.0);
    VTK::Cylinder cy(/*X0*/Vec3_t(0.0,0.0,0.0), /*X1*/Vec3_t(1.0,0.0,0.0), /*R*/0.1, /*Cap*/false, /*Res*/30);
    VTK::Disk     di(/*X0*/Vec3_t(0.0,0.0,0.0), /*X2*/Vec3_t(1.0,1.0,1.0), /*Rin*/0.1, /*Rout*/0.2, /*Cap*/false);
    VTK::Disk     d1(/*X0*/Vec3_t(1.0,0.0,0.0), /*X2*/Vec3_t(1.0,1.0,1.0), /*Rin*/0.1, /*Rout*/0.2, /*Cap*/false);
    VTK::Disk     d2(/*X0*/Vec3_t(0.0,1.0,0.0), /*X2*/Vec3_t(1.0,1.0,1.0), /*Rin*/0.1, /*Rout*/0.2, /*Cap*/false);
    VTK::Disk     d3;
    VTK::Text     tx(/*X*/Vec3_t(0.8,0.3,0.5), "VTK");
    VTK::Plane    pl(/*Ori*/Vec3_t(0,0,0), /*P1*/Vec3_t(2,0,0), /*P2*/Vec3_t(0,2,0), /*normal*/Vec3_t(0,0,1));
    VTK::Plane    p2(/*Cen*/Vec3_t(1,1,1), /*n*/Vec3_t(1,1,1));

    p2.SetColor ("green");
    cu.SetColor ("peacock", 0.1);
    di.SetColor ("blue",    0.3);
    d1.SetColor ("blue",    0.3);
    d2.SetColor ("blue",    0.3);

    d3.SetRadiusIn  (0.05);
    d3.SetRadiusOut (0.08);
    d3.SetPoints    (/*X0*/Vec3_t(2.0,0.0,0.0), /*X2*/Vec3_t(0.0,2.0,0.0));

    Array<Vec3_t> X(3);
    Array<double> R(3);
    X[0] = 0.0, 0.0, 1.0;
    X[1] = 1.0, 0.0, 1.0;
    X[2] = 0.0, 1.0, 1.0;
    R    = 0.1, 0.2, 0.3;
    VTK::Spheres ss(X, R);
    ss.SetCenter (0, Vec3_t(1,1,1));
    ss.SetRadius (0, 0.05);

    Array<double> x(3),y(3),z(3);
    x = 0.0, 1.0, 0.0;
    y = 0.0, 0.0, 1.0;
    z = 0.5, 0.5, 0.5;
    VTK::Spheres s2;
    s2.SetSpheres (x, y, z, &R);
    s2.SetColor   ("aquamarine", 0.5);
    s2.ShowIds    ();

    VTK::Win win;
    ax.AddTo (win);
    ar.AddTo (win);
    sp.AddTo (win);
    cu.AddTo (win);
    cy.AddTo (win);
    di.AddTo (win);
    d1.AddTo (win);
    d2.AddTo (win);
    d3.AddTo (win);
    tx.AddTo (win);
    pl.AddTo (win);
    p2.AddTo (win);
    ss.AddTo (win);
    s2.AddTo (win);
    win.Show();
    win.WritePNG("objects.png");
    cout << "file <[1;34mobjects.png[0m> written" << endl;

    return 0;
}
MECHSYS_CATCH
