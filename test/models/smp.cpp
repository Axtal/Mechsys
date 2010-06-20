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
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/fatal.h>
#include <mechsys/vtk/axes.h>
#include <mechsys/vtk/arrows.h>
#include <mechsys/vtk/sgrid.h>

using std::cout;
using std::endl;
using Util::PI;
using Util::SQ3;

void FuncOCT (Vec3_t X, double & F, Vec3_t & V, void * UserData)
{
    F = 0;
    V = -1.0/SQ3, -1.0/SQ3, -1.0/SQ3;
}

void FuncSMP (Vec3_t X, double & F, Vec3_t & V, void * UserData)
{
    Vec3_t noct(-1.0/SQ3, -1.0/SQ3, -1.0/SQ3);
    F = 0;
    //if (X(0)<0.0 && X(1)<0.0 && X(2)<0.0)
    {
        double I2  = X(0)*X(1) + X(1)*X(2) + X(2)*X(0);
        double I3  = X(0)*X(1)*X(2);
        double cof = sqrt(-I3/I2);
        //V = -cof/sqrt(-X(0)), -cof/sqrt(-X(1)), -cof/sqrt(-X(2));
        //V = -1.0/sqrt(-X(0)), -1.0/sqrt(-X(1)), -1.0/sqrt(-X(2));
        V = 1.0/X(0), 1.0/X(1), 1.0/X(2);
        //V = 1.0/(1.0+X(0)), 1.0/(1.0+X(1)), 1.0/(1.0+X(2));
        V = V/norm(V);

        Vec3_t p(dot(V,noct)*noct);
        V = V - p;

    }
    //else V = 0.0, 0.0, 0.0;
}

int main(int argc, char **argv) try
{
    // number:  nx ny nz
    Array<int> N(5, 3, 51);
    double scale = 6.0;
    if (argc>1) N[0]  = atoi(argv[1]);
    if (argc>2) N[1]  = atoi(argv[2]);
    if (argc>3) N[2]  = atoi(argv[3]);
    if (argc>4) scale = atof(argv[4]);
    if (N[0]<2) throw new Fatal("nx must be greater than 1");
    if (N[1]<2) throw new Fatal("ny must be greater than 1");
    if (N[2]<2) throw new Fatal("nz must be greater than 1");

    // limits
    Array<double> L(0.1,10,  0,15,  -PI,PI);

    // grid
    VTK::SGrid gri_oct(N, L);
    VTK::SGrid gri_smp(N, L);
    gri_oct.ShowPoints ();
    gri_smp.ShowPoints ();

    // rotate
    Vec3_t x, l;
    for (int i=0; i<gri_oct.Size(); ++i)
    {
        gri_oct.GetPoint (i, x);
        pqth2L           (x(0), x(1), x(2), l, "cam");
        gri_oct.SetPoint (i, l);
        gri_smp.SetPoint (i, l);
    }

    // set function
    gri_oct.SetFunc (&FuncOCT);
    gri_smp.SetFunc (&FuncSMP);

    // arrows
    VTK::Arrows arr_oct(gri_oct);
    VTK::Arrows arr_smp(gri_smp);
    arr_oct.SetScale (scale);
    arr_smp.SetScale (scale);
    arr_smp.SetColor ("red");

    // window and axes
    VTK::Win  win;
    VTK::Axes axe(/*scale*/20, /*hydroline*/true, /*reverse*/true);
    win.SetViewDefault (/*reverse*/true);
    axe    .AddTo (win);
    gri_oct.AddTo (win);
    arr_oct.AddTo (win);
    gri_smp.AddTo (win);
    arr_smp.AddTo (win);
    win    .Show  ();

    // end
    return 0;
}
MECHSYS_CATCH
