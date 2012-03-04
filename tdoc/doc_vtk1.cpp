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
#include <mechsys/vtk/sphere.h>
#include <mechsys/vtk/win.h>
#include <mechsys/util/fatal.h>

int main(int argc, char **argv) try {
    VTK::Axes   ax(/*scale*/1, /*hydroline*/true);
    VTK::Arrow  ar(/*x*/Vec3_t(0.0,0.0,0.0), /*v*/Vec3_t(0.5,0.5,0.5));
    VTK::Sphere sp(/*x*/Vec3_t(0.75,0.75,0.75), /*R*/0.25);
    VTK::Cube   cu(/*Cen*/Vec3_t(0.0,0.0,0.0), /*lx*/1.0, /*ly*/1.0, /*lz*/1.0);
    VTK::Win win;
    ax.AddTo (win);
    ar.AddTo (win);
    sp.AddTo (win);
    cu.AddTo (win);
    win.Show();
} MECHSYS_CATCH
