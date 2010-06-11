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
#include <mechsys/vtk/win.h>
#include <mechsys/vtk/arrow.h>
#include <mechsys/vtk/axes.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    VTK::Axes ax(2.0, true);

    VTK::Arrow a1(Vec3_t(0.0,0.0,0.0), Vec3_t(1.0,0.0,0.0));
    VTK::Arrow a2(Vec3_t(1.0,1.0,1.0), Vec3_t(1.0,0.0,0.0));
    VTK::Arrow a3(Vec3_t(0.0,0.0,0.0), Vec3_t(0.0,1.0,0.0));
    VTK::Arrow a4(Vec3_t(0.0,0.0,0.0), Vec3_t(0.0,0.0,1.0));
    VTK::Arrow a5(Vec3_t(0.0,0.0,0.0), Vec3_t(1.0,1.0,1.0));
    VTK::Arrow a6(Vec3_t(0.0,0.0,0.0), Vec3_t(2.0,2.0,2.0));
    VTK::Arrow a7(Vec3_t(0.8,0.8,0.2), Vec3_t(0.8,0.2,0.8));

    a6.SetGeometry (0.05, 0.025, 0.01);

    a1.SetColor ("red"    );
    a2.SetColor ("green"  );
    a3.SetColor ("blue"   );
    a4.SetColor ("cyan"   );
    a5.SetColor ("magenta");
    a6.SetColor ("black"  );
    a7.SetColor ("peacock");

    VTK::Win win;
    ax.AddTo (win);
    a1.AddTo (win);
    a2.AddTo (win);
    a3.AddTo (win);
    a4.AddTo (win);
    a5.AddTo (win);
    a6.AddTo (win);
    a7.AddTo (win);
    win.Show ();

    return 0;
}
MECHSYS_CATCH
