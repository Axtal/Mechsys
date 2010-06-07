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
#include <mechsys/vtk/vtkwin.h>
#include <mechsys/vtk/arrow.h>
#include <mechsys/vtk/axes.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    Axes ax(2.0, true);

    Arrow a1(Vec3_t(0.0,0.0,0.0), Vec3_t(1.0,0.0,0.0));
    Arrow a2(Vec3_t(1.0,1.0,1.0), Vec3_t(1.0,0.0,0.0));
    Arrow a3(Vec3_t(0.0,0.0,0.0), Vec3_t(0.0,1.0,0.0));
    Arrow a4(Vec3_t(0.0,0.0,0.0), Vec3_t(0.0,0.0,1.0));
    Arrow a5(Vec3_t(0.0,0.0,0.0), Vec3_t(1.0,1.0,1.0));
    Arrow a6(Vec3_t(0.0,0.0,0.0), Vec3_t(2.0,2.0,2.0), /*BodRad*/0.005, /*TipRad*/0.01, /*TipLen*/0.50, /*Resolution*/18);
    Arrow a7(Vec3_t(0.5,0.5,0.5), Vec3_t(0.5,0.5,0.5));

    a1.SetColor ("red"    );
    a2.SetColor ("green"  );
    a3.SetColor ("blue"   );
    a4.SetColor ("cyan"   );
    a5.SetColor ("magenta");
    a6.SetColor ("black"  );
    a7.SetColor ("peacock");

    VTKWin win;
    ax.AddActorsTo (win);
    a1.AddActorsTo (win);
    a2.AddActorsTo (win);
    a3.AddActorsTo (win);
    a4.AddActorsTo (win);
    a5.AddActorsTo (win);
    a6.AddActorsTo (win);
    a7.AddActorsTo (win);
    win.Show       ();

    return 0;
}
MECHSYS_CATCH
