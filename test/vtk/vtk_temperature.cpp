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
#include <mechsys/vtk/sgrid.h>

using std::cout;
using std::endl;
using Util::PI;
using Util::SQ3;

void T (Vec3_t const & X, double & F, Vec3_t & V, void * UserData)
{
    F = X(0);
}

int main(int argc, char **argv) try
{
    // number:  nx ny nz
    Array<int> N(10, 5, 2);
    if (argc>1) N[0] = atoi(argv[1]);
    if (argc>2) N[1] = atoi(argv[2]);
    if (argc>3) N[2] = atoi(argv[3]);
    if (N[0]<2) throw new Fatal("nx must be greater than 1");
    if (N[1]<2) throw new Fatal("ny must be greater than 1");
    if (N[2]<2) throw new Fatal("nz must be greater than 1");

    // limits
    Array<double> L(0,10, 0,5, 0,1); // xmin,xmax  ymin,ymax  zmin,zmax

    // grid
    //VTK::SGrid gri(N, L, T);  bool use_T = true;  // <<<<  Using Function
    VTK::SGrid gri(N, L);       bool use_T = false; // <<<<  Setting each point

    // Show surface
    gri.ShowSurface ();

    // Seting points
    if (!use_T)
    {
        for (int k=0; k<N[2]; ++k)
        for (int j=0; j<N[1]; ++j)
        for (int i=0; i<N[0]; ++i)
        {
            if (i>3 && i<6 && j>1 && j<3) gri.SetF (i,j,k, 20);
            else                          gri.SetF (i,j,k, i+j);
        }
        gri.RescaleCMap();
    }

    // Set colormap
    gri.SetCMap ("Rainbow");
    
    // window and axes
    VTK::Win  win;
    win.Camera (/*xup*/0, /*yup*/1, /*zup*/0, /*xfoc*/L[0]/2, /*yfoc*/L[1]/2, /*zfoc*/0, 
                                              /*xpos*/L[0]/2, /*ypos*/L[1]/2, /*zpos*/1); 
    win.Parallel(true);
    gri.AddTo (win);
    win.Show  ();

    // end
    return 0;
}
MECHSYS_CATCH
