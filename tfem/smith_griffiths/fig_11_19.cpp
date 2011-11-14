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

/*  Smith & Griffiths (2004): Figure 11.19 p500 
 *  ===========================================  */

// MechSys
#include <mechsys/fem/fem.h>
#include <mechsys/inpfile.h>

int main(int argc, char **argv) try
{
    // init
    String inpfn("fig_11_19.inp");
    String matfn("materials.mat");
    bool   verbose  = true;
    bool   forcegty = true;
    INIT_MAT_INP (argc, argv, inpfn, matfn, verbose, forcegty,   mat, inp);

    // fnkey
    if (inp.rk) inp.fnkey += "_rk";

    // mesh
    size_t nc = 6;             // number of cells
    size_t nv = 3*(nc+1)+2*nc; // number of vertices
    double H  = 1.0;           // height
    double h  = 0.5;           // half height
    double L  = 15.0;          // length
    double l  = L/nc;          // element length
    double m  = l/2.;          // half length (for mid nodes)
    double x  = 0.0;
    size_t k  = 0;
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize (nv/*verts*/, nc/*cells*/);
    mesh.SetVert (k,    0, x, H  );  k++;
    mesh.SetVert (k,    0, x, h  );  k++;
    mesh.SetVert (k, -100, x, 0.0);  k++;
    for (size_t i=0; i<nc; ++i)
    {
        int tag = (i==nc-1 ? -200 : 0);
        mesh.SetVert (k,    0, x+m, H  );  k++;
        mesh.SetVert (k, -200, x+m, 0.0);  k++;
        mesh.SetVert (k,  tag, x+l, H  );  k++;
        mesh.SetVert (k,  tag, x+l, h  );  k++;
        mesh.SetVert (k, -200, x+l, 0.0);  k++;
        mesh.SetCell (i, -1, Array<int>(i*5, i*5+2, i*5+7, i*5+5,  i*5+1, i*5+4, i*5+6, i*5+3));
        mesh.SetBryTag (i, 3, -10);
        x += l;
    }
    mesh.WriteMPY (inp.fnkey.CStr());
    
    // allocate and solve
    FEM_ALLOC_DOMAIN(mesh, mat, inp, dom);
    FEM_ALLOC_SOLVER(inp, dom, sol);
    FEM_SOLVE(verbose, inp, sol, dom);

    // end
    return 0.0;
}
MECHSYS_CATCH
