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
#include <mechsys/util/maps.h>
#include <mechsys/vtk/win.h>
#include <mechsys/vtk/axes.h>
#include <mechsys/vtk/sgrid.h>
#include <mechsys/vtk/spheres.h>
#include <mechsys/util/colors.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    // input
    String filename;
    bool   show_ids     = false;
    bool   with_control = true;
    if (argc>1) filename     =      argv[1];
    if (argc>2) with_control = atoi(argv[2]);
    if (argc>3) show_ids     = atoi(argv[3]);
    if (argc<2) throw new Fatal("filename is needed as argument");

    // file key
    String fkey;
    if (with_control) fkey = filename;
    else
    {
        fkey.resize(filename.size()-4);
        for (size_t i=0; i<fkey.size(); ++i) fkey[i] = filename[i];
    }

    // parse file
    int           nout = 1; // number of time output
    Array<int>    N(3);     // number of cells along each axis
    Array<double> L(6);     // limits
    String buf;
    if (with_control)
    {
        buf.Printf("%s_control.res",fkey.CStr());
        std::fstream fi(buf.CStr(), std::ios::in);
        if (!fi.is_open()) throw new Fatal("Could not open file <%s>",buf.CStr());
        String line,key;
        double val;
        while (!fi.eof())
        {
            std::getline (fi,line);
            std::istringstream iss(line);
            iss >> key >> val;
            if      (key=="fkey") {}
            else if (key=="nout") nout = val;
            else if (key=="nx")   N[0] = val;
            else if (key=="ny")   N[1] = val;
            else if (key=="nz")   N[2] = val;
            else if (key=="lxmi") L[0] = val;
            else if (key=="lxma") L[1] = val;
            else if (key=="lymi") L[2] = val;
            else if (key=="lyma") L[3] = val;
            else if (key=="lzmi") L[4] = val;
            else if (key=="lzma") L[5] = val;
            else throw new Fatal("Key %s is wrong",key.CStr());
        }
        fi.close();
    }

    // window
    VTK::Win win;
    win.GetCamera()->SetViewUp     (0,1,0);
    win.GetCamera()->SetPosition   (0.5,0.5,3);
    win.GetCamera()->SetFocalPoint (0.5,0.5,0);
    win.GetCamera()->ParallelProjectionOn();

    // axes
    VTK::Axes axe(1.1);
    axe.XLabel().SetOrientation (0,0,0);
    axe.YLabel().SetOrientation (0,0,0);
    axe.ZLabel().SetOrientation (0,0,0);
    axe.AddTo (win);

    // grid
    if (with_control)
    {
        for (size_t i=0; i<N.Size(); ++i) N[i]++;
        VTK::SGrid grd(N.GetPtr(), L.GetPtr());
        grd.SetColor ("black", 0.2);
        grd.ShowIds  (90,90,45,0.003,8,false);
        grd.AddTo    (win);
    }

    // spheres
    for (int stp_out=0; stp_out<nout; ++stp_out)
    {
        // read data
        if (with_control) buf.Printf("%s_%08d.res", fkey.CStr(), stp_out);
        else              buf = filename;
        Table tab;
        tab.Read (buf.CStr());
        //Array<int>    const & id = tab("id");
        Array<double> const & xc = tab("xc");
        Array<double> const & yc = tab("yc");
        Array<double> const & zc = tab("zc");
        Array<double> const & ra = tab("ra");
        //Array<double> const & vx = tab("vx");
        //Array<double> const & vy = tab("vy");
        //Array<double> const & vz = tab("vz");

        // spheres
        Array<Vec3_t> X(xc.Size());
        for (size_t i=0; i<xc.Size(); ++i) X[i] = xc[i], yc[i], zc[i];
        VTK::Spheres spheres(X, ra);
        if (show_ids) spheres.ShowIds  (0,0,0,0.003,10,false);
        else          spheres.SetColor ("red",1.0);
        spheres.AddTo (win);
        win.Show();
    }

    // end
    return 0;
}
MECHSYS_CATCH
