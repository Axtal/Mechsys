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
#include <mechsys/vtk/cube.h>
#include <mechsys/util/colors.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    // input
    String filename;
    bool   show_ids     = true;
    bool   with_control = true;
    bool   shadow       = true;
    if (argc>1) filename     =      argv[1];
    if (argc>2) with_control = atoi(argv[2]);
    if (argc>3) show_ids     = atoi(argv[3]);
    if (argc>4) shadow       = atoi(argv[4]);
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
    int            nout = 1; // number of time output
    Array<int>     N(3);     // number of cells along each axis
    Array<double>  L(6);     // limits
    Array<int>     proc;     // processor ID of each cell
    Array<String>  clrs;     // colors of each processor
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
            if (iss >> key >> val)
            {
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
                else if (key=="proc")
                {
                    if (proc.Size()>0) throw new Fatal("Error with input file");
                    int ncells = val;
                    proc.Resize (ncells);
                    int max_proc = 0;
                    for (int i=0; i<ncells; ++i)
                    {
                        iss >> proc[i];
                        if (proc[i]>max_proc) max_proc = proc[i];
                    }
                    int nprocs = max_proc+1;
                    Colors::GetRandom (nprocs, clrs);
                }
                else throw new Fatal("Key %s is wrong",key.CStr());
            }
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
        Array<int> Nlines(3);
        for (size_t i=0; i<3; ++i) Nlines[i] = N[i]+1;
        VTK::SGrid grd(Nlines.GetPtr(), L.GetPtr());
        grd.SetColor ("black", 0.2);
        //grd.ShowIds  (90,90,45,0.003,8,false);
        grd.AddTo    (win);
        if (true)
        {
            double dx = (L[1]-L[0])/static_cast<double>(N[0]);
            double dy = (L[3]-L[2])/static_cast<double>(N[1]);
            double dz = (L[5]-L[4])/static_cast<double>(N[2]);
            int    nx = N[0];
            int    ny = N[1];
            int    nc = N[0]*N[1]*N[2];
            for (int n=0; n<nc; ++n)
            {
                int    K    =  n / (nx*ny);
                int    J    = (n % (nx*ny)) / nx;
                int    I    = (n % (nx*ny)) % nx;
                double xmin = L[0] +  I   *dx;
                double xmax = L[0] + (I+1)*dx;
                double ymin = L[2] +  J   *dy;
                double ymax = L[2] + (J+1)*dy;
                double zmin = L[4] +  K   *dz;
                double zmax = L[4] + (K+1)*dz;
                Vec3_t cen((xmin+xmax)/2.0, (ymin+ymax)/2.0, (zmin+zmax)/2.0);
                VTK::Cube cube(cen, xmax-xmin, ymax-ymin, zmax-zmin);
                //printf("n=%d, proc[n]=%d, clrs[proc[n]]=%s\n",n,proc[n],clrs[proc[n]].CStr());
                cube.SetColor (clrs[proc[n]].CStr(), 0.1);
                cube.AddTo    (win);
            }
        }
    }

    // spheres
    Array<VTK::Spheres*> sph;
    for (int stp_out=0; stp_out<nout; ++stp_out)
    {
        // read data
        if (with_control) buf.Printf("%s_%08d.res", fkey.CStr(), stp_out);
        else              buf = filename;
        Table tab;
        tab.Read (buf.CStr());
        Array<double> const & id = tab("id");
        Array<double> const & xc = tab("xc");
        Array<double> const & yc = tab("yc");
        Array<double> const & zc = tab("zc");
        Array<double> const & ra = tab("ra");
        //Array<double> const & vx = tab("vx");
        //Array<double> const & vy = tab("vy");
        //Array<double> const & vz = tab("vz");

        Array<int> ids(id.Size());
        for (size_t i=0; i<ids.Size(); ++i) ids[i] = static_cast<int>(id[i]);

        // spheres
        Array<Vec3_t> X(xc.Size());
        for (size_t i=0; i<xc.Size(); ++i) X[i] = xc[i], yc[i], zc[i];
        if (shadow)
        {
            if (sph.Size()>0) sph.Last()->SetColor ("black",0.1);
            sph.Push (new VTK::Spheres(X,ra));
        }
        else
        {
            if (sph.Size()>0) sph.Last()->DelFrom (win);
            sph.Push (new VTK::Spheres(X,ra));
        }
        sph.Last()->Ids = ids.GetPtr();
        if (show_ids) sph.Last()->ShowIds  (0,0,0,0.003,10,false);
        else          sph.Last()->SetColor ("red",1.0);
        sph.Last()->AddTo (win, /*rstcam*/(stp_out==0));
        win.Show();
    }
    for (size_t i=0; i<sph.Size(); ++i) delete sph[i];

    // end
    return 0;
}
MECHSYS_CATCH
