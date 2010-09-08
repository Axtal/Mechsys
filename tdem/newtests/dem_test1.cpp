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
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/maps.h>
#include <mechsys/dem/domain.h>

using std::cout;
using std::endl;
using Util::PI;

void Report (DEM::Domain & Dom, void * UserData)
{
    int * idx_out = static_cast<int*>(UserData);

    // header
    Array<String> keys("id", "xc", "yc", "zc", "ra", "vx", "vy", "vz", "ct");
    std::ostringstream oss;
    oss << Util::_6 << keys[0];
    for (size_t i=1; i<keys.Size()-1; ++i) { oss << Util::_8s << keys[i]; }
    oss << Util::_6 << keys.Last() << "\n";

    // values
    for (size_t i=0; i<Dom.Particles.Size(); ++i)
    {
        oss << Util::_6  << i;
        oss << Util::_8s << Dom.Particles[i]->x(0);
        oss << Util::_8s << Dom.Particles[i]->x(1);
        oss << Util::_8s << Dom.Particles[i]->x(2);
        oss << Util::_8s << Dom.Particles[i]->Props.R;
        oss << Util::_8s << Dom.Particles[i]->v(0);
        oss << Util::_8s << Dom.Particles[i]->v(1);
        oss << Util::_8s << Dom.Particles[i]->v(2);
        oss << Util::_6  << -1;
        oss << "\n";
    }

    // open file and save data
    String buf;   buf.Printf("dem_test1_%08d.res",(*idx_out));
    std::ofstream of(buf.CStr(), std::ios::out);
    of << oss.str();
    of.close();

    (*idx_out) += 1;
}

int main(int argc, char **argv) try
{
    // input
    double dt = 0.001;
    if (argc>1) dt = atof(argv[1]);

    // read data
    Table tab;
    tab.Read ("parts4.dat");
    Array<double> const & xc = tab("xc");
    Array<double> const & yc = tab("yc");
    Array<double> const & zc = tab("zc");
    Array<double> const & ra = tab("ra");
    Array<double> const & vx = tab("vx");
    Array<double> const & vy = tab("vy");
    Array<double> const & vz = tab("vz");

    // allocate particles
    int idx_out = 0;
    DEM::Domain dom(&idx_out);
    double Ekin0 = 0.0;
    double mass  = 1.0;
    for (size_t i=0; i<xc.Size(); ++i)
    {
        double vol = 4.0*PI*pow(ra[i],3.0)/3.0;
        double rho = mass/vol;
        dom.AddSphere (-1, Vec3_t(xc[i],yc[i],zc[i]), ra[i], rho);
        Particle & p = (*dom.Particles.Last());
        p.v = vx[i], vy[i], vz[i];
        p.Initialize();
        Ekin0 += 0.5*p.Props.m*dot(p.v,p.v);
    }

    // properties
    Dict prps;
    prps.Set (-1,"Kn Kt Gn Gt Mu Beta Eta", 1000.0, 0., 0., 0., 0. ,0.,0.);
    dom.SetProps (prps);

    // solve
    double tf    = 5.0;
    double dtout = 0.1;
    dom.Solve (tf, dt, dtout, NULL, &Report);

    // energy
    double Ekin1 = 0.0;
    for (size_t i=0; i<dom.Particles.Size(); ++i)
    {
        Particle const & p = (*dom.Particles[i]);
        Ekin1 += 0.5*p.Props.m*dot(p.v,p.v);
    }
    printf("\nEkin (before) = %16.8e\n", Ekin0);
    printf("Ekin (after)  = %16.8e\n\n", Ekin1);

    // write control file
    std::ofstream of("dem_test1_control.res", std::ios::out);
    of << "fkey  " << "dem_test1" << "\n";
    of << "nout  " << idx_out     << "\n";
    of << "nx    " << 1           << "\n";
    of << "ny    " << 1           << "\n";
    of << "nz    " << 1           << "\n";
    of << "lxmi  " << -2.0        << "\n";
    of << "lxma  " <<  2.0        << "\n";
    of << "lymi  " << -2.0        << "\n";
    of << "lyma  " <<  2.0        << "\n";
    of << "lzmi  " <<  0.0        << "\n";
    of << "lzma  " <<  0.1        << "\n";
    of << "dt    " <<  dt         << "\n";
    of.close();

    // end
    return 0;
}
MECHSYS_CATCH
