/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using DEM::Domain;

// Random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main(int argc, char **argv) try
{
    // constants
    const double x_min=-1, x_max=1;
    const double y_min=-1, y_max=1;
    const double z_min=-1, z_max=1;

    // number of blocks that the container is divided into
    const size_t n_x=3, n_y=3, n_z=3;

    // create voro container.allocate space for eight particles within each computational block
    container con1 (x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z, true,true,true,8);
    container con2 (x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z, true,true,true,8);

    // randomly add particles into the container
    size_t n = 0;
    for (size_t i=0; i<n_x; i++)
    {
        double x = x_min+(i+rnd())*(x_max-x_min)/n_x;
        for (size_t j=0; j<n_y; j++)
        {
            double y = y_min+(j+rnd())*(y_max-y_min)/n_y;
            for (size_t k=0; k<n_z; k++)
            {
                double z = z_min+(k+rnd())*(z_max-z_min)/n_z;
                con1.put (n,x,y,z);
                n++;
            }
        }
    }

    // randomly add particles into the container
    n = 0;
    for (size_t i=0; i<n_x; i++)
    {
        double x = x_min+(i+rnd())*(x_max-x_min)/n_x;
        for (size_t j=0; j<n_y; j++)
        {
            double y = y_min+(j+rnd())*(y_max-y_min)/n_y;
            for (size_t k=0; k<n_z; k++)
            {
                double z = z_min+(k+rnd())*(z_max-z_min)/n_z;
                con2.put (n,x,y,z);
                n++;
            }
        }
    }

    // domain
    Domain d;
    d.GenFromVoro (-1,con1, /*R*/0.1,3.0);
    for (size_t i=0; i < d.Particles.Size() ; ++i)
    {
        Vec3_t trans(3,0,0);
        d.Particles[i]->Translate(trans);
        d.Particles[i]->v = Vec3_t(-0.6,0,0);
    }

    d.GenFromVoro(-1,con2,0.1,3.0);

    // output
    d.WriteBPY ("test_voro02");

    // output the particle positions in gnuplot format
    //con.draw_particles("test_voro02_p.gnu");

    // Output the Voronoi cells in gnuplot format
    //con.draw_cells_gnuplot("test_voro02_v.gnu");
    d.CamPos = 0.0,5.0,3.0;
    d.Solve (/*tf*/30, 0.001, /*dtOut*/0.1, "test_voro02", true);
}
MECHSYS_CATCH
