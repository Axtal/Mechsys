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
#include "dem/domain.h"
#include "util/fatal.h"

using std::cout;
using std::endl;

// Random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main(int argc, char **argv) try
{
    // constants
    const double x_min=-1, x_max=1;
    const double y_min=-1, y_max=1;
    const double z_min=-1, z_max=1;
    const double cvol = (x_max-x_min)*(y_max-y_min)*(x_max-x_min);

    // number of blocks that the container is divided into
    const int n_x=6, n_y=6, n_z=6;

    // create voro container.allocate space for eight particles within each computational block
    size_t num_particles = 2;
    container con (x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z, false,false,false,8);

    // randomly add particles into the container
    for (size_t i=0; i<num_particles; i++)
    {
        double x = x_min+rnd()*(x_max-x_min);
        double y = y_min+rnd()*(y_max-y_min);
        double z = z_min+rnd()*(z_max-z_min);
        con.put (i,x,y,z);
    }

    // domain
    Domain d;
    d.GenFromVoro (con, /*R*/0.1);

    // output
    d.WriteBPY ("test_voro02");

    // output the particle positions in gnuplot format
    con.draw_particles("test_voro02_p.gnu");

    // Output the Voronoi cells in gnuplot format
    con.draw_cells_gnuplot("test_voro02_v.gnu");
}
MECHSYS_CATCH
