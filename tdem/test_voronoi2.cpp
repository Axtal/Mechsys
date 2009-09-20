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
using std::ofstream;

const double x_min=-1,x_max=1;
const double y_min=-1,y_max=1;
const double z_min=-1,z_max=1;
const double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);

// Set up the number of blocks that the container is divided into
const int n_x=6,n_y=6,n_z=6;

// Set the number of particles that are going to be randomly introduced
const int particles=2;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}
int main(int argc, char **argv) try
{

        int i;
        double x,y,z;

        // Create a container with the geometry given above, and make it
        // non-periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block
        container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                        false,false,false,8);

        // Randomly add particles into the container
        for(i=0;i<particles;i++) {
                x=x_min+rnd()*(x_max-x_min);
                y=y_min+rnd()*(y_max-y_min);
                z=z_min+rnd()*(z_max-z_min);
                con.put(i,x,y,z);
        }
        Domain D;
        D.AddVoronoiContainer(con,0.1);

        std::ofstream  of2("test_voronoi2.bpy",std::ios::out);
        BPYHeader  (of2);
        D.WriteBPY (of2);
        of2.close  ();
        // Output the particle positions in gnuplot format
        con.draw_particles("random_points_p.gnu");

        // Output the Voronoi cells in gnuplot format
        con.draw_cells_gnuplot("random_points_v.gnu");
}
MECHSYS_CATCH
