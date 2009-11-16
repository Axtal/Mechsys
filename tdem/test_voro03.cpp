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
using DEM::Domain;

// Random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main(int argc, char **argv) try
{
    // constants
    const double x_min=-1, x_max=1;
    const double y_min=-1, y_max=1;
    const double z_min=-1, z_max=1;
    //const double cvol = (x_max-x_min)*(y_max-y_min)*(x_max-x_min);

    // number of blocks that the container is divided into
    const int n_x=2, n_y=2, n_z=2;

    // create voro container.allocate space for eight particles within each computational block
    //size_t num_particles = n_x*n_y*n_z;
    double qinter = 0.0,radius = 0.05;
    //cin >> qinter;
    //cin >> radius;
    container con1 (x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z, false,false,false,8);

    // randomly add particles into the container
    size_t n = 0;
    for (int i=0; i<n_x; i++)
    {
        double x = x_min+(i+qinter+(1-2*qinter)*rnd())*(x_max-x_min)/n_x;
        for (int j=0; j<n_y; j++)
        {
            double y = y_min+(j+qinter+(1-2*qinter)*rnd())*(y_max-y_min)/n_y;
            for (int k=0; k<n_z; k++)
            {
                double z = z_min+(k+qinter+(1-2*qinter)*rnd())*(z_max-z_min)/n_z;
                con1.put (n,x,y,z);
                n++;
            }
        }
    }

    // domain
    Domain d;
    d.GenFromVoro (-1,con1, /*R*/radius,1.0);
    for (size_t i=0; i < d.Particles.Size() ; ++i)
    {
        Vec3_t trans(x_min,y_min,z_min+2.0);
        d.Particles[i]->Translate(trans);
        d.Particles[i]->v = Vec3_t(0,0,0);
    }

    Vec3_t r(0.0,0.0,-1.0);
    // TODO this is not compiling: d.AddPlane(-2,r,0.1,10.0,1.0);

    r = 0.0,0.0,5.0;
    // TkODO this is not compiling: d.AddPlane(-3,r,0.1,10.0,1.0);
    // Dictionary of parameters
    Dict B;
    B.Set(-2,"vx vy vz",0.0,0.0,0.0);
    B.Set(-3,"fx fy fz",0.0,0.0,-1.0);


    // output
    d.WriteBPY ("test_voro03");

    // output the particle positions in gnuplot format
    //con.draw_particles("test_voro02_p.gnu");

    // Output the Voronoi cells in gnuplot format
    //con.draw_cells_gnuplot("test_voro02_v.gnu");
    d.SetBC(B);
    d.CamPos = 0.0,15.0,0.0;
    d.Solve (/*tf*/30, 0.001, /*dtOut*/0.1, "test_voro03", true);
}
MECHSYS_CATCH
