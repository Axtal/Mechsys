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

double rnd() {return double(rand())/RAND_MAX;}
int main(int argc, char **argv) try
{

	double x,y,z,rsq,r;
	voronoicell v,u;
	v.init(-1,1,-1,1,-1,1);
	for(int i=0;i<20;i++) 
    {
		x=2*rnd()-1;
		y=2*rnd()-1;
		z=2*rnd()-1;
		rsq=x*x+y*y+z*z;
		if(rsq>0.01&&rsq<1) 
        {
			r=1/sqrt(rsq);x*=r;y*=r;z*=r;
			v.plane(x,y,z,1);
		}
	}

	u.init(-1,1,-1,1,-1,1);
	for(int i=0;i<20;i++) 
    {
		x=2*rnd()-1;
		y=2*rnd()-1;
		z=2*rnd()-1;
		rsq=x*x+y*y+z*z;
		if(rsq>0.01&&rsq<1) 
        {
			r=1/sqrt(rsq);x*=r;y*=r;z*=r;
			u.plane(x,y,z,1);
		}
	}
    

    Domain D;
    D.AddVoronoiCell(v,0.1);
    D.AddVoronoiCell(u,0.1);
    
    Vec3_t trans(2.5,0,0);
    D.Particles[0]->Translation(trans);
    trans = -trans;
    D.Particles[1]->Translation(trans);
    D.Particles[0]->CalcMassProperties();
    D.Particles[1]->CalcMassProperties();
    D.Particles[0]->v = Vec3_t(-0.5,0,0);

    double dt =0.001;
    D.Initialize(dt);
    D.Solve(0,30,dt,0.1,"test_voronoi");


    std::ofstream of2("test_voronoi.bpy",std::ios::out);
    BPYHeader  (of2);
    D.WriteBPY (of2);
    of2.close  ();
	// Output the Voronoi cell to a file, in the gnuplot format
	//v.draw_gnuplot("single_cell.gnu",0,0,0);
    //cout << v.pts[3][1] << " " << v.pts[4][1] << " " << v.pts[5][1] <<endl;
}
MECHSYS_CATCH
