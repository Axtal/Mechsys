/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
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

// Std lib
#include <math.h>

// GSL
#include <gsl/gsl_linalg.h>

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/dem/distance.h>
#include <mechsys/util/fatal.h>
#include <mechsys/mesh/structured.h>

using std::cout;
using std::endl;
using std::map;
using DEM::Domain;

int main(int argc, char **argv) try
{
    
    /////////////////////////////////////////////////////////////////////////////////////////// Domain /////
    
    Domain d;

    //////////////////////////////////////////////////////////////////////////////////// Figures to draw ///
    
    //d.AddVoroPack (-1, 0.1, 2,2,2, 2,2,5, 1.0, true, 1000, 1.0, 0.8);
    d.GenRice(-1,1,1,0.2,1.0,80000,1.0);
    d.CamPos= 14.0,13.0,13.0;
    d.WriteBPY ("figure");

    /////////////////////////////////////////////////////////////////////////////////////////// Domain /////
    
    Domain e;

    //////////////////////////////////////////////////////////////////////////////////// Figures to draw ///

    e.AddCube(1,Vec3_t(-1.0,0.0,0.0),0.1,1.8,1.0,0.0,&OrthoSys::e0);
    e.AddCube(2,Vec3_t( 1.0,0.0,0.0),0.1,1.8,1.0,0.0,&OrthoSys::e0);
    e.Particles[0]->FixVeloc();
    e.Particles[1]->FixVeloc();
    e.AddRice(3,OrthoSys::O,0.2,1.8,1.0,M_PI/2.0,&OrthoSys::e1);
    e.CamPos =  1.0,5.0,2.0;
    e.WritePOV("figure2");


    return 0;    
}
MECHSYS_CATCH
