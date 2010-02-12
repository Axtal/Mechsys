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
    
    //Domain d;

    //////////////////////////////////////////////////////////////////////////////////// Figures to draw ///
    
    //d.AddVoroPack (-1, 0.1, 10,10,10, 10,10,10, 1.0, true, 1000, 1.0);
    //d.FreeParticles = d.Particles; //Neccesary for the particles to be drawn
    //d.CamPos= 14.0,13.0,13.0;
    //d.WritePOV ("figure");

    /////////////////////////////////////////////////////////////////////////////////////////// Domain /////
    
    Domain e;

    //////////////////////////////////////////////////////////////////////////////////// Figures to draw ///

    Vec3_t axis(2,0,1);
    e.AddCube(1,Vec3_t(-5.0,0.0,0.0),0.01,2.0,1.0,30,&axis);
    e.AddCube(1,Vec3_t(-0.0,0.0,0.0),0.01,1.4,1.0,30,&axis);
    e.AddCube(1,Vec3_t( 5.0,0.0,0.0),0.30,1.4,1.0,30,&axis);
    e.CamPos =  0.0,15.0,0.0;
    e.FreeParticles = e.Particles; //Neccesary for the particles to be drawn
    e.WritePOV("figure2");


    return 0;    
}
MECHSYS_CATCH
