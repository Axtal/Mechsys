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
#include <mechsys/mesh/unstructured.h>

using std::cout;
using std::endl;
using std::map;
using DEM::Domain;

int main(int argc, char **argv) try
{
    
    /////////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    Mesh::Unstructured mesh(/*NDim*/3);
    mesh.GenBox  (/*O2*/false,/*V*/0.05);
    mesh.WriteVTU ("test_tets");
    cout << " File <test_tets.vtu> generated\n";

    /////////////////////////////////////////////////////////////////////////////////////////// Domain /////
    
    Domain d;
    d.GenFromMesh (-1,mesh,/*R*/0.01,/*rho*/1.0);

    //////////////////////////////////////////////////////////////////////////////////// First timestep /////
    
    d.CamPos= 4.0,3.0,3.0;
    d.WritePOV ("test_tets");
    d.WriteBPY ("test_tets");

    return 0;    
}
MECHSYS_CATCH
