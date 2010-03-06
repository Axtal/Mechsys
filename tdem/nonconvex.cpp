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

using std::cout;
using std::endl;
using std::map;
using std::pair;
using DEM::Domain;

int main(int argc, char **argv) try
{
    // camera coordinates
    double cx=1., cy=1., cz=1.;
    if (argc>1) cx = atof(argv[1]);
    if (argc>2) cy = atof(argv[2]);
    if (argc>3) cz = atof(argv[3]);

    // nonconvex particle
    Mesh::Generic mesh(/*NDim*/3);
    mesh.ReadMesh("nonconvex", /*IsShell*/true); // read file nonconvex.mesh
    Particle * p0 = new Particle(-1, mesh, 0.1);

    // domain
	Domain dom;
    dom.Particles.Push (p0);
    dom.FreeParticles = dom.Particles;
    dom.CamPos= cx,cy,cz;
    dom.WritePOV ("nonconvex");
    dom.WriteBPY ("nonconvex");
    return 0;
}
MECHSYS_CATCH
