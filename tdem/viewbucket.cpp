/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
 * Copyright (C) 2013 William Oquendo                                   *
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

// VTK
#include <vtkLightKit.h>

// C++ STL
#include <fstream>

// MechSys
#include <mechsys/dem/graph.h>
#include <mechsys/dem/domain.h>
#include <mechsys/dem/visualise.h>
#include <mechsys/vtk/win.h>
#include <mechsys/vtk/axes.h>

int main(int argc, char **argv) try
{
    // Declare the DEM domain
    DEM::Domain d;

    // Introduce the mesh of a particle
    double R = 0.5;
    Vec3_t ex (-10,0.0,0.0);
    d.Particles.Push(new DEM::Particle(-2,"Rawbucket",R,3.0,10.0));// Changed tag to -2

    std::ofstream of1("viewbucket.py", std::ios::out);
    DEM::BPYHeader(of1);
    d.Particles[0]->Draw(of1, "Red", true);
    of1.close();

    std::ofstream of2("viewbucket.pov", std::ios::out);
    DEM::POVHeader(of2);
    d.Particles[0]->Draw(of2, "Red", false);
    of2.close();

    bool showvert = true;
    bool showedge = true;
    DEM::Visualise vis(d, /*parts*/Array<int>(-2,true), /*walls*/Array<int>(-1,true), showvert, showedge);
	VTK::Axes axe(1.0, /*hline*/false, /*reverse*/false, /*full*/true);

	VTK::Win win;

	/*
	vis.PartFaces.SetMaterial (0.2, 0.1, 1.0, 200.0);
	vis.PartFaces.SetColor    ("cobalt", 1.0);
	double elevation = 60.0;
	double azimuth   = -45.0;./o./
	vtkLightKit *lgt = NULL;
	lgt->AddLightsToRenderer (win.GetRen());
	lgt->SetKeyLightAngle (elevation, azimuth);
	lgt->SetKeyToHeadRatio (0.01);
	lgt->SetKeyToFillRatio (1.2);
	*/

	if (false)
	{
		vis.AddTo(win);
		axe.AddTo(win);
		win.Render();
		win.Show();
	}

    return 0;
}
MECHSYS_CATCH
