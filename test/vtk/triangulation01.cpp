/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
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

// Std Lib
#include <iostream>

// VTK
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkShrinkFilter.h>
#include <vtkDataSetMapper.h>
#include <vtkOutlineFilter.h>
#include <vtkProperty.h>
#include <vtkPolyDataMapper.h>

// MechSys
#include <mechsys/vtk/win.h>
#include <mechsys/vtk/sgrid.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    size_t ndiv = 6;
    if (argc>1) ndiv = atoi(argv[1]);

    // Create a meshgrid
    VTK::SGrid gri(ndiv,0.0,1.0,  // X
                   ndiv,0.0,1.0,  // Y
                   ndiv,0.0,1.0); // Z

    // Create a 3D triangulation
    //   - The Tolerance is the distance that nearly coincident points are merged together.
    //   - Delaunay does better if points are well spaced.
    //   - The alpha value is the radius of circumcircles, circumspheres.
    //     Any mesh entity whose circumcircle is smaller than this value is output.
    vtkPolyData   * vertices = vtkPolyData::New();
    vtkDelaunay3D * delaunay = vtkDelaunay3D::New();
    vertices -> SetPoints(gri.GetPoints());
    delaunay -> SetInput(vertices);
    delaunay -> SetTolerance(0.01);
    delaunay -> SetAlpha(0.2);
    delaunay -> BoundingTriangulationOff();

    // Shrink the result to help see it better.
    vtkShrinkFilter * shrink = vtkShrinkFilter::New();
    shrink -> SetInputConnection(delaunay->GetOutputPort());
    shrink -> SetShrinkFactor(0.9);

    // Mapper and actor
    vtkDataSetMapper * mapper = vtkDataSetMapper::New();
    vtkActor         * actor  = vtkActor::New();
    mapper -> SetInputConnection (shrink->GetOutputPort());
    actor  -> SetMapper          (mapper);
    actor  -> GetProperty        () -> SetColor(0,1,0);

    // Create a box around the points (Outline)
    vtkOutlineFilter  * outline        = vtkOutlineFilter::New();
    vtkPolyDataMapper * outline_mapper = vtkPolyDataMapper::New();
    vtkActor          * outline_actor  = vtkActor::New();
    outline        -> SetInputConnection(shrink->GetOutputPort());
    outline_mapper -> SetInput    (outline->GetOutput());
    outline_actor  -> SetMapper   (outline_mapper);
    outline_actor  -> GetProperty () -> SetColor(0,0,0);

    // Window
    VTK::Win win;
    win.AddActor (actor);
    win.AddActor (outline_actor);
    win.Show     ();

    // Clean up
    vertices       -> Delete();
    delaunay       -> Delete();
    shrink         -> Delete();
    mapper         -> Delete();
    actor          -> Delete();
    outline        -> Delete();
    outline_mapper -> Delete();
    outline_actor  -> Delete();

    return 0;
}
MECHSYS_CATCH
