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

#ifndef MECHSYS_TRIANGULATE_H
#define MECHSYS_TRIANGULATE_H

// Std Lib
#include <vector>

// VTK
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkUnstructuredGridWriter.h>

// MechSys
#include <mechsys/vtk/meshgrid.h>

class Triangulate
{
public:
    // Constructors
    Triangulate (double const * X,
                 double const * Y,
                 double const * Z, int Size, char const * Filename);

    Triangulate (std::vector<double> & X,
                 std::vector<double> & Y,
                 std::vector<double> & Z, char const * Filename);

    // Methods
    void Start ();

private:
    bool                  _use_vector;
    int                   _npoints;
    char const          * _filename;
    double const        * _X; 
    double const        * _Y; 
    double const        * _Z; 
    std::vector<double>   _vX;
    std::vector<double>   _vY;
    std::vector<double>   _vZ;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Triangulate::Triangulate (double const * X,
                                 double const * Y,
                                 double const * Z, int Size, char const * Filename)
    : _use_vector(false), _npoints(Size), _filename(Filename), _X(X), _Y(Y), _Z(Z)
{
    if (Size<1) throw new Fatal("Triangulate::Triangulate: Size of arrays must be greater than zero");
}

inline Triangulate::Triangulate (std::vector<double> & X,
                                 std::vector<double> & Y,
                                 std::vector<double> & Z, char const * Filename)
    : _use_vector(true), _npoints(X.size()), _filename(Filename), _X(NULL), _Y(NULL), _Z(NULL), _vX(X), _vY(Y), _vZ(Z)
{
    if (X.size()<1)         throw new Fatal("Triangulate::Triangulate: Size of arrays must be greater than zero");
    if (X.size()!=Y.size()) throw new Fatal("Triangulate::Triangulate: X, Y, and Z arrays must have the same size");
    if (X.size()!=Z.size()) throw new Fatal("Triangulate::Triangulate: X, Y, and Z arrays must have the same size");
}

inline void Triangulate::Start ()
{
    // Create VTK points
    vtkPoints * points = vtkPoints::New();
    points->Allocate(_npoints);
    if (_use_vector)
    {
        for (int i=0; i<_npoints; ++i)
        {
            double P[3] = {_vX[i], _vY[i], _vZ[i]};
            points->InsertPoint(i,P);
        }
    }
    else 
    {
        for (int i=0; i<_npoints; ++i)
        {
            double P[3] = {_X[i], _Y[i], _Z[i]};
            points->InsertPoint(i,P);
        }
    }

    // Create a 3D triangulation
    //   - The Tolerance is the distance that nearly coincident points are merged together.
    //   - Delaunay does better if points are well spaced.
    //   - The alpha value is the radius of circumcircles, circumspheres.
    //     Any mesh entity whose circumcircle is smaller than this value is output.
    vtkPolyData   * vertices = vtkPolyData::New();
    vtkDelaunay3D * delaunay = vtkDelaunay3D::New();
    vertices -> SetPoints(points);
    delaunay -> SetInput(vertices);
    delaunay -> SetTolerance(0.01);
    delaunay -> SetAlpha(0.2);
    delaunay -> BoundingTriangulationOff();

    // Write file
    vtkUnstructuredGridWriter * writer = vtkUnstructuredGridWriter::New();
    writer -> SetInputConnection(delaunay->GetOutputPort());
    writer -> SetFileName(_filename);
    writer -> Write();

    // Clean up
    points   -> Delete();
    vertices -> Delete();
    delaunay -> Delete();
    writer   -> Delete();
}

#endif // MECHSYS_TRIANGULATE_H
