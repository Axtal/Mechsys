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

#ifndef MECHSYS_CYLINDER_H
#define MECHSYS_CYLINDER_H

// VTK
#include <vtkGlyph3D.h>
#include <vtkCylinderSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkLODActor.h>
#include <vtkProperty.h>

// MechSys
#include <mechsys/vtk/vtkwin.h>
#include <mechsys/util/colors.h>
#include <mechsys/linalg/matvec.h>

class Cylinder
{
public:
    // Constructor & Destructor
     Cylinder (Vec3_t const & X0, Vec3_t const & X1, double Radius, bool Capping=true, int Resolution=20);
    ~Cylinder ();

    // Methods
    void AddActorsTo (VTKWin & Win) { Win.AddActor(_cylin_actor); }

    // Set methods
    Cylinder & SetColor  (char const * Name="peacock", double Opacity=0.8);
    void       SetPoints (Vec3_t const & X0, Vec3_t const & X1);

private:
    vtkGlyph3D        * _cylin;
    vtkPolyDataMapper * _cylin_mapper;
    vtkLODActor       * _cylin_actor;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Cylinder::Cylinder (Vec3_t const & X0, Vec3_t const & X1, double Radius, bool Capping, int Resolution)
{
    // cylinder length
    Vec3_t V(X1-X0);

    // create cylinder
    double len = Norm(V);
    double cen = len/2.0;
    vtkCylinderSource * cylin = vtkCylinderSource ::New(); // will be parallel to y-axis
    cylin -> SetCenter     (0.0, cen, 0.0);
    cylin -> SetRadius     (Radius);
    cylin -> SetHeight     (len);
    cylin -> SetCapping    (Capping);
    cylin -> SetResolution (Resolution);

    // rotate around z axis because glyph3D needs the cylinder to be parallel to the x-axis
    vtkTransform       * rotate        = vtkTransform       ::New();
    vtkTransformFilter * cylin_along_x = vtkTransformFilter ::New();
    rotate        -> RotateZ      (-90.0);
    cylin_along_x -> SetTransform (rotate);
    cylin_along_x -> SetInput     (cylin->GetOutput());

    // setup data for the glyph
    vtkPoints     * points   = vtkPoints     ::New();
    vtkFloatArray * vectors  = vtkFloatArray ::New();
    vtkPolyData   * polydata = vtkPolyData   ::New();
    points   -> SetNumberOfPoints     (1);
    points   -> InsertPoint           (0, X0(0), X0(1), X0(2));
    vectors  -> SetNumberOfComponents (3);
    vectors  -> SetNumberOfTuples     (1);
    vectors  -> InsertTuple3          (0, V(0), V(1), V(2));
    polydata -> SetPoints             (points);
    polydata -> GetPointData() -> SetVectors(vectors);

    // create the glyph
    _cylin =  vtkGlyph3D                   ::New();
    _cylin -> SetInput                     (polydata);
    _cylin -> SetSource                    (cylin_along_x->GetPolyDataOutput());
    _cylin -> SetVectorModeToUseVector     ();
    _cylin -> SetScaleModeToDataScalingOff ();

    // create mapper and actor
    _cylin_mapper = vtkPolyDataMapper   ::New();
    _cylin_actor  = vtkLODActor         ::New();
    _cylin_mapper -> SetInputConnection (_cylin->GetOutputPort());
    _cylin_actor  -> SetMapper          (_cylin_mapper);
    _cylin_actor  -> VisibilityOn       ();
    SetColor ();

    // clean up
    cylin         -> Delete();
    rotate        -> Delete();
    cylin_along_x -> Delete();
    points        -> Delete();
    vectors       -> Delete();
    polydata      -> Delete();
}

inline Cylinder::~Cylinder ()
{
    _cylin        -> Delete();
    _cylin_mapper -> Delete();
    _cylin_actor  -> Delete();
}

inline Cylinder & Cylinder::SetColor (char const * Name, double Opacity)
{
    Vec3_t c(Colors::Get(Name));
    _cylin_actor->GetProperty()->SetColor   (c(0), c(1), c(2));
    _cylin_actor->GetProperty()->SetOpacity (Opacity);
    return (*this);
}

inline void Cylinder::SetPoints (Vec3_t const & X0, Vec3_t const & X1)
{
    Vec3_t V(X1-X0);
    vtkPoints     * points   = vtkPoints     ::New();
    vtkFloatArray * vectors  = vtkFloatArray ::New();
    vtkPolyData   * polydata = static_cast<vtkPolyData*>(_cylin->GetInput());
    points   -> SetNumberOfPoints     (1);
    points   -> InsertPoint           (0, X0(0), X0(1), X0(2));
    vectors  -> SetNumberOfComponents (3);
    vectors  -> SetNumberOfTuples     (1);
    vectors  -> InsertTuple3          (0, V(0), V(1), V(2));
    polydata -> SetPoints             (points);
    polydata -> GetPointData() -> SetVectors(vectors);
    points   -> Delete();
    vectors  -> Delete();
}

#endif // MECHSYS_CYLINDER_H
