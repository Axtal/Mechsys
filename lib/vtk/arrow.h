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

#ifndef MECHSYS_ARROW3D_H
#define MECHSYS_ARROW3D_H

// VTK
#include <vtkConeSource.h>
#include <vtkCylinderSource.h>
#include <vtkAppendPolyData.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkLODActor.h>

// MechSys
#include <mechsys/vtk/vtkwin.h>
#include <mechsys/util/colors.h>
#include <mechsys/linalg/matvec.h>

class Arrow
{
public:
    /*                +-------------------------------+ 
     *                |            length             |
     *                +-----------------------+-------+
     *                |        bod_len        |tip_len|
     *                |                       |       |
     *                                        `.      ----+  
     *                                        | ``.       |
     *             +  +-----------------------|    ``.    |
     *     bod_rad |  |           +           |   +   >   | tip_rad   
     *             +  +-----------|-----------|   |_-'    |
     *                |           |           | _-|       |
     *                |           |           ''  |     --+  
     *                |           |               |
     *                +-----------+---------------+-------> y axis
     *                |           |               |    
     *                y0      y_bod_cen      y_tip_cen
     */

    // Constructor & Destructor
     Arrow (Vec3_t const & X, Vec3_t const & V, double BodRad=0.015, double TipRad=0.03, double TipLen=0.12, int Resolution=18);
    ~Arrow ();

    // Set methods
    Arrow & SetColor (char const * Name="yellow", double Opacity=1.0);

    // Methods
    void AddActorsTo (VTKWin & Win) { Win.AddActor(_vector_actor); }

    // Data
    double BodRad;     // body (cylinder) radius
    double TipRad;     // tip (cone) radius
    double TipLen;     // tip (cone) length/height
    int    Resolution; // number of facets

private:
    // Data
    vtkGlyph3D        * _vector;
    vtkPolyDataMapper * _vector_mapper;
    vtkLODActor       * _vector_actor;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Arrow::Arrow (Vec3_t const & X, Vec3_t const & V, double BR, double TR, double TL, int Res)
    : BodRad(BR), TipRad(TR), TipLen(TL), Resolution(Res)
{
    // Dimensions
    double length  = Norm(V);
    double bod_len = length-TipLen;
    double bod_cen = bod_len/2.0;
    double tip_cen = bod_len+TipLen/2.0;

    // Create the body (cylinder) and tip (cone) of the arrow
    vtkConeSource     * cone  = vtkConeSource     ::New();
    vtkCylinderSource * cylin = vtkCylinderSource ::New();
    vtkAppendPolyData * arrow = vtkAppendPolyData ::New();
    cylin -> SetCenter     (0.0, bod_cen, 0.0);
    cylin -> SetHeight     (bod_len);
    cylin -> SetRadius     (BodRad);
    cylin -> SetResolution (Resolution);
    cone  -> SetCenter     (0.0, tip_cen, 0.0);
    cone  -> SetDirection  (0,1,0); // because cylinder direction is along y-axis
    cone  -> SetHeight     (TipLen);
    cone  -> SetRadius     (TipRad);
    cone  -> SetResolution (Resolution);
    arrow -> AddInput      (cone->GetOutput());
    arrow -> AddInput      (cylin->GetOutput());

    // Rotate arrow around z axis because glyph3D creates arrows parallel to the x-axis
    vtkTransform       * rotate        = vtkTransform       ::New();
    vtkTransformFilter * arrow_along_x = vtkTransformFilter ::New();
    rotate        -> RotateZ      (-90.0);
    arrow_along_x -> SetTransform (rotate);
    arrow_along_x -> SetInput     (arrow->GetOutput());

    // Setup data for the glyph
    vtkPoints     * points   = vtkPoints     ::New();
    vtkFloatArray * vectors  = vtkFloatArray ::New();
    vtkPolyData   * polydata = vtkPolyData   ::New();
    points   -> SetNumberOfPoints     (1);
    points   -> InsertPoint           (0, X(0), X(1), X(2));
    vectors  -> SetNumberOfComponents (3);
    vectors  -> SetNumberOfTuples     (1);
    vectors  -> InsertTuple3          (0, V(0), V(1), V(2));
    polydata -> SetPoints             (points);
    polydata -> GetPointData() -> SetVectors(vectors);

    // Create the glyph
    _vector =  vtkGlyph3D                   ::New();
    _vector -> SetInput                     (polydata);
    _vector -> SetSource                    (arrow_along_x->GetPolyDataOutput());
    _vector -> SetVectorModeToUseVector     ();
    _vector -> SetScaleModeToDataScalingOff ();

    // Create mapper and actor
    _vector_mapper = vtkPolyDataMapper   ::New();
    _vector_actor  = vtkLODActor         ::New();
    _vector_mapper -> SetInputConnection (_vector->GetOutputPort());
    _vector_actor  -> SetMapper          (_vector_mapper);
    _vector_actor  -> VisibilityOn       ();
    SetColor ();

    // Clean up
    cone          -> Delete();
    cylin         -> Delete();
    arrow         -> Delete();
    rotate        -> Delete();
    arrow_along_x -> Delete();
    points        -> Delete();
    vectors       -> Delete();
    polydata      -> Delete();
}

inline Arrow::~Arrow ()
{
    _vector        -> Delete();
    _vector_mapper -> Delete();
    _vector_actor  -> Delete();
}

inline Arrow & Arrow::SetColor (char const * Name, double Opacity)
{
    Vec3_t c(Colors::Get(Name));
    _vector_actor->GetProperty()->SetColor   (c(0), c(1), c(2));
    _vector_actor->GetProperty()->SetOpacity (Opacity);
    return (*this);
}

#endif // MECHSYS_ARROW3D_H
