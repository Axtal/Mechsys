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

#ifndef MECHSYS_PLANE_H
#define MECHSYS_PLANE_H

// Std Lib
#include <string>

// VTK
#include <vtkPlaneSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

// MechSys
#include <mechsys/vtk/vtkwin.h>
#include <mechsys/util/colors.h>
#include <mechsys/linalg/matvec.h>

class Plane
{
public:
    // Constructors
    Plane (Vec3_t const & Cen, Vec3_t const & n);                                         // centre, and normal
    Plane (Vec3_t const & Ori, Vec3_t const & Pt1, Vec3_t const & Pt2, Vec3_t const & n); // origin, point1, point2, normal

    // Destructor
    ~Plane ();

    // Methods
    void AddActorsTo (VTKWin & Win) { Win.AddActor(_plane_actor); }

    // Set methods
    Plane & SetColor (char const * Name="light_salmon", double Opacity=0.5);

private:
    vtkPlaneSource    * _plane;
    vtkPolyDataMapper * _plane_mapper;
    vtkActor          * _plane_actor;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Plane::Plane (Vec3_t const & Cen, Vec3_t const & n)
{
    _plane        = vtkPlaneSource    ::New();
    _plane_mapper = vtkPolyDataMapper ::New();
    _plane_actor  = vtkActor          ::New();
    _plane        -> SetCenter          (Cen(0), Cen(1), Cen(2));
    _plane        -> SetNormal          (  n(0),   n(1),   n(2));
    _plane_mapper -> SetInputConnection (_plane->GetOutputPort());
    _plane_actor  -> SetMapper          (_plane_mapper);
    SetColor ();
}

inline Plane::Plane (Vec3_t const & Ori, Vec3_t const & Pt1, Vec3_t const & Pt2, Vec3_t const & n)
{
    _plane        = vtkPlaneSource    ::New();
    _plane_mapper = vtkPolyDataMapper ::New();
    _plane_actor  = vtkActor          ::New();
    _plane        -> SetOrigin          (Ori(0), Ori(1), Ori(2));
    _plane        -> SetNormal          (  n(0),   n(1),   n(2));
    _plane        -> SetPoint1          (Pt1(0), Pt1(1), Pt1(2));
    _plane        -> SetPoint2          (Pt2(0), Pt2(1), Pt2(2));
    _plane_mapper -> SetInputConnection (_plane->GetOutputPort());
    _plane_actor  -> SetMapper          (_plane_mapper);
    SetColor ();
}

inline Plane::~Plane ()
{
    _plane        -> Delete();
    _plane_mapper -> Delete();
    _plane_actor  -> Delete();
}

inline Plane & Plane::SetColor (char const * Name, double Opacity)
{
    Vec3_t c(Colors::Get(Name));
    _plane_actor->GetProperty()->SetColor   (c(0), c(1), c(2));
    _plane_actor->GetProperty()->SetOpacity (Opacity);
    return (*this);
}

#endif // MECHSYS_PLANE_H
