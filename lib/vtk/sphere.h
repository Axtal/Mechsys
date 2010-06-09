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

#ifndef MECHSYS_SPHERE_H
#define MECHSYS_SPHERE_H

// VTK
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

// MechSys
#include <mechsys/vtk/vtkwin.h>
#include <mechsys/util/colors.h>
#include <mechsys/linalg/matvec.h>

class Sphere
{
public:
    // Constructor & Destructor
     Sphere (Vec3_t const & x, double R, int ThetaRes=20, int PhiRes=20);
    ~Sphere ();

    // Methods
    void AddActorsTo (VTKWin & Win) { Win.AddActor(_sphere_actor); }
    void SetCenter   (Vec3_t const & x) { _sphere->SetCenter (x(0), x(1), x(2)); }

    // Set methods
    Sphere & SetColor (char const * Name="brown", double Opacity=0.8);

private:
    vtkSphereSource   * _sphere;
    vtkPolyDataMapper * _sphere_mapper;
    vtkActor          * _sphere_actor;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Sphere::Sphere (Vec3_t const & x, double R, int ThetaRes, int PhiRes)
{
    _sphere        = vtkSphereSource   ::New();
    _sphere_mapper = vtkPolyDataMapper ::New();
    _sphere_actor  = vtkActor          ::New();
    _sphere        -> SetCenter          (x(0), x(1), x(2));
    _sphere        -> SetRadius          (R);
    _sphere        -> SetThetaResolution (ThetaRes);
    _sphere        -> SetPhiResolution   (PhiRes);
    _sphere_mapper -> SetInputConnection (_sphere->GetOutputPort());
    _sphere_actor  -> SetMapper          (_sphere_mapper);
    SetColor ();
}

inline Sphere::~Sphere ()
{
    _sphere        -> Delete();
    _sphere_mapper -> Delete();
    _sphere_actor  -> Delete();
}

inline Sphere & Sphere::SetColor (char const * Name, double Opacity)
{
    Vec3_t c(Colors::Get(Name));
    _sphere_actor->GetProperty()->SetColor   (c(0), c(1), c(2));
    _sphere_actor->GetProperty()->SetOpacity (Opacity);
    return (*this);
}

#endif // MECHSYS_SPHERE_H
