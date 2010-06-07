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

#ifndef MECHSYS_CUBE_H
#define MECHSYS_CUBE_H

// VTK
#include <vtkCubeSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

// MechSys
#include <mechsys/vtk/vtkwin.h>
#include <mechsys/util/colors.h>

class Cube
{
public:
    // Constructor & Destructor
     Cube (Vec3_t const & Cen, double Lx=1.0, double Ly=1.0, double Lz=1.0);
    ~Cube ();

    // Set methods
    Cube & SetColor     (char const * Name="yellow", double Opacity=1.0);
    Cube & SetWireColor (char const * Name="black");
    Cube & SetWireWidth (int          Width=1.0);

    // Methods
    void AddActorsTo (VTKWin & Win) { Win.AddActor(_cube_actor);  Win.AddActor(_wire_actor); }

private:
    vtkCubeSource     * _cube;
    vtkPolyDataMapper * _cube_mapper;
    vtkActor          * _cube_actor;
    vtkPolyDataMapper * _wire_mapper;
    vtkActor          * _wire_actor;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Cube::Cube (Vec3_t const & Cen, double Lx, double Ly, double Lz)
{
    // create object
    _cube        = vtkCubeSource     ::New();
    _cube_mapper = vtkPolyDataMapper ::New();
    _cube_actor  = vtkActor          ::New();
    _cube        -> SetCenter          (Cen(0), Cen(1), Cen(2));
    _cube        -> SetXLength         (Lx);
    _cube        -> SetYLength         (Ly);
    _cube        -> SetZLength         (Lz);
    _cube_mapper -> SetInputConnection (_cube->GetOutputPort());
    _cube_actor  -> SetMapper          (_cube_mapper);
    SetColor ();

    // borders
    _wire_mapper = vtkPolyDataMapper ::New();
    _wire_actor  = vtkActor          ::New();
    _wire_mapper -> SetInput            (_cube->GetOutput());
    _wire_mapper -> ScalarVisibilityOff ();
    _wire_actor  -> SetMapper           (_wire_mapper);
    _wire_actor  -> GetProperty         ()->SetRepresentationToWireframe();

    _cube_mapper->SetResolveCoincidentTopologyPolygonOffsetParameters(0,1);
    _cube_mapper->SetResolveCoincidentTopologyToPolygonOffset();
    _wire_mapper->SetResolveCoincidentTopologyPolygonOffsetParameters(1,1);
    _wire_mapper->SetResolveCoincidentTopologyToPolygonOffset();

    // same color for inside and outside edges
    _wire_mapper->ScalarVisibilityOff();
    _wire_actor->GetProperty()->SetAmbient(1.0);
    _wire_actor->GetProperty()->SetDiffuse(0.0);
    _wire_actor->GetProperty()->SetSpecular(0.0);
    SetWireWidth ();
    SetWireColor ();
}

inline Cube::~Cube ()
{
    _cube        -> Delete();
    _cube_mapper -> Delete();
    _cube_actor  -> Delete();
    _wire_mapper -> Delete();
    _wire_actor  -> Delete();
}

inline Cube & Cube::SetColor (char const * Name, double Opacity)
{
    Vec3_t c(Colors::Get(Name));
    _cube_actor->GetProperty()->SetColor   (c(0), c(1), c(2));
    _cube_actor->GetProperty()->SetOpacity (Opacity);
    return (*this);
}

inline Cube & Cube::SetWireColor (char const * Name) 
{
    Vec3_t c(Colors::Get(Name));
    _wire_actor->GetProperty()->SetColor (c(0),c(1),c(2));
    return (*this); 
}

inline Cube & Cube::SetWireWidth (int Width)
{
    _wire_actor->GetProperty()->SetLineWidth(Width);
    return (*this); 
}

#endif // MECHSYS_CUBE_H
