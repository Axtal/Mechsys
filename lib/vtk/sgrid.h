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

#ifndef MECHSYS_SGRID_H
#define MECHSYS_SGRID_H

// VTK
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkTextActor3D.h>
#include <vtkProperty.h>

// MechSys
#include <mechsys/vtk/win.h>
#include <mechsys/util/array.h>
#include <mechsys/util/colors.h>
#include <mechsys/util/string.h>
#include <mechsys/linalg/matvec.h>

namespace VTK
{

class SGrid
{
public:
    // Constructor & Destructor
     SGrid () { _create(); }
    ~SGrid ();

    // Alternative constructors
    SGrid (int N[3], double L[6])             { _create();  Resize(N,L); }
    SGrid (int N[3], Array<double> const & L) { _create();  Resize(N,L); }
    SGrid (int Nx, double Xmin, double Xmax,
           int Ny, double Ymin, double Ymax,
           int Nz, double Zmin, double Zmax) { _create();  Resize(Nx,Xmin,Xmax, Ny,Ymin,Ymax, Nz,Zmin,Zmax); }

    // Set methods
    SGrid & Resize   (int N[3], double L[6])             { return Resize(N[0],L[0],L[1], N[1],L[2],L[3], N[2],L[4],L[5]); }
    SGrid & Resize   (int N[3], Array<double> const & L) { return Resize(N[0],L[0],L[1], N[1],L[2],L[3], N[2],L[4],L[5]); }
    SGrid & Resize   (int Nx, double Xmin, double Xmax,
                      int Ny, double Ymin, double Ymax,
                      int Nz, double Zmin, double Zmax);
    SGrid & SetColor (char const * Name="black", double Opacity=1.0);

    // Methods
    void ShowWire   ()             { _sgrid_actor->GetProperty()->SetRepresentationToWireframe(); }
    void ShowPoints (int PtSize=4) { _sgrid_actor->GetProperty()->SetRepresentationToPoints();  _sgrid_actor->GetProperty()->SetPointSize(PtSize); }
    void ShowIds    (double OriX=90, double OriY=90, double OriZ=45, double Scale=0.003, int SizePt=14, bool Shadow=true, char const * Color="blue");
    void AddTo      (VTK::Win & win);

private:
    vtkPoints              * _points;
    vtkStructuredGrid      * _sgrid;
    vtkDataSetMapper       * _sgrid_mapper;
    vtkActor               * _sgrid_actor;
    Array<vtkTextActor3D*>   _text;
    void _create ();
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline SGrid & SGrid::Resize (int Nx, double Xmin, double Xmax, int Ny, double Ymin, double Ymax, int Nz, double Zmin, double Zmax)
{
    if (Nx<2) throw new Fatal("SGrid::Resize: Nx==N[0]=%d must be greater than 1",Nx);
    if (Ny<2) throw new Fatal("SGrid::Resize: Ny==N[1]=%d must be greater than 1",Ny);
    if (Nz<2) throw new Fatal("SGrid::Resize: Nz==N[2]=%d must be greater than 1",Nz);
    _points -> Reset         ();
    _points -> Allocate      (Nx*Ny*Nz);
    _sgrid  -> SetDimensions (Nx,Ny,Nz);
    double dx = (Xmax-Xmin)/(Nx-1.0);
    double dy = (Ymax-Ymin)/(Ny-1.0);
    double dz = (Zmax-Zmin)/(Nz-1.0);
    size_t m  = 0;
    for (int k=0; k<Nz; ++k)
    for (int j=0; j<Ny; ++j)
    for (int i=0; i<Nx; ++i)
        _points->InsertPoint (m++, Xmin+i*dx, Ymin+j*dy, Zmin+k*dz);
    return (*this);
}

inline SGrid::~SGrid ()
{
    _points       -> Delete();
    _sgrid        -> Delete();
    _sgrid_mapper -> Delete();
    _sgrid_actor  -> Delete();
    for (size_t i=0; i<_text.Size(); ++i) _text[i] -> Delete();
}

inline SGrid & SGrid::SetColor (char const * Name, double Opacity)
{
    Vec3_t c = Colors::Get(Name);
    _sgrid_actor->GetProperty()->SetColor   (c.data());
    _sgrid_actor->GetProperty()->SetOpacity (Opacity);
    return (*this);
}

inline void SGrid::ShowIds (double OriX, double OriY, double OriZ, double Scale, int SizePt, bool Shadow, char const * Color)
{
    Vec3_t c(Colors::Get(Color));
    for (size_t i=0; i<_text.Size(); ++i) _text[i] -> Delete();
    String buf;
    _text.Resize (_points->GetNumberOfPoints());
    for (int i=0; i<_points->GetNumberOfPoints(); ++i)
    {
        buf.Printf ("%d",i);
        _text[i] = vtkTextActor3D                   ::New();
        _text[i] -> SetInput                        (buf.CStr());
        _text[i] -> SetPosition                     (_points->GetPoint(i));
        _text[i] -> SetOrientation                  (OriX, OriY, OriZ);
        _text[i] -> SetScale                        (Scale);
        _text[i] -> GetTextProperty()-> SetFontSize (SizePt);
        _text[i] -> GetTextProperty()-> SetShadow   (Shadow);
        _text[i] -> GetTextProperty()-> SetColor    (c.data());
    }
}

inline void SGrid::AddTo (VTK::Win & win)
{
    win.AddActor (_sgrid_actor); 
    for (size_t i=0; i<_text.Size(); ++i) win.AddActor (reinterpret_cast<vtkActor*>(_text[i]));
}

inline void SGrid::_create ()
{
    _points       = vtkPoints         ::New();
    _sgrid        = vtkStructuredGrid ::New();
    _sgrid_mapper = vtkDataSetMapper  ::New();
    _sgrid_actor  = vtkActor          ::New();
    _sgrid        -> SetPoints        (_points);
    _sgrid_mapper -> SetInput         (_sgrid);
    _sgrid_actor  -> SetMapper        (_sgrid_mapper);
    _sgrid_actor  -> GetProperty() -> SetPointSize (4);
    ShowWire ();
    SetColor ();
}

}; // namespace VTK

#endif // MECHSYS_SGRID_H
