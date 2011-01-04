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
#include <vtkDoubleArray.h>
#include <vtkStructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkTextActor3D.h>
#include <vtkProperty.h>
#include <vtkTextProperty.h>
#include <vtkPointData.h>
#include <vtkStructuredGridWriter.h>

// MechSys
#include <mechsys/vtk/win.h>
#include <mechsys/util/array.h>
#include <mechsys/util/colors.h>
#include <mechsys/util/string.h>
#include <mechsys/linalg/matvec.h>

namespace VTK
{

typedef void (*GridCallBack) (Vec3_t const & X, double & F, Vec3_t & V, void * UserData);

class SGrid
{
public:
    // Constructor & Destructor
     SGrid () : _func(NULL), _udat(NULL) { _create(); }
    ~SGrid ();

    // Alternative constructors
    SGrid (int N[3], double L[6], GridCallBack Func=NULL, void * UserData=NULL);
    SGrid (Array<int> const & N, Array<double> const & L, GridCallBack Func=NULL, void * UserData=NULL);
    SGrid (int Nx, double Xmin, double Xmax,
           int Ny, double Ymin, double Ymax,
           int Nz, double Zmin, double Zmax, GridCallBack Func=NULL, void * UserData=NULL);

    // Set methods
    SGrid & Resize   (int N[3], double L[6])                         { return Resize(N[0],L[0],L[1], N[1],L[2],L[3], N[2],L[4],L[5]); }
    SGrid & Resize   (Array<int> const & N, Array<double> const & L) { return Resize(N[0],L[0],L[1], N[1],L[2],L[3], N[2],L[4],L[5]); }
    SGrid & Resize   (int Nx, double Xmin, double Xmax,
                      int Ny, double Ymin, double Ymax,
                      int Nz, double Zmin, double Zmax);
    SGrid & SetColor (char const * Name="black", double Opacity=1.0);
    SGrid & SetFunc  (GridCallBack Func, void * UserData=NULL) { _func=Func;  _udat=UserData;  _calc_f();  return (*this); }

    // Access methods
    int                 Size       ()                  const { return _points->GetNumberOfPoints(); }
    void                GetPoint   (int i, Vec3_t & x) const { _points->GetPoint (i, x.data());     }
    void                SetPoint   (int i, Vec3_t & x)       { _points->SetPoint (i, x.data());     }
    vtkStructuredGrid * GetGrid    ()                        { return _sgrid;                       }
    vtkPoints         * GetPoints  ()                        { return _points;                      }
    vtkDoubleArray    * GetScalars ()                        { return _scalars;                     }
    vtkDoubleArray    * GetVectors ()                        { return _vectors;                     }
    
    // Methods
    void ShowWire   ()             { _sgrid_actor->GetProperty()->SetRepresentationToWireframe(); }
    void ShowPoints (int PtSize=4) { _sgrid_actor->GetProperty()->SetRepresentationToPoints();  _sgrid_actor->GetProperty()->SetPointSize(PtSize); }
    void ShowIds    (double OriX=90, double OriY=90, double OriZ=45, double Scale=0.003, int SizePt=14, bool Shadow=true, char const * Color="blue");
    void AddTo      (VTK::Win & win);
    void WriteVTK   (char const * Filekey);
    void FilterV    (double F=0.0, double Tol=1.0e-3);

private:
    GridCallBack             _func;
    void                   * _udat;
    vtkPoints              * _points;
    vtkDoubleArray         * _scalars;
    vtkDoubleArray         * _vectors;
    vtkStructuredGrid      * _sgrid;
    vtkDataSetMapper       * _sgrid_mapper;
    vtkActor               * _sgrid_actor;
    Array<vtkTextActor3D*>   _text;
    void _create ();
    void _calc_f ();
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline SGrid::SGrid (int N[3], double L[6], GridCallBack Func, void * UserData)
    : _func(Func), _udat(UserData) 
{
    _create ();
    Resize  (N,L);
}

inline SGrid::SGrid (Array<int> const & N, Array<double> const & L, GridCallBack Func, void * UserData)
    : _func(Func), _udat(UserData) 
{
    _create ();
    Resize  (N,L);
}

inline SGrid::SGrid (int Nx, double Xmin, double Xmax, int Ny, double Ymin, double Ymax, int Nz, double Zmin, double Zmax, GridCallBack Func, void * UserData)
    : _func(Func), _udat(UserData) 
{
    _create ();
    Resize  (Nx,Xmin,Xmax, Ny,Ymin,Ymax, Nz,Zmin,Zmax);
}

inline SGrid & SGrid::Resize (int Nx, double Xmin, double Xmax, int Ny, double Ymin, double Ymax, int Nz, double Zmin, double Zmax)
{
    if (Nx<2) throw new Fatal("SGrid::Resize: Nx==N[0]=%d must be greater than 1",Nx);
    if (Ny<2) throw new Fatal("SGrid::Resize: Ny==N[1]=%d must be greater than 1",Ny);
    if (Nz<2) throw new Fatal("SGrid::Resize: Nz==N[2]=%d must be greater than 1",Nz);
    _points   -> Reset                 ();
    _points   -> Allocate              (Nx*Ny*Nz);
    _scalars  -> Reset                 ();
    _scalars  -> Allocate              (Nx*Ny*Nz);
    _vectors  -> Reset                 ();
    _vectors  -> SetNumberOfComponents (3);
    _vectors  -> SetNumberOfTuples     (Nx*Ny*Nz);
    _sgrid    -> SetDimensions         (Nx,Ny,Nz);
    double dx  = (Xmax-Xmin)/(Nx-1.0);
    double dy  = (Ymax-Ymin)/(Ny-1.0);
    double dz  = (Zmax-Zmin)/(Nz-1.0);
    size_t idx = 0;
    double f   = 0.0;
    Vec3_t x, v;
    for (int k=0; k<Nz; ++k)
    for (int j=0; j<Ny; ++j)
    for (int i=0; i<Nx; ++i)
    {
        x = Xmin+i*dx, Ymin+j*dy, Zmin+k*dz;
        _points -> InsertPoint (idx, x.data());
        if (_func==NULL)
        {
            _scalars -> InsertTuple1 (idx, 0);
            _vectors -> InsertTuple3 (idx, 0,0,0);
        }
        else
        {
            (*_func) (x, f, v, _udat);
            _scalars -> InsertTuple1 (idx, f);
            _vectors -> InsertTuple3 (idx, v(0), v(1), v(2));
        }
        idx++;
    }
    _sgrid -> GetPointData() -> SetScalars (_scalars);
    _sgrid -> GetPointData() -> SetVectors (_vectors);
    return (*this);
}

inline SGrid::~SGrid ()
{
    _points       -> Delete();
    _scalars      -> Delete();
    _vectors      -> Delete();
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

inline void SGrid::WriteVTK (char const * Filekey)
{
    String buf(Filekey);
    buf.append(".vtk");
    vtkStructuredGridWriter * writer = vtkStructuredGridWriter::New();
    writer -> SetInput    (_sgrid);
    writer -> SetFileName (buf.CStr());
    writer -> Write       ();
    writer -> Delete      ();
}

inline void SGrid::FilterV (double F, double Tol)
{
    for (int i=0; i<_scalars->GetNumberOfTuples(); ++i)
    {
        double   f = _scalars->GetTuple1(i);
        //double * v = _vectors->GetTuple3(i);
        if (fabs(f-F)>Tol) _vectors->SetTuple3 (i, 0.0, 0.0, 0.0);
    }
}

inline void SGrid::_create ()
{
    _points       = vtkPoints         ::New();
    _scalars      = vtkDoubleArray    ::New();
    _vectors      = vtkDoubleArray    ::New();
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

inline void SGrid::_calc_f ()
{
    double f;
    Vec3_t x, v;
    for (int i=0; i<Size(); ++i)
    {
        GetPoint (i, x);
        if (_func==NULL)
        {
            _scalars -> SetTuple1 (i, 0);
            _vectors -> SetTuple3 (i, 0,0,0);
        }
        else
        {
            (*_func) (x, f, v, _udat);
            _scalars -> SetTuple1 (i, f);
            _vectors -> SetTuple3 (i, v(0), v(1), v(2));
        }
    }
}

}; // namespace VTK

#endif // MECHSYS_SGRID_H
