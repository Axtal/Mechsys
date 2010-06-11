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

#ifndef MECHSYS_ISOSURF_H
#define MECHSYS_ISOSURF_H

// VTK
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkMarchingContourFilter.h>
#include <vtkStructuredGridWriter.h>
#include <vtkStructuredGridOutlineFilter.h>
#include <vtkHedgeHog.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkLookupTable.h>
#include <vtkProperty.h>

// MechSys
#include <mechsys/vtk/win.h>
#include <mechsys/util/array.h>
#include <mechsys/util/meshgrid.h>
#include <mechsys/linalg/matvec.h>

namespace VTK
{

typedef void (*GridCallBack) (Vec3_t X, double & F, Vec3_t & V, void * UserData);

class IsoSurf
{
public:
    // Constructor & Destructor
    IsoSurf (MeshGrid const & MG, GridCallBack Func, void * UserData=NULL, double PointSize=3.0);

    // Destructor
    ~IsoSurf ();

    // Set methods
    void SetIsoColors (Array<String> const & Names, double Opacity=1.0);
    void SetIsoColor  (char const * Name="blue",    double Opacity=1.0) { SetIsoColors(Array<String>(Name,true), Opacity); }
    void SetIsoValue  (double F=0.0)                                    { _isosurf->SetValue       (0,F);                  }
    void GenIsoValues (int NSurfs, double fMin, double fMax)            { _isosurf->GenerateValues (NSurfs,fMin,fMax);     }
    void SetVecScale  (double Factor=1.0)                               { _hedgehog->SetScaleFactor(Factor);               }
    void SetVecF      (double F=0.0, double Tol=1.0e-3);

    // Methods
    void AddTo     (VTK::Win & win);
    void WriteFile (char const * Filename);

    // Data
    bool ShowPoints;
    bool ShowIsoSurf;
    bool ShowVectors;
    bool ShowOutline;

private:
    vtkPoints                      * _points;
    vtkDoubleArray                 * _scalars;
    vtkDoubleArray                 * _vectors;
    vtkStructuredGrid              * _sgrid;
    vtkDataSetMapper               * _sgrid_mapper;
    vtkActor                       * _sgrid_actor;
    vtkMarchingContourFilter       * _isosurf;
    vtkPolyDataMapper              * _isosurf_mapper;
    vtkActor                       * _isosurf_actor;
	vtkLookupTable                 * _isosurf_lt;
    vtkHedgeHog                    * _hedgehog;
    vtkPolyDataMapper              * _hedgehog_mapper;
    vtkActor                       * _hedgehog_actor;
    vtkStructuredGridOutlineFilter * _outline;
    vtkPolyDataMapper              * _outline_mapper;
    vtkActor                       * _outline_actor;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline IsoSurf::IsoSurf (MeshGrid const & MG, GridCallBack Func, void * UserData, double PointSize)
    : ShowPoints  (false),
      ShowIsoSurf (true),
      ShowVectors (true),
      ShowOutline (false)
{
    // create VTK points, scalars and vectors
    _points  = vtkPoints              ::New();
    _scalars = vtkDoubleArray         ::New();
    _vectors = vtkDoubleArray         ::New();
    _points  -> Allocate              (MG.Length());
    _scalars -> Allocate              (MG.Length());
    _vectors -> SetNumberOfComponents (3);
    _vectors -> SetNumberOfTuples     (MG.Length());

    // fill structures
    double f;
    Vec3_t x, v;
    for (int i=0; i<MG.Length(); ++i)
    {
        x = MG.X[i], MG.Y[i], MG.Z[i];
        (*Func) (x, f, v, UserData);
        _points  -> InsertPoint  (i, x(0), x(1), x(2));
        _scalars -> InsertTuple1 (i, f);
        _vectors -> InsertTuple3 (i, v(0), v(1), v(2));
    }

    // create VTK StructuredGrid
    _sgrid        = vtkStructuredGrid ::New();
    _sgrid_mapper = vtkDataSetMapper  ::New();
    _sgrid_actor  = vtkActor          ::New();
    _sgrid        -> SetDimensions    (MG.nX(), MG.nY(), MG.nZ());
    _sgrid        -> SetPoints        (_points);
    _sgrid        -> GetPointData     () -> SetScalars(_scalars);
    _sgrid        -> GetPointData     () -> SetVectors(_vectors);
    _sgrid_mapper -> SetInput         (_sgrid);
    _sgrid_actor  -> SetMapper        (_sgrid_mapper);
    _sgrid_actor  -> GetProperty() -> SetRepresentationToPoints ();
    _sgrid_actor  -> GetProperty() -> SetPointSize (PointSize);

    // isosurf
    _isosurf        = vtkMarchingContourFilter ::New();
    _isosurf_mapper = vtkPolyDataMapper        ::New();
    _isosurf_actor  = vtkActor                 ::New();
	_isosurf_lt     = vtkLookupTable           ::New();
    _isosurf        -> SetInput                (_sgrid);
    _isosurf_mapper -> SetInputConnection      (_isosurf->GetOutputPort());
    _isosurf_mapper -> SetLookupTable          (_isosurf_lt);
    _isosurf_actor  -> SetMapper               (_isosurf_mapper);
    SetIsoColor ();
    SetIsoValue ();

    // hedgehog
    _hedgehog        = vtkHedgeHog         ::New();
    _hedgehog_mapper = vtkPolyDataMapper   ::New();
    _hedgehog_actor  = vtkActor            ::New();
    _hedgehog        -> SetInput           (_sgrid);
    _hedgehog_mapper -> SetInputConnection (_hedgehog->GetOutputPort());
    _hedgehog_actor  -> SetMapper          (_hedgehog_mapper);
    SetVecScale ();

    // outline
    _outline        = vtkStructuredGridOutlineFilter ::New();
    _outline_mapper = vtkPolyDataMapper              ::New();
    _outline_actor  = vtkActor                       ::New();
    _outline        -> SetInput                      (_sgrid);
    _outline_mapper -> SetInputConnection            (_outline->GetOutputPort());
    _outline_actor  -> SetMapper                     (_outline_mapper);
    _outline_actor  -> GetProperty                   () -> SetColor(0,0,1); 
}

inline IsoSurf::~IsoSurf ()
{
    _points          -> Delete();
    _scalars         -> Delete();
    _vectors         -> Delete();
    _sgrid           -> Delete();
    _sgrid_mapper    -> Delete();
    _sgrid_actor     -> Delete();
    _outline         -> Delete();
    _outline_mapper  -> Delete();
    _outline_actor   -> Delete();
    _hedgehog        -> Delete();
    _hedgehog_mapper -> Delete();
    _hedgehog_actor  -> Delete();
    _isosurf         -> Delete();
    _isosurf_mapper  -> Delete();
    _isosurf_actor   -> Delete();
    _isosurf_lt      -> Delete();
}

inline void IsoSurf::SetIsoColors (Array<String> const & Colors, double Opacity)
{
    _isosurf_lt->SetNumberOfColors (Colors.Size());
    _isosurf_lt->Build             ();
    for (size_t i=0; i<Colors.Size(); ++i)
    {
        Vec3_t c = Colors::Get(Colors[i]);
        _isosurf_lt->SetTableValue (i, c(0), c(1), c(2));
    }
	_isosurf_actor->GetProperty()->SetOpacity(Opacity);
}

inline void IsoSurf::WriteFile (char const * Filename)
{
    vtkStructuredGridWriter * writer = vtkStructuredGridWriter::New();
    writer -> SetInput    (_sgrid);
    writer -> SetFileName (Filename);
    writer -> Write       ();
    writer -> Delete      ();
}

inline void IsoSurf::AddTo (VTK::Win & win)
{
    if (ShowPoints)  win.AddActor (_sgrid_actor);
    if (ShowIsoSurf) win.AddActor (_isosurf_actor);
    if (ShowVectors) win.AddActor (_hedgehog_actor);
    if (ShowOutline) win.AddActor (_outline_actor);
}

inline void IsoSurf::SetVecF (double F, double Tol)
{
    vtkDoubleArray * scalars = static_cast<vtkDoubleArray*>(_sgrid->GetPointData()->GetScalars());
    vtkDoubleArray * vectors = static_cast<vtkDoubleArray*>(_sgrid->GetPointData()->GetVectors());
    for (int i=0; i<scalars->GetNumberOfTuples(); ++i)
    {
        double   f = scalars->GetTuple1(i);
        //double * v = vectors->GetTuple3(i);
        if (fabs(f)>Tol) vectors->SetTuple3 (i, 0.0, 0.0, 0.0);
    }
}

}; // namespace VTK

#endif // MECHSYS_ISOSURF_H
