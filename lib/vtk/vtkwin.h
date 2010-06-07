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

#ifndef MECHSYS_VTKWIN_H
#define MECHSYS_VTKWIN_H

// VTK
#include <vtkCamera.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkActor.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>

// MechSys
#include <mechsys/util/string.h>
#include <mechsys/util/colors.h>

class VTKWin
{
public:
    // Constructors
    VTKWin (int Width=300, int Height=300);
    VTKWin (vtkRenderWindowInteractor * vtkRWI, int Width=300, int Height=300);

    // Destructor
    ~VTKWin ();

    // Methods
    void AddActor (vtkActor * TheActor);
    void DelActor (vtkActor * TheActor) { _renderer -> RemoveActor(TheActor); }
    void AddLight (vtkLight * Light)    { _renderer -> AddLight(Light);       }
    void Show     ();
    void WritePNG (char const * Filename);

    // Access methods
    vtkRenderer * GetRen    () { return _renderer; }
    vtkCamera   * GetCamera () { return _camera;   }

    // Set methods
    VTKWin & SetViewDefault ();
    VTKWin & SetViewPIplane ();
    VTKWin & SetBgColor     (char const * Name="white") { _bg_clr = Colors::Get(Name);  return (*this); }

private:
    Vec3_t                      _bg_clr;
    vtkCamera                 * _camera;
    vtkRenderer               * _renderer;
    vtkRenderWindow           * _ren_win;
    vtkRenderWindowInteractor * _interactor;
    vtkInteractorStyleSwitch  * _int_switch;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline VTKWin::VTKWin (int Width, int Height)
{
    _camera = vtkCamera::New();
    SetViewDefault();

    _renderer = vtkRenderer::New();
    _ren_win  = vtkRenderWindow::New();
    _ren_win->AddRenderer(_renderer);
    _ren_win->SetSize(Width,Height);

    _interactor = vtkRenderWindowInteractor::New();
    _int_switch = vtkInteractorStyleSwitch::New();
    _interactor -> SetRenderWindow(_ren_win);
    _interactor -> SetInteractorStyle(_int_switch);
    _int_switch -> SetCurrentStyleToTrackballCamera();

    SetBgColor ();
}

inline VTKWin::VTKWin (vtkRenderWindowInteractor * vtkRWI, int Width, int Height)
{
    _camera = vtkCamera::New();
    SetViewDefault();

    _renderer = vtkRenderer::New();
    _ren_win  = vtkRenderWindow::New();
    _ren_win->AddRenderer(_renderer);
    _ren_win->SetSize(Width,Height);

    _interactor = vtkRWI;
    _int_switch = vtkInteractorStyleSwitch::New();
    _interactor -> SetRenderWindow(_ren_win);
    _interactor -> SetInteractorStyle(_int_switch);
    _int_switch -> SetCurrentStyleToTrackballCamera();

    SetBgColor ();
}

inline VTKWin::~VTKWin ()
{
    _camera     -> Delete();
    _renderer   -> Delete();
    _ren_win    -> Delete();
    _interactor -> Delete();
    _int_switch -> Delete();
}

inline void VTKWin::AddActor (vtkActor * TheActor)
{
    _renderer -> AddActor        (TheActor);
    _renderer -> SetActiveCamera (_camera);
    _renderer -> ResetCamera     ();
}

inline void VTKWin::Show ()
{
    _renderer   -> SetBackground (_bg_clr(0), _bg_clr(1), _bg_clr(2));
    _ren_win    -> Render        ();
    _interactor -> Start         ();
}

inline void VTKWin::WritePNG (char const * Filename)
{
    // Window to image filter
    vtkWindowToImageFilter * win_to_img = vtkWindowToImageFilter::New();
    win_to_img -> SetInput(_ren_win);
    win_to_img -> Update();

    // PNG writer
    vtkPNGWriter * writer = vtkPNGWriter::New();
    writer -> SetInput(win_to_img->GetOutput());
    writer -> SetFileName(Filename);
    writer -> Write();

    // Clean up
    win_to_img -> Delete();
    writer     -> Delete();
}

inline VTKWin & VTKWin::SetViewDefault ()
{
    _camera->SetViewUp     (0,0,1);
    _camera->SetPosition   (2,1,1);
    _camera->SetFocalPoint (0,0,0);
    return (*this);
}

inline VTKWin & VTKWin::SetViewPIplane ()
{
    _camera->SetViewUp     (0,0,1);
    _camera->SetPosition   (1,1,1);
    _camera->SetFocalPoint (0,0,0);
    return (*this);
}

#endif // MECHSYS_VTKWIN_H
