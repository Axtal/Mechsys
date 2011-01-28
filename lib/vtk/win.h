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

namespace VTK
{

class Win
{
public:
    // Constructors
    Win (int Width=300, int Height=300);
    Win (vtkRenderWindowInteractor * vtkRWI, int Width=300, int Height=300);

    // Destructor
    ~Win ();

    // Methods
    void AddActor (vtkActor * TheActor, bool RstCam=true);
    void DelActor (vtkActor * TheActor) { _renderer -> RemoveActor(TheActor); }
    void AddLight (vtkLight * Light)    { _renderer -> AddLight(Light);       }
    void Render   ();
    void Show     ();
    void WritePNG (char const * Filename);
    void Camera   (double xUp, double yUp, double zUp, double xFoc, double yFoc, double zFoc, double xPos, double yPos, double zPos);
    void Parallel (bool ParallelProjection=true) { _camera->SetParallelProjection(ParallelProjection); }

    // (low level) Access methods
    vtkRenderer * GetRen    () { return _renderer; }
    vtkCamera   * GetCamera () { return _camera;   }

    // Set methods
    Win & SetViewDefault (bool RevCam=false);
    Win & SetViewPIplane (bool RevCam=false);
    Win & SetBgColor     (char const * Name="white") { _bg_clr = Colors::Get(Name);  return (*this); }

private:
    Vec3_t                      _bg_clr;
    vtkCamera                 * _camera;
    vtkRenderer               * _renderer;
    vtkRenderWindow           * _ren_win;
    vtkRenderWindowInteractor * _interactor;
    vtkInteractorStyleSwitch  * _int_switch;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Win::Win (int Width, int Height)
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

inline Win::Win (vtkRenderWindowInteractor * vtkRWI, int Width, int Height)
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

inline Win::~Win ()
{
    _camera     -> Delete();
    _renderer   -> Delete();
    _ren_win    -> Delete();
    _interactor -> Delete();
    _int_switch -> Delete();
}

inline void Win::AddActor (vtkActor * TheActor, bool RstCam)
{
    _renderer -> AddActor (TheActor);
    if (RstCam)
    {
        _renderer -> SetActiveCamera (_camera);
        _renderer -> ResetCamera     ();
    }
}

inline void Win::Render ()
{
    _renderer -> SetBackground (_bg_clr(0), _bg_clr(1), _bg_clr(2));
    _ren_win  -> Render        ();
}

inline void Win::Show ()
{
    Render ();
    _interactor->Start ();
}

inline void Win::WritePNG (char const * Filename)
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

inline void Win::Camera (double xUp, double yUp, double zUp, double xFoc, double yFoc, double zFoc, double xPos, double yPos, double zPos)
{
    _camera->SetViewUp     (xUp, yUp, zUp);
    _camera->SetFocalPoint (xFoc,yFoc,zFoc);
    _camera->SetPosition   (xPos,yPos,zPos);
    _renderer->ResetCamera ();
}

inline Win & Win::SetViewDefault (bool RevCam)
{
    double c = (RevCam ? -1 : 1);
    _camera->SetViewUp     (0,0,c);
    _camera->SetPosition   (2*c,c,c);
    _camera->SetFocalPoint (0,0,0);
    return (*this);
}

inline Win & Win::SetViewPIplane (bool RevCam)
{
    double c = (RevCam ? -1 : 1);
    _camera->SetViewUp     (0,0,c);
    _camera->SetPosition   (c,c,c);
    _camera->SetFocalPoint (0,0,0);
    return (*this);
}

}; // namespace VTK

#endif // MECHSYS_VTKWIN_H
