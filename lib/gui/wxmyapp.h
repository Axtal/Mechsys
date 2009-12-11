/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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

#ifndef MECHSYS_GUI_WXMYAPP_H
#define MECHSYS_GUI_WXMYAPP_H

// wxWidgets
#include "wx/app.h"
#include "wx/utils.h"
#include "wx/msgdlg.h"
#include "wx/log.h"

// MechSys
#include "util/fatal.h"

class MyApp : public wxApp
{
public:
    // program startup
    virtual bool OnInit();

    // 2nd-level exception handling: we get all the exceptions occurring in any
    // event handler here
    virtual bool OnExceptionInMainLoop();

    // 3rd, and final, level exception handling: whenever an unhandled
    // exception is caught, this function is called
    virtual void OnUnhandledException();

    // and now for something different: this function is called in case of a
    // crash (e.g. dereferencing null pointer, division by 0, ...)
    virtual void OnFatalException();
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


bool MyApp::OnInit()
{
    MyFrame * frame = new MyFrame(MYAPP_TITLE);
    frame->Show   (true);
    SetTopWindow  (frame);
    return true;
}

bool MyApp::OnExceptionInMainLoop()
{
    String msg;
    try
    {
        msg = "MyApp::OnExceptionInMainLoop:: Internal Error";
        throw;
    }
    catch (Fatal      * e)     { msg.Printf("Fatal: %s",e->Msg().CStr());  delete e; }
    catch (char const * m)     { msg.Printf("Fatal: %s",m); }
    catch (std::exception & e) { msg.Printf("Fatal: %s",e.what()); }

    wxMessageBox(msg.CStr(), _("Exception caught."), wxOK | wxICON_ERROR);

    return true;
}

void MyApp::OnUnhandledException()
{
    // this shows how we may let some exception propagate uncaught
    try
    {
        throw;
    }
    catch ( ... )
    {
        wxMessageBox(_("Unhandled exception caught, program will terminate."),
                     _("Exception caught."), wxOK | wxICON_ERROR);
    }
}

void MyApp::OnFatalException()
{
    wxMessageBox(_("Program has crashed and will terminate."),
                 _("Exception caught."), wxOK | wxICON_ERROR);
}

DECLARE_APP   (MyApp)
IMPLEMENT_APP (MyApp)

#endif // MECHSYS_GUI_WXMYAPP_H
