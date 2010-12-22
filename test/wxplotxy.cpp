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

#define USE_WXWIDGETS

// MechSys
#include <mechsys/gui/plotxy.h>
#include <mechsys/gui/pixmaps/icon.xpm>

class MyFrame: public wxFrame
{
public:
    // Constructor and Destructor
     MyFrame (wxString const & Title);
    ~MyFrame () { Aui.UnInit(); }

    // Methods
    void ReBuild ();

    // Data
    wxAuiManager   Aui;
    GUI::PlotXY  * P0;
    GUI::PlotXY  * P1;
    Array<double>  X, Y0, Y1, Y2, Y3, Y4, Y5;

    // Controls' data
    bool ShowGrid,  EqScales;
    int  xNumTicks, yNumTicks;

    // Events
    void OnQuit    (wxCommandEvent & event) { Close (true); }
    void OnAbout   (wxCommandEvent & event) { WXMSG ("Testing PlotXY"); }
    void OnShowAll (wxCommandEvent & Event) { SHOW_ALL_WXPANES (Aui); }
    void OnRun     (wxCommandEvent & Event) { WXMSG ("Method not implemented"); }
    void OnAdd     (wxCommandEvent & Event);
    void OnDel     (wxCommandEvent & Event);
    void OnReBuild (wxCommandEvent & Event) { ReBuild(); }
    DECLARE_EVENT_TABLE()
};

#define MYAPP_TITLE "Testing PlotXY"

// MechSys
#include <mechsys/gui/wxmyapp.h>


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


enum
{
    ID_MNU_RUN = wxID_HIGHEST+1,
    ID_MNU_SHOWALL,
    ID_XNTCKS,
    ID_YNTCKS,
    ID_GRID,
    ID_EQSF,
    ID_ADD,
    ID_DEL,
};

BEGIN_EVENT_TABLE (MyFrame, wxFrame)
    EVT_MENU     (wxID_EXIT,      MyFrame::OnQuit)
    EVT_MENU     (wxID_ABOUT,     MyFrame::OnAbout)
    EVT_MENU     (ID_MNU_SHOWALL, MyFrame::OnShowAll)
    EVT_MENU     (ID_MNU_RUN,     MyFrame::OnRun)
    EVT_BUTTON   (ID_ADD,         MyFrame::OnAdd)
    EVT_BUTTON   (ID_DEL,         MyFrame::OnDel)
    EVT_CHECKBOX (ID_GRID,        MyFrame::OnReBuild)
    EVT_CHECKBOX (ID_EQSF,        MyFrame::OnReBuild)
END_EVENT_TABLE ()

MyFrame::MyFrame (wxString const & Title)
   : wxFrame(NULL, wxID_ANY, Title), ShowGrid(true), EqScales(false)
{
    // force validation of child controls
    SetExtraStyle (wxWS_EX_VALIDATE_RECURSIVELY);

    // tell wxAuiManager to manage this frame
    Aui.SetManagedWindow (this);
    SetMinSize (wxSize(640,480));
    SetIcon    (wxIcon(icon_xpm));

    // menu and status bars
    ADD_WXMENUSTATUS (mnu, mnu_fil, mnu_wnd, mnu_run, mnu_hlp);

    // control panel
    ADD_WXPANEL     (pnl, szt, szr, 4, 2);
    ADD_WXNUMINPUT2 (pnl, szr, ID_XNTCKS, c0, "x num ticks",  xNumTicks);
    ADD_WXNUMINPUT2 (pnl, szr, ID_YNTCKS, c1, "y num ticks",  yNumTicks);
    ADD_WXCHECKBOX2 (pnl, szr, ID_GRID,   c2, "Show grid",    ShowGrid);
    ADD_WXCHECKBOX2 (pnl, szr, ID_EQSF,   c3, "Equal scales", EqScales);
    ADD_WXBUTTON    (pnl, szt, ID_ADD,    c4, "Add curves");
    ADD_WXBUTTON    (pnl, szt, ID_DEL,    c5, "Delete curves");

    // P0 and P1
    P0 = new GUI::PlotXY (this, "First Graph", "x", "y(x)");
    P1 = new GUI::PlotXY (this, "Second Graph", "x", "y(x)");

    // commit all changes to wxAuiManager
    Aui.AddPane (pnl, wxAuiPaneInfo().Name("pnl").Caption("Control panel").DestroyOnClose(false).Left().MinSize(wxSize(180,100)));
    Aui.AddPane (P0,  wxAuiPaneInfo().Name("P0") .Caption("First graph")  .DestroyOnClose(false).Centre());
    Aui.AddPane (P1,  wxAuiPaneInfo().Name("P1") .Caption("Second graph") .DestroyOnClose(false).Centre().Position(1));
    Aui.Update  ();

    // data for plots
    const int np = 20;
    X .Resize(np+1);
    Y0.Resize(np+1);
    Y1.Resize(np+1);
    Y2.Resize(np+1);
    Y3.Resize(np+1);
    Y4.Resize(np+1);
    Y5.Resize(np+1);
    for (int i=0; i<np+1; ++i)
    {
        X [i] = static_cast<double>(i)/np;
        Y0[i] = 2.0*X[i];
        Y1[i] = X[i]*X[i];
        Y2[i] = X[i];
        Y3[i] = 0.5*X[i];
        Y4[i] = log(1.0+X[i]);
        Y5[i] = sqrt(X[i]);
    }
}

void MyFrame::OnAdd (wxCommandEvent & Event)
{
    // P0
    P0->AddCurve (&X, &Y0, "2*x").Pen.Set("red", "solid",2);  P0->C[0].Pch=1;
    P0->AddCurve (&X, &Y1, "x^2").Pen.Set("blue","dot",  1);  P0->C[1].Typ=GUI::CT_BOTH;  P0->C[1].Pch=2;

    // P1
    P1->AddCurve (&X, &Y2, "x")       .Pen.Set("red",    "dash" , 2);  P1->C[0].Typ=GUI::CT_LINES;
    P1->AddCurve (&X, &Y3, "x/2")     .Pen.Set("dgreen", "solid", 3);  P1->C[1].Typ=GUI::CT_BOTH;  P1->C[1].Pch=100;
    P1->AddCurve (&X, &Y4, "log(1+x)").Pen.Set("blue",   "solid", 1);  P1->C[2].Typ=GUI::CT_BOTH;  P1->C[2].Pch=3;
    P1->AddCurve (&X, &Y5, "sqrt(x)") .Pen.Set("orange", "solid", 2);  P1->C[3].Typ=GUI::CT_BOTH;  P1->C[3].Pch=4;

    ReBuild ();
}

void MyFrame::OnDel (wxCommandEvent & Event)
{
    P0->DelCurves ();
    P1->DelCurves ();
    ReBuild ();
}

void MyFrame::ReBuild ()
{
    // update control's data
    TransferDataFromWindow ();

    // number of ticks
    if (xNumTicks>0) { P0->BNumTck = xNumTicks;   P1->BNumTck = xNumTicks; }
    if (yNumTicks>0) { P0->LNumTck = yNumTicks;   P1->LNumTck = yNumTicks; }

    // show grid
    if (ShowGrid)
    {
        P0->Grid = true;
        P1->Grid = true;
    }
    else
    {
        P0->Grid = false;
        P1->Grid = false;
    }

    // equal scales
    if (EqScales)
    {
        P0->EqSF = true;
        P1->EqSF = true;
    }
    else
    {
        P0->EqSF = false;
        P1->EqSF = false;
    }

    // redraw
    P0->Redraw ();
    P1->Redraw ();
}
