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

// wxWidgets
#include "wx/aui/aui.h"
#include "wx/menu.h"
#include "wx/msgdlg.h"
#include "wx/artprov.h"

#define USE_WXWIDGETS

// MechSys
#include "gui/plotxy.h"
#include "gui/pixmaps/icon.xpm"

class MyFrame: public wxFrame
{
public:
    // Constructor and Destructor
     MyFrame (const wxString & Title);
    ~MyFrame () { Aui.UnInit(); }

    // Methods
    void OnQuit    (wxCommandEvent & event) { Close(true); }
    void OnAbout   (wxCommandEvent & event) { wxMessageBox( _("Testing PlotXY"), _("About wxplotxy"), wxOK | wxICON_INFORMATION, this ); }
    void OnShowAll (wxCommandEvent & Event);
    void OnRun     (wxCommandEvent & Event);

    // Data
    wxAuiManager   Aui;
    GUI::PlotXY  * P0;
    GUI::PlotXY  * P1;
    Array<double>  X, Y0, Y1, Y2, Y3, Y4, Y5;

    // Events
    DECLARE_EVENT_TABLE()
};

#define MYAPP_TITLE "Testing PlotXY"
#include "gui/wxmyapp.h"


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


enum
{
    ID_RUN = wxID_HIGHEST+1,
    ID_SHOW_ALL,
};

BEGIN_EVENT_TABLE (MyFrame, wxFrame)
    EVT_MENU (wxID_EXIT,   MyFrame::OnQuit)
    EVT_MENU (wxID_ABOUT,  MyFrame::OnAbout)
    EVT_MENU (ID_SHOW_ALL, MyFrame::OnShowAll)
    EVT_MENU (ID_RUN,      MyFrame::OnRun)
END_EVENT_TABLE ()


MyFrame::MyFrame (const wxString & Title)
   : wxFrame(NULL, wxID_ANY, Title)
{
    // settings
    Aui.SetManagedWindow (this);   // tell wxAuiManager to manage this frame
    SetMinSize (wxSize(400,300));  // minimum size
    SetIcon    (wxIcon(icon_xpm)); // set frame icon

    // menu and status bars
    wxMenuBar * mnu      = new wxMenuBar;
    wxMenu    * mnu_file = new wxMenu;
    wxMenu    * mnu_wnd  = new wxMenu;
    wxMenu    * mnu_run  = new wxMenu;
    wxMenu    * mnu_help = new wxMenu;
    mnu_file -> Append (wxID_EXIT   , _("&Quit\tCtrl-Q") , _("Quit this program"));
    mnu_wnd  -> Append (ID_SHOW_ALL , _("&Show All")     , _("Show all windows"));
    mnu_run  -> Append (ID_RUN      , _("&Run\tCtrl-R")  , _("Run simulation"));
    mnu_help -> Append (wxID_ABOUT  , _("&About\tF1")    , _("Show about dialog"));
    mnu      -> Append (mnu_file    , _("&File"));
    mnu      -> Append (mnu_wnd     , _("&Window"));
    mnu      -> Append (mnu_run     , _("&Run"));
    mnu      -> Append (mnu_help    , _("&Help"));
    SetMenuBar      (mnu);
    CreateStatusBar (2);
    SetStatusText   (_("Welcome"));

    // toolbar
    wxAuiToolBar * tb = new wxAuiToolBar(this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxAUI_TB_DEFAULT_STYLE | wxAUI_TB_OVERFLOW);
    tb->SetToolBitmapSize (wxSize(48,48));
    tb->AddTool (wxID_EXIT, _("Quit"), wxArtProvider::GetBitmap(wxART_ERROR));
    tb->Realize();

    // panel
    wxPanel * pnl = new wxPanel (this, wxID_ANY, wxDefaultPosition, wxDefaultSize);

    // plots
    const int np = 20;
    X .Resize(np);
    Y0.Resize(np);
    Y1.Resize(np);
    Y2.Resize(np);
    Y3.Resize(np);
    Y4.Resize(np);
    Y5.Resize(np);
    for (int i=0; i<np; ++i)
    {
        X [i] = static_cast<double>(i)/np;
        Y0[i] = 2.0*X[i];
        Y1[i] = X[i]*X[i];
        Y2[i] = X[i];
        Y3[i] = 0.5*X[i];
        Y4[i] = log(1.0+X[i]);
        Y5[i] = sqrt(X[i]);
    }

    // P0
    P0 = new GUI::PlotXY (this);
    P0->EqScales (false).RecalcSF(true);

    // P1
    P1 = new GUI::PlotXY (this);
    P1->EqScales (false).RecalcSF(true);

    // panes
    Aui.AddPane (tb,  wxAuiPaneInfo().Name(_("tb" )).Caption(_("Toolbar"      )).DestroyOnClose(false).ToolbarPane().Top().LeftDockable(false).RightDockable(false));
    Aui.AddPane (pnl, wxAuiPaneInfo().Name(_("pnl")).Caption(_("Control panel")).DestroyOnClose(false).Left());
    Aui.AddPane (P0,  wxAuiPaneInfo().Name(_("P0")) .Caption(_("First graph"))  .DestroyOnClose(false).Centre());
    Aui.AddPane (P1,  wxAuiPaneInfo().Name(_("P1")) .Caption(_("Second graph")) .DestroyOnClose(false).Centre().Position(1));

    // commit all changes to wxAuiManager
    Aui.Update();
}

void MyFrame::OnShowAll (wxCommandEvent & Event)
{
    wxAuiPaneInfoArray & all_panes = Aui.GetAllPanes();
    for (size_t i=0; i<all_panes.GetCount(); ++i)
    {
        all_panes.Item(i).Show(true);
        all_panes.Item(i).Dock();
    }
    Aui.Update();
}

void MyFrame::OnRun (wxCommandEvent & Event)
{
    P0->AddCurve (&X, &Y0, "2*x")     .Pen.Set("red",     "solid", 2);  P0->C[0].Pch=16;
    P0->AddCurve (&X, &Y1, "x^2")     .Pen.Set("blue",    "dot",   1);  P0->C[1].Typ=GUI::CT_BOTH;
    P1->AddCurve (&X, &Y2, "x")       .Pen.Set("red",     "dash" , 2);  P1->C[0].Typ=GUI::CT_LINES;
    P1->AddCurve (&X, &Y3, "x/2")     .Pen.Set("dgreen",  "solid", 3);  P1->C[1].Typ=GUI::CT_BOTH;
    P1->AddCurve (&X, &Y4, "log(1+x)").Pen.Set("blue",    "solid", 1);  P1->C[2].Pch=16;
    P1->AddCurve (&X, &Y5, "sqrt(x)") .Pen.Set("dmagenta","solid", 2);  P1->C[3].Typ=GUI::CT_BOTH;
}
