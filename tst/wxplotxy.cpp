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
#include "wx/app.h"
#include "wx/aui/aui.h"
#include "wx/menu.h"
#include "wx/msgdlg.h"
#include "wx/artprov.h"

// MechSys
#include "icon.xpm"

class MyApp: public wxApp { virtual bool OnInit(); };

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

    // Data
    wxAuiManager   Aui;
    wxPanel      * Pnl;

    // Events
    DECLARE_EVENT_TABLE()
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


bool MyApp::OnInit()
{
    MyFrame * frame = new MyFrame(_("Hello World"));
    frame->Show   (true);
    SetTopWindow  (frame);
    return true;
}

DECLARE_APP   (MyApp)
IMPLEMENT_APP (MyApp)


enum
{
    ID_RUN = wxID_HIGHEST+1,
    ID_SHOW_ALL,
};

BEGIN_EVENT_TABLE (MyFrame, wxFrame)
    EVT_MENU (wxID_EXIT,   MyFrame::OnQuit)
    EVT_MENU (wxID_ABOUT,  MyFrame::OnAbout)
    EVT_MENU (ID_SHOW_ALL, MyFrame::OnShowAll)
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

    // panes
    Aui.AddPane (pnl, wxAuiPaneInfo().Name(_("Pnl")).Caption(_("Control panel")).DestroyOnClose(false));
    Aui.AddPane (tb,  wxAuiPaneInfo().Name(_("tb" )).Caption(_("Toolbar"      )).DestroyOnClose(false).ToolbarPane().Top().LeftDockable(false).RightDockable(false));

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
