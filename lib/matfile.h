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

#ifndef MECHSYS_MATFILE_H
#define MECHSYS_MATFILE_H

// Std Lib
#include <fstream>

// wxWidgets
#ifdef HAS_WXW
  #include <mechsys/gui/wxdict.h>
  #include <mechsys/gui/common.h>
#endif

// MechSys
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/fem/fem.h>

#ifdef USE_WXWIDGETS
class MatFile : public wxWindow
#else
class MatFile
#endif
{
public:
    // Constructor
#ifdef USE_WXWIDGETS
     MatFile (wxFrame * Parent);
    ~MatFile () { Aui.UnInit(); }
#else
    MatFile () { Defaults(); }
#endif

    // Methods
    void Defaults ();
    void Read     (char const * FileName);

    // Data
    Dict     Prms; ///< Tag => prms

#ifdef USE_WXWIDGETS
    wxAuiManager     Aui;
    wxTextCtrl     * TxtFname;   // control with filename
    Array<String>    Names;      // names
    Array<WxDict*>   Dicts;      // dictionaries with parameters
    wxString         LstDir;     // last directory
    void OnLoad (wxCommandEvent & Event);
    void OnSave (wxCommandEvent & Event);
    DECLARE_EVENT_TABLE();
private:
    void _refresh ();
#endif
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void MatFile::Defaults ()
{
}

inline void MatFile::Read (char const * FileName)
{
    // parse material file
    std::fstream mat_file(FileName, std::ios::in);
    if (!mat_file.is_open()) throw new Fatal("MatFile::Read: Could not open file <%s>",FileName);
}

std::ostream & operator<< (std::ostream & os, MatFile const & IF)
{
    return os;
}

#ifdef USE_WXWIDGETS

enum
{
    ID_MAT_LOAD = wxID_HIGHEST+2000,
    ID_MAT_SAVE ,
};

BEGIN_EVENT_TABLE(MatFile, wxWindow)
    EVT_BUTTON (ID_MAT_LOAD, MatFile::OnLoad)
    EVT_BUTTON (ID_MAT_SAVE, MatFile::OnSave)
END_EVENT_TABLE()

inline MatFile::MatFile (wxFrame * Parent)
    : wxWindow (Parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE)
{
    Defaults();

    // flags for sizers
    long f1 = wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL|wxALL;
    long f2 = wxALIGN_LEFT|wxALL|wxEXPAND;

    // tell wxAuiManager to manage this window
    Aui.SetManagedWindow (this);

    // control panel
    wxPanel         * pnl = new wxPanel    (this, wxID_ANY);
    wxBoxSizer      * szt = new wxBoxSizer (wxVERTICAL);
    wxFlexGridSizer * sz0 = new wxFlexGridSizer (/*rows*/1,/*cols*/3,/*vgap*/0,/*hgap*/0);
    TxtFname = new wxTextCtrl (pnl, wxID_ANY, wxEmptyString,wxDefaultPosition,wxSize(400,10),wxTE_READONLY);
    sz0->Add (  new wxButton   (pnl, ID_MAT_LOAD, "Load"), 0,f2,2);
    sz0->Add (  new wxButton   (pnl, ID_MAT_SAVE, "Save"), 0,f2,2);
    sz0->Add (TxtFname, 0,f2,2);
    szt->Add          (sz0,0,f2,2);
    pnl->SetSizer     (szt);
    szt->Fit          (pnl);
    szt->SetSizeHints (pnl);

    // models
    Array<wxScrolled<wxPanel>*> P; // panels
    Array<wxFlexGridSizer*>     S; // sizers
    for (ModelFactory_t::iterator it=ModelFactory.begin(); it!=ModelFactory.end(); ++it)
    {
        CREATE_WXPANEL (p, s, 2, 2);
        P.Push (p);
        S.Push (s);
        Names.Push (it->first);
        if (Names.Last()=="CamClay")
        {
            WxDict * c = new WxDict (p);
            s->Add     (c, 0,f2,2);
            Dicts.Push (c);
            c->SetZero (-1, MODEL_PRM_NAMES["CamClay"]);
            c->ReBuild ();
        }
    }

    // notebook
    long nbk_style = (wxAUI_NB_DEFAULT_STYLE | wxAUI_NB_TAB_EXTERNAL_MOVE | wxNO_BORDER) & ~(wxAUI_NB_CLOSE_BUTTON | wxAUI_NB_CLOSE_ON_ACTIVE_TAB | wxAUI_NB_CLOSE_ON_ALL_TABS);
    wxAuiNotebook * nbk0 = new wxAuiNotebook(this, wxID_ANY, wxDefaultPosition, wxDefaultSize, nbk_style);
    for (size_t i=0; i<Names.Size(); ++i) nbk0->AddPage(P[i], Names[i], false);

    // set panes in Aui
    Aui.AddPane (pnl,  wxAuiPaneInfo().Name("cpnl").Caption("cpnl").Top().MinSize(wxSize(100,40)).DestroyOnClose(false).CaptionVisible(false) .CloseButton(false));
    Aui.AddPane (nbk0, wxAuiPaneInfo().Name("nbk0").Caption("nbk0").Centre().Position(0).DestroyOnClose(false).CaptionVisible(false).CloseButton(false));

    // commit all changes to wxAuiManager
    Aui.Update();
}

inline void MatFile::OnLoad (wxCommandEvent & Event)
{
    wxFileDialog fd(this, "Load material (.mat) file", LstDir, "", "*.mat");
    if (fd.ShowModal()==wxID_OK)
    {
        Read (fd.GetPath().ToStdString().c_str());
        TxtFname->SetValue (fd.GetFilename());
        LstDir = fd.GetDirectory ();
        _refresh ();
    }
}

inline void MatFile::OnSave (wxCommandEvent & Event)
{
    wxFileDialog fd(this, "Save material (.mat) file", LstDir, "", "*.mat", wxFD_SAVE|wxFD_OVERWRITE_PROMPT);
    if (fd.ShowModal()==wxID_OK)
    {
        std::fstream of(fd.GetPath().ToStdString().c_str(), std::ios::out);
        of.close();
    }
}

inline void MatFile::_refresh ()
{
}

#endif

#endif // MECHSYS_MATFILE_H
