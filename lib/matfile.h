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
#ifdef USE_WXWIDGETS
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
    MatFile () {}
#endif

    // Methods
    void Read (char const * FileName);

    // Data
    Dict Prms; ///< Main data structure: maps tag to parameters

#ifdef USE_WXWIDGETS
    // Methods
    void Sync () { TransferDataFromWindow(); } ///< Synchronise (validate/transfer) data in controls

    // Data
    wxAuiManager          Aui;      ///< Aui manager
    wxString              LstDir;   ///< Last accessed directory
    wxTextCtrl          * TxtFName; ///< Control with input file filename
    String                FName;    ///< Material file (.mat) filename
    Array<String>         MdlNames; ///< Available model names
    Array<GUI::WxDict*>   Dicts;    ///< dictionaries with parameters

    // Events
    void OnLoad (wxCommandEvent & Event);
    void OnSave (wxCommandEvent & Event);
    DECLARE_EVENT_TABLE();
#endif
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void MatFile::Read (char const * FileName)
{
    // parse material file
    std::fstream mat_file(FileName, std::ios::in);
    if (!mat_file.is_open()) throw new Fatal("MatFile::Read: Could not open file <%s>",FileName);
}

std::ostream & operator<< (std::ostream & os, MatFile const & MF)
{
    os << MF.Prms << std::endl;
    return os;
}

#ifdef USE_WXWIDGETS

enum
{
    ID_MATFILE_LOAD = wxID_HIGHEST+2000,
    ID_MATFILE_SAVE ,
};

BEGIN_EVENT_TABLE(MatFile, wxWindow)
    EVT_BUTTON (ID_MATFILE_LOAD, MatFile::OnLoad)
    EVT_BUTTON (ID_MATFILE_SAVE, MatFile::OnSave)
END_EVENT_TABLE()

inline MatFile::MatFile (wxFrame * Parent)
    : wxWindow (Parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE)
{
    // force validation of child controls
    SetExtraStyle (wxWS_EX_VALIDATE_RECURSIVELY);

    // tell wxAuiManager to manage this window
    Aui.SetManagedWindow (this);

    // control panel
    ADD_WXPANEL     (pnl, szt, szr, 1, 3);
    ADD_WXBUTTON    (pnl, szr, ID_MATFILE_LOAD, c0, "Load");
    ADD_WXBUTTON    (pnl, szr, ID_MATFILE_SAVE, c1, "Save");
    ADD_WXTEXTCTRL_ (pnl, szr, wxID_ANY, TxtFName, "", FName);
    TxtFName->SetMinSize (wxSize(200,20));

    // models
    for (ModelFactory_t::iterator it=ModelFactory.begin(); it!=ModelFactory.end(); ++it)
    {
        MdlNames.Push (it->first);
        Dicts.Push    (new GUI::WxDict(this));
        Dicts.Last()->ShowSK = false;
        Dicts.Last()->SameSK = true;
        Str2ArrayStr_t::const_iterator mprms = MODEL_PRM_NAMES.find(MdlNames.Last());
        if (mprms!=MODEL_PRM_NAMES.end())
        {
            Dicts.Last()->Tab->SetZero (-1, mprms->second);
            Dicts.Last()->ReBuild      (false);
        }
        else WxError("MatFile::MatFile: __internal_error__ Model named <%s> is not in map: MODEL_PRM_NAMES",MdlNames.Last().CStr());
    }

    // notebook
    ADD_WXNOTEBOOK (this, nbk);
    for (size_t i=0; i<MdlNames.Size(); ++i) nbk->AddPage (Dicts[i], MdlNames[i], false);

    // commit all changes to wxAuiManager
    Aui.AddPane (pnl, wxAuiPaneInfo().Name("cpnl").Caption("cpnl").Top().MinSize(wxSize(100,40)).DestroyOnClose(false).CaptionVisible(false) .CloseButton(false));
    Aui.AddPane (nbk, wxAuiPaneInfo().Name("nbk0").Caption("nbk0").Centre().Position(0).DestroyOnClose(false).CaptionVisible(false).CloseButton(false));
    Aui.Update  ();
}

inline void MatFile::OnLoad (wxCommandEvent & Event)
{
    wxFileDialog fd(this, "Load material (.mat) file", LstDir, "", "*.mat");
    if (fd.ShowModal()==wxID_OK)
    {
        Read (fd.GetPath().ToStdString().c_str());
        TxtFName->SetValue (fd.GetFilename());
        LstDir = fd.GetDirectory ();
        for (size_t i=0; i<Dicts.Size(); ++i)
        {
            (*Dicts[i]->Tab) = Prms;
            Dicts[i]->ReBuild ();
        }
        TransferDataToWindow ();
    }
}

inline void MatFile::OnSave (wxCommandEvent & Event)
{
    Sync ();
    wxFileDialog fd(this, "Save material (.mat) file", LstDir, "", "*.mat", wxFD_SAVE|wxFD_OVERWRITE_PROMPT);
    if (fd.ShowModal()==wxID_OK)
    {
        std::fstream of(fd.GetPath().ToStdString().c_str(), std::ios::out);
        of.close();
    }
}

#endif

#endif // MECHSYS_MATFILE_H
