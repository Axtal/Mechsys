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

#ifndef MECHSYS_WXDICT_H
#define MECHSYS_WXDICT_H

// wxWidgets
#ifdef HAS_WXW
  #include <wx/grid.h>
  #include <mechsys/gui/common.h>
#endif

// MechSys
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

namespace GUI
{

class WxDictTable : public wxGridTableBase, public Dict
{
public:
    // Constructor
    WxDictTable () : ShowSubKeys(true), SameSubKeys(false) {}

    // Methods
    int      GetNumberRows ();
    int      GetNumberCols () { return Keys.Size(); }
    wxString GetValue      (int row, int col);
    void     SetValue      (int row, int col, wxString const & Str);
    bool     IsEmptyCell   (int row, int col) { return false; }

    // Formatting methods
    wxString GetColLabelValue (int col);
    wxString GetRowLabelValue (int row);

    // operators
    void operator= (Dict const & R); ///< Assignment operator

    // Data
    bool ShowSubKeys; ///< Show subkeys in cell ?
    bool SameSubKeys; ///< Show same subkeys for all columns ?
};

class WxDict : public wxWindow
{
public:
    // Constructor & Destructor
     WxDict (wxWindow * Parent);
    ~WxDict () { Aui.UnInit(); }

    // Methods
    void ReBuild (bool ReadControls=true); ///< Rebuild grid

    // Data
    wxAuiManager   Aui;      ///< Aui
    WxDictTable  * Tab;      ///< The table
    wxGrid       * Grd;      ///< The grid
    bool           FitCol;   ///< Fit data to columns
    bool           SameSK;   ///< Same subkeys
    bool           ShowSK;   ///< Show subkeys
    bool           HideCol0; ///< Hide column # 0

    // Events
    void OnReBuild (wxCommandEvent & Event) { ReBuild (); }
    DECLARE_EVENT_TABLE()
};


//////////////////////////////////////////////////////////////////////// Implementation ///// WxDictTable///////


inline int WxDictTable::GetNumberRows ()
{
    size_t nrows = 0;
    for (Dict_t::const_iterator it=begin(); it!=end(); ++it)
    {
        if (it->second.size()>nrows) nrows = it->second.size();
    }
    return nrows;
}

inline wxString WxDictTable::GetValue (int row, int col)
{
    SDPair const & pair = (*this)(Keys[col]);
    wxString buf;
    if (pair.Keys.Size()>(size_t)row)
    {
        if (ShowSubKeys) buf.Printf("%s: %g", pair.Keys[row].CStr(), pair(pair.Keys[row]));
        else             buf.Printf("%g",                            pair(pair.Keys[row]));
    }
    return buf;
}

inline void WxDictTable::SetValue (int row, int col, wxString const & Str)
{
    SDPair & pair = (*this)(Keys[col]);
    if (pair.Keys.Size()>(size_t)row)
    {
        double val;
        if (Str.ToDouble(&val)) pair(pair.Keys[row]) = val;
    }
}

inline wxString WxDictTable::GetColLabelValue (int col)
{
    if (Keys.Size()>(size_t)col)
    {
        wxString buf;  buf.Printf("%d", Keys[col]);
        return buf;
    }
    return wxEmptyString;
}

inline wxString WxDictTable::GetRowLabelValue (int row)
{
    if (SameSubKeys)
    {
        if (Keys.Size()>0 && begin()!=end())
        {
            SDPair const & pair = begin()->second;
            if (pair.Keys.Size()>(size_t)row) return wxString(pair.Keys[row]);
        }
        return wxEmptyString;
    }
    else return wxGridTableBase::GetRowLabelValue (row);
}

inline void WxDictTable::operator= (Dict const & R)
{
    clear();
    Keys.Resize (R.Keys.Size());
    for (size_t i=0; i<R.Keys.Size(); ++i)
    {
        Keys[i]          = R.Keys[i];
        (*this)[Keys[i]] = R(Keys[i]);
    }
}


//////////////////////////////////////////////////////////////////////// Implementation ////// WxDict //////////


enum
{
    ID_WXDICT_FITCOL = wxID_HIGHEST+3000,
    ID_WXDICT_SAMESK,
    ID_WXDICT_SHOWSK,
};

BEGIN_EVENT_TABLE (WxDict, wxWindow)
    EVT_CHECKBOX (ID_WXDICT_FITCOL, WxDict::OnReBuild)
    EVT_CHECKBOX (ID_WXDICT_SAMESK, WxDict::OnReBuild)
    EVT_CHECKBOX (ID_WXDICT_SHOWSK, WxDict::OnReBuild)
END_EVENT_TABLE ()

WxDict::WxDict (wxWindow * Parent)
    : wxWindow (Parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE),
      FitCol   (false),
      SameSK   (false),
      ShowSK   (true),
      HideCol0 (false)
{
    // force validation of child controls
    SetExtraStyle (wxWS_EX_VALIDATE_RECURSIVELY);

    // tell Aui to manage this frame
    Aui.SetManagedWindow (this);

    // control panel
    ADD_WXPANEL    (pnl, szt, szr, 1, 3);
    ADD_WXCHECKBOX (pnl, szr, ID_WXDICT_FITCOL, c0, "Fit columns to data", FitCol);
    ADD_WXCHECKBOX (pnl, szr, ID_WXDICT_SAMESK, c1, "Same subkeys",        SameSK);
    ADD_WXCHECKBOX (pnl, szr, ID_WXDICT_SHOWSK, c2, "Show subkeys",        ShowSK);

    // create table and grid
    Grd = new wxGrid           (this, wxID_ANY);
    Tab = new GUI::WxDictTable ();
    ReBuild ();

    // set Aui
    Aui.AddPane (pnl, wxAuiPaneInfo().Name("cpnl").Caption("Control").Top().MinSize(wxSize(100,30)).DestroyOnClose(false).CaptionVisible(false).CloseButton(false));
    Aui.AddPane (Grd, wxAuiPaneInfo().Name("grid").Caption("Grid").Centre().DestroyOnClose(false).CaptionVisible(false).CloseButton(false));
    Aui.Update  ();
}

void WxDict::ReBuild (bool ReadControls)
{
    // synchronise data from/to controls
    if (ReadControls) TransferDataFromWindow ();
    else              TransferDataToWindow   ();

    // rebuild grid based on updated table
    Grd->SetTable     (Tab, /*take_ownership*/false);
    Grd->ForceRefresh ();
    Tab->SameSubKeys = SameSK;
    Tab->ShowSubKeys = ShowSK;

    // reshape window
    if (FitCol) Grd->Fit ();
    Layout ();
    if (HideCol0 && Grd->GetNumberCols()>0) Grd->HideCol (0);
}

}; // namespace GUI

#endif // MECHSYS_WXDICT_H
