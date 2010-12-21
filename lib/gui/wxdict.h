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
    // Constructor & Destructor
    WxDictTable () : Val2Name(NULL) {}

    // Derived Methods
    int      GetNumberRows ();
    int      GetNumberCols () { return Keys.Size(); }
    wxString GetValue      (int row, int col);
    void     SetValue      (int row, int col, wxString const & Str);
    bool     IsEmptyCell   (int row, int col) { return false; }

    // Derived formatting methods
    wxString GetColLabelValue (int col);

    // operators
    void operator= (Dict const & R); ///< Assignment operator

    // Methods
    void DeleteCols (wxArrayInt const & ColsToDelete);

    // Data
    SDPair const * Val2Name; ///< Maps 'name' value to description of subkey 'name'
};

class WxDict : public wxWindow
{
public:
    // Constructor & Destructor
     WxDict (wxWindow * Parent, WxDictTable * tab=NULL);

    // Methods
    void ReBuild (bool ReadControls=true); ///< Rebuild grid

    // Data
    WxDictTable * Tab;    ///< The table
    wxGrid      * Grd;    ///< The grid
    bool          FitCol; ///< Fit data to columns

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
        double val = pair(pair.Keys[row]);
        if (pair.Keys[row]=="name" && Val2Name!=NULL)
        {
            String name;
            Val2Name->Val2Key (val, name);
            buf = name;
        }
        else buf.Printf("%s: %g", pair.Keys[row].CStr(), val);
    }
    return buf;
}

inline void WxDictTable::SetValue (int row, int col, wxString const & Str)
{
    SDPair & pair = (*this)(Keys[col]);
    if (pair.Keys.Size()>(size_t)row)
    {
        if (pair.Keys[row]=="name" && Val2Name!=NULL) return;
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

inline void WxDictTable::DeleteCols (wxArrayInt const & ColsToDelete)
{
    wxArrayInt keys_to_delete(ColsToDelete.size());
    for (size_t i=0; i<ColsToDelete  .size(); ++i) keys_to_delete[i] = Keys[ColsToDelete[i]];
    for (size_t i=0; i<keys_to_delete.size(); ++i) Del(keys_to_delete[i]);
}


//////////////////////////////////////////////////////////////////////// Implementation ////// WxDict //////////


enum
{
    ID_WXDICT_FITCOL = wxID_HIGHEST+3000,
};

BEGIN_EVENT_TABLE (WxDict, wxWindow)
    EVT_CHECKBOX (ID_WXDICT_FITCOL, WxDict::OnReBuild)
END_EVENT_TABLE ()

WxDict::WxDict (wxWindow * Parent, WxDictTable * tab)
    : wxWindow (Parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE),
      FitCol   (false)
{
    // force validation of child controls
    SetExtraStyle (wxWS_EX_VALIDATE_RECURSIVELY);

    // controls
    wxBoxSizer * sz = new wxBoxSizer(wxVERTICAL);
    SetSizer         (sz);
    sz->Fit          (this);                                           
    sz->SetSizeHints (this);
    ADD_WXCHECKBOX   (this, sz, ID_WXDICT_FITCOL, c0, "Fit columns to data", FitCol);

    // create table and grid
    Tab = (tab==NULL ? new GUI::WxDictTable() : tab);
    Grd = new wxGrid (this, wxID_ANY);
    sz->Add (Grd, 0,wxALIGN_LEFT|wxALL|wxEXPAND,2);
    ReBuild ();
}

void WxDict::ReBuild (bool ReadControls)
{
    // synchronise data from/to controls
    if (ReadControls) TransferDataFromWindow ();
    else              TransferDataToWindow   ();

    // rebuild grid based on updated table
    Grd->SetTable     (Tab, /*take_ownership*/false);
    Grd->ForceRefresh ();

    // reshape window
    if (FitCol) Grd->Fit ();
    Layout ();
}

}; // namespace GUI

#endif // MECHSYS_WXDICT_H
