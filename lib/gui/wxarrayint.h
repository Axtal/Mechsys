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

#ifndef MECHSYS_WXARRAYINT_H
#define MECHSYS_WXARRAYINT_H

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

class WxArrayIntTable : public wxGridTableBase, public Array<int>
{
public:
    // Constructor
    WxArrayIntTable () : Transposed(false) {}

    // Derived Methods
    int      GetNumberRows () { return (Transposed ? 1 : Size()); }
    int      GetNumberCols () { return (Transposed ? Size() : 1); }
    wxString GetValue      (int row, int col);
    void     SetValue      (int row, int col, wxString const & Str);
    bool     IsEmptyCell   (int row, int col) { return false; }

    // Methods
    void DeleteCols (wxArrayInt const & ColsToDelete);
    void DeleteRows (wxArrayInt const & RowsToDelete);
    void Set        (Array<int> const & R);

    // Data
    bool Transposed; ///< Show values along lines instead columns
};

class WxArrayInt : public wxWindow
{
public:
    // Constructor
    WxArrayInt (wxWindow * Parent, WxArrayIntTable * tab=NULL);

    // Methods
    void ReBuild (bool ReadControls=true); ///< Rebuild grid

    // Data
    WxArrayIntTable * Tab;    ///< The table
    wxGrid          * Grd;    ///< The grid
    bool              FitCol; ///< Fit data to columns
    bool              Transp; ///< Tranposed

    // Events
    void OnReBuild (wxCommandEvent & Event) { ReBuild (); }
    void OnAdd     (wxCommandEvent & Event);
    void OnDel     (wxCommandEvent & Event);
    DECLARE_EVENT_TABLE()
};


//////////////////////////////////////////////////////////////////////// Implementation ///// WxArrayIntTable///////


inline wxString WxArrayIntTable::GetValue (int Row, int Col)
{
    int row = Row;
    int col = Col;
    if (Transposed) { row=Col; col=Row; }
    wxString buf;
    buf.Printf ("%d", (*this)[row]);
    return buf;
}

inline void WxArrayIntTable::SetValue (int Row, int Col, wxString const & Str)
{
    int row = Row;
    int col = Col;
    if (Transposed) { row=Col; col=Row; }
    long val;
    Str.ToLong (&val);
    (*this)[row] = val;
}

inline void WxArrayIntTable::DeleteCols (wxArrayInt const & ColsToDelete)
{
    if (!Transposed) return;
    Array<int> idxs_to_delete(ColsToDelete.size());
    for (size_t i=0; i<ColsToDelete.size(); ++i) idxs_to_delete[i] = ColsToDelete[i];
    DelItems (idxs_to_delete);
}

inline void WxArrayIntTable::DeleteRows (wxArrayInt const & RowsToDelete)
{
    if (Transposed) return;
    Array<int> idxs_to_delete(RowsToDelete.size());
    for (size_t i=0; i<RowsToDelete.size(); ++i) idxs_to_delete[i] = RowsToDelete[i];
    DelItems (idxs_to_delete);
}

inline void WxArrayIntTable::Set (Array<int> const & R)
{
    Resize (R.Size());
    for (size_t i=0; i<R.Size(); ++i) (*this)[i] = R[i];
}


//////////////////////////////////////////////////////////////////////// Implementation ////// WxArrayInt //////////


enum
{
    ID_WXARRAYINT_FITCOL = wxID_HIGHEST+3000,
    ID_WXARRAYINT_TRANSP,
    ID_WXARRAYINT_ADD,
    ID_WXARRAYINT_DEL,
};

BEGIN_EVENT_TABLE (WxArrayInt, wxWindow)
    EVT_CHECKBOX (ID_WXARRAYINT_FITCOL, WxArrayInt::OnReBuild)
    EVT_CHECKBOX (ID_WXARRAYINT_TRANSP, WxArrayInt::OnReBuild)
    EVT_BUTTON   (ID_WXARRAYINT_ADD,    WxArrayInt::OnAdd)
    EVT_BUTTON   (ID_WXARRAYINT_DEL,    WxArrayInt::OnDel)
END_EVENT_TABLE ()

WxArrayInt::WxArrayInt (wxWindow * Parent, WxArrayIntTable * tab)
    : wxWindow (Parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE),
      FitCol   (false),
      Transp   (false)
{
    // force validation of child controls
    SetExtraStyle (wxWS_EX_VALIDATE_RECURSIVELY);

    // controls
    wxBoxSizer * szt = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer * szr = new wxBoxSizer(wxHORIZONTAL);
    SetSizer          (szt);
    szt->Fit          (this);                                           
    szt->SetSizeHints (this);
    ADD_WXCHECKBOX    (this, szr, ID_WXARRAYINT_TRANSP, c0, "Transposed",          Transp);
    ADD_WXCHECKBOX    (this, szr, ID_WXARRAYINT_FITCOL, c1, "Fit columns to data", FitCol);
    ADD_WXBUTTON      (this, szr, ID_WXARRAYINT_ADD,    c2, "Add"                        );
    ADD_WXBUTTON      (this, szr, ID_WXARRAYINT_DEL,    c3, "Del"                        );
    szt->Add          (szr);

    // create table and grid
    Tab = (tab==NULL ? new GUI::WxArrayIntTable() : tab);
    Grd = new wxGrid (this, wxID_ANY);
    szt->Add (Grd, 0,wxALIGN_LEFT|wxALL|wxEXPAND,2);
    Transp = Tab->Transposed;
    ReBuild (false);
}

void WxArrayInt::ReBuild (bool ReadControls)
{
    // synchronise data from/to controls
    if (ReadControls) TransferDataFromWindow ();
    else              TransferDataToWindow   ();

    // rebuild grid based on updated table
    Tab->Transposed = Transp;
    Grd->SetTable     (Tab, /*take_ownership*/false);
    Grd->ForceRefresh ();

    // reshape window
    if (FitCol) Grd->Fit ();
    Layout ();
}

inline void WxArrayInt::OnAdd (wxCommandEvent & Event)
{
    Tab->Push (0);
    ReBuild   ();
}

inline void WxArrayInt::OnDel (wxCommandEvent & Event)
{
    wxArrayInt sel = (Tab->Transposed ? Grd->GetSelectedCols() : Grd->GetSelectedRows());
    if (Tab->Transposed) Tab->DeleteCols (sel);
    else                 Tab->DeleteRows (sel);
    ReBuild ();
}

}; // namespace GUI

#endif // MECHSYS_WXARRAYINT_H
