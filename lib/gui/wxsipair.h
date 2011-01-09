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

#ifndef MECHSYS_WXSIPAIR_H
#define MECHSYS_WXSIPAIR_H

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

class WxSIPairTable : public wxGridTableBase, public SIPair
{
public:
    // Constructor & Destructor
    WxSIPairTable () : Val2Name(NULL), Transposed(false) {}

    // Derived Methods
    int      GetNumberRows () { return (Transposed ? 1 : Keys.Size()); }
    int      GetNumberCols () { return (Transposed ? Keys.Size() : 1); }
    wxString GetValue      (int row, int col);
    void     SetValue      (int row, int col, wxString const & Str);
    bool     IsEmptyCell   (int row, int col) { return false; }

    // Methods
    void DeleteCols (wxArrayInt const & ColsToDelete);
    void DeleteRows (wxArrayInt const & RowsToDelete);

    // Data
    SDPair const * Val2Name;   ///< Maps 'name' value to description of subkey 'name'
    bool           Transposed; ///< Show keys along lines instead columns
};

class WxSIPair : public wxWindow
{
public:
    // Constructor & Destructor
     WxSIPair (wxWindow * Parent, WxSIPairTable * tab=NULL);

    // Methods
    void ReBuild (bool ReadControls=true); ///< Rebuild grid

    // Data
    WxSIPairTable * Tab;    ///< The table
    wxGrid        * Grd;    ///< The grid
    bool            FitCol; ///< Fit data to columns
    bool            Transp; ///< Tranposed

    // Events
    void OnReBuild (wxCommandEvent & Event) { ReBuild (); }
    void OnAdd     (wxCommandEvent & Event);
    void OnDel     (wxCommandEvent & Event);
    DECLARE_EVENT_TABLE()
};


//////////////////////////////////////////////////////////////////////// Implementation ///// WxSIPairTable///////


inline wxString WxSIPairTable::GetValue (int Row, int Col)
{
    int row = Row;
    int col = Col;
    if (Transposed) { row=Col; col=Row; }
    wxString buf;
    int val = (*this)(Keys[row]);
    if (Keys[row]=="name" && Val2Name!=NULL)
    {
        String name;
        Val2Name->Val2Key (val, name);
        buf = name;
    }
    else buf.Printf("%s: %d", Keys[row].CStr(), val);
    return buf;
}

inline void WxSIPairTable::SetValue (int Row, int Col, wxString const & Str)
{
    int row = Row;
    int col = Col;
    if (Transposed) { row=Col; col=Row; }
    if (Keys[row]=="name" && Val2Name!=NULL) return;
    long val;
    if (Str.ToLong(&val)) (*this)(Keys[row]) = val;
}

inline void WxSIPairTable::DeleteCols (wxArrayInt const & ColsToDelete)
{
    if (!Transposed) return;
    Array<String> keys_to_delete(ColsToDelete.size());
    for (size_t i=0; i<ColsToDelete  .size(); ++i) keys_to_delete[i] = Keys[ColsToDelete[i]];
    for (size_t i=0; i<keys_to_delete.size(); ++i) Del (keys_to_delete[i].CStr());
}

inline void WxSIPairTable::DeleteRows (wxArrayInt const & RowsToDelete)
{
    if (Transposed) return;
    Array<String> keys_to_delete(RowsToDelete.size());
    for (size_t i=0; i<RowsToDelete  .size(); ++i) keys_to_delete[i] = Keys[RowsToDelete[i]];
    for (size_t i=0; i<keys_to_delete.size(); ++i) Del (keys_to_delete[i].CStr());
}


//////////////////////////////////////////////////////////////////////// Implementation ////// WxSIPair //////////


enum
{
    ID_WXSIPAIR_FITCOL = wxID_HIGHEST+3000,
    ID_WXSIPAIR_TRANSP,
    ID_WXSIPAIR_ADD,
    ID_WXSIPAIR_DEL,
};

BEGIN_EVENT_TABLE (WxSIPair, wxWindow)
    EVT_CHECKBOX (ID_WXSIPAIR_FITCOL, WxSIPair::OnReBuild)
    EVT_CHECKBOX (ID_WXSIPAIR_TRANSP, WxSIPair::OnReBuild)
    EVT_BUTTON   (ID_WXSIPAIR_ADD,    WxSIPair::OnAdd)
    EVT_BUTTON   (ID_WXSIPAIR_DEL,    WxSIPair::OnDel)
END_EVENT_TABLE ()

WxSIPair::WxSIPair (wxWindow * Parent, WxSIPairTable * tab)
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
    ADD_WXCHECKBOX    (this, szr, ID_WXSIPAIR_TRANSP, c0, "Transposed",          Transp);
    ADD_WXCHECKBOX    (this, szr, ID_WXSIPAIR_FITCOL, c1, "Fit columns to data", FitCol);
    ADD_WXBUTTON      (this, szr, ID_WXSIPAIR_ADD,    c2, "Add"                        );
    ADD_WXBUTTON      (this, szr, ID_WXSIPAIR_DEL,    c3, "Del"                        );
    szt->Add          (szr);

    // create table and grid
    Tab = (tab==NULL ? new GUI::WxSIPairTable() : tab);
    Grd = new wxGrid (this, wxID_ANY);
    szt->Add (Grd, 0,wxALIGN_LEFT|wxALL|wxEXPAND,2);
    Transp = Tab->Transposed;
    ReBuild (false);
}

void WxSIPair::ReBuild (bool ReadControls)
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

inline void WxSIPair::OnAdd (wxCommandEvent & Event)
{
    String buf;
    for (size_t i=0; i<Tab->Keys.Size(); ++i)
    {
        int idx = atoi(Tab->Keys[i].CStr()) + 1;
        buf.Printf ("%d", idx);
        if (!Tab->HasKey(buf.CStr()))
        {
            Tab->Set (buf.CStr(), 0);
            break;
        }
    }
    if (Tab->Keys.Size()==0) Tab->Set ("0", 0);
    ReBuild ();
}

inline void WxSIPair::OnDel (wxCommandEvent & Event)
{
    wxArrayInt sel = (Tab->Transposed ? Grd->GetSelectedCols() : Grd->GetSelectedRows());
    if (Tab->Transposed) Tab->DeleteCols (sel);
    else                 Tab->DeleteRows (sel);
    ReBuild ();
}

}; // namespace GUI

#endif // MECHSYS_WXSIPAIR_H
