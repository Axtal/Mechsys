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
  #include <wx/aui/aui.h>
  #include <wx/checkbox.h>
  #include <wx/valtext.h>
  #include <mechsys/gui/common.h>
  #include <mechsys/gui/wxfloatvalidator.h>
#endif

// MechSys
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>


using std::cout;
using std::endl;

namespace GUI
{

class WxDictTable : public wxGridTableBase, public Dict
{
public:
    // Constructor
    WxDictTable () : ShowSubKeys(true) {}

    // Methods
    int      GetNumberRows ();
    int      GetNumberCols () { return Keys.Size(); }
    wxString GetValue      (int row, int col);
    void     SetValue      (int row, int col, wxString const & Str);
    bool     IsEmptyCell   (int row, int col) { return false; }

    // Data
    bool ShowSubKeys; ///< Show subkeys in cell ?
};

class WxDict : public wxWindow
{
public:
    // static
    static const int WidFirstCol = 50;
    static const int WidCols     = 80;
    static const int HeiFirstRow = 20;
    static const int HeiRows     = 22;

    // Constructor & Destructor
     WxDict (wxWindow * Parent);
    ~WxDict () { Aui.UnInit(); }

    // Methods
    void Update  () { TransferDataFromWindow(); } ///< Update (validate/transfer) data in controls
    void ReBuild ();                              ///< Rebuild grid

    // Data
    wxAuiManager   Aui;    ///< Aui
    WxDictTable  * Tab;    ///< The table
    wxGrid       * Grd;    ///< The grid
    bool           FitCol; ///< Fit data to columns
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


//////////////////////////////////////////////////////////////////////// Implementation ////// WxDict //////////


enum
{
    ID_WXDICT_FITCOL = wxID_HIGHEST+3000,
};

WxDict::WxDict (wxWindow * Parent)
    : wxWindow (Parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE),
      FitCol   (false)
{
    // force validation of child controls
    SetExtraStyle (wxWS_EX_VALIDATE_RECURSIVELY);

    // tell Aui to manage this frame
    Aui.SetManagedWindow (this);

    // control panel
    ADD_WXPANEL    (pnl, szr, 1, 1)
    ADD_WXCHECKBOX (pnl, szr, ID_WXDICT_FITCOL, c0, "Fit columns to data", FitCol);

    // create table and grid
    Grd = new wxGrid           (this, wxID_ANY);
    Tab = new GUI::WxDictTable ();
    ReBuild ();

    // set Aui
    Aui.AddPane (pnl, wxAuiPaneInfo().Name("cpnl").Caption("Control").Top().MinSize(wxSize(100,40)).DestroyOnClose(false).CaptionVisible(false).CloseButton(false));
    Aui.AddPane (Grd, wxAuiPaneInfo().Name("grid").Caption("Grid").Centre().DestroyOnClose(false).CaptionVisible(false).CloseButton(false));
    Aui.Update  ();
}

void WxDict::ReBuild ()
{
    Update            ();
    Grd->SetTable     (Tab, /*take_ownership*/false);
    Grd->ForceRefresh ();
    if (FitCol) Grd->Fit ();

/*
    for (Dict_t::const_iterator it=begin(); it!=end(); ++it) if (it->second.size()>nrows) nrows = it->second.size();
    for (size_t i=0; i<Keys.Size(); ++i)
    {
        wxString key;  key.Printf("%d",Keys[i]);
        SetColLabelValue (i, key);
        SetColSize       (i, WidCols);
        if (!SamePairs)
        {
            throw new Fatal("WxDict::ReBuild: not yet (!SamePairs)");
            SDPair const & pair = (*this)(Keys[i]);
            for (size_t j=0; j<pair.Keys.Size(); ++j)
            {
            }
        }
    }
    if (SamePairs)
    {
        wxString buf;
        for (size_t i=0; i<begin()->second.Keys.Size(); ++i)
        {
            SetRowLabelValue (i, begin()->second.Keys[i]);
            SetRowSize       (i, HeiRows);
            SDPair const & pair = (*this)(Keys[i]);
            for (size_t j=0; j<ncols; ++j)
            {
                buf.Printf("%g",pair(pair.Keys[j]));
                SetCellValue (i,j,buf);
            }
        }
    }
*/
}

}; // namespace GUI

#endif // MECHSYS_WXDICT_H
