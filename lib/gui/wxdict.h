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

// Std Lib
#include <fstream>

// wxWidgets
#ifdef HAS_WXW
  #include <wx/grid.h>
  #include <wx/log.h>
  #include <mechsys/gui/common.h>
#endif

// MechSys
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

class WxDict : public wxGrid, public Dict
{
public:
    static const int WidFirstCol = 50;
    static const int WidCols     = 80;
    static const int HeiFirstRow = 20;
    static const int HeiRows     = 22;

    // Constructor
    WxDict (wxWindow * Parent);

    // Methods
    void ReBuild ();

    // Data
    bool SamePairs; ///< SDPairs contain the same keys

    // Events
    void OnCellValueChanging (wxGridEvent & Event);
    void OnCellChange        (wxGridEvent & Event);
    DECLARE_EVENT_TABLE();
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


BEGIN_EVENT_TABLE(WxDict, wxGrid)
    EVT_GRID_CELL_CHANGING (WxDict::OnCellValueChanging)
    EVT_GRID_CELL_CHANGED  (WxDict::OnCellChange)
END_EVENT_TABLE()

WxDict::WxDict (wxWindow * Parent)
    : wxGrid (Parent,wxID_ANY),
      SamePairs (true)
{
    CreateGrid       (1,1);
    SetRowLabelSize  (WidFirstCol);
    SetColLabelSize  (HeiFirstRow);
    //for (size_t i=0; i<3;  ++i) SetColSize (i,wid_cols);
    //for (size_t i=0; i<10; ++i) SetRowSize (i,hei_rows);
}

void WxDict::ReBuild ()
{
    DeleteRows(0,GetTable()->GetNumberRows());
    DeleteCols(0,GetTable()->GetNumberCols());
    int nrows = 0;
    int ncols = Keys.Size();
    for (Dict_t::const_iterator it=begin(); it!=end(); ++it) if (it->second.size()>nrows) nrows = it->second.size();
    AppendRows (nrows);
    AppendCols (ncols);
    for (size_t i=0; i<Keys.Size(); ++i)
    {
        wxString key;  key.Printf("%d",Keys[i]);
        SetColLabelValue (i, key);
        SetColSize       (i, WidCols);
        if (!SamePairs)
        {
            SDPair const & pair = (*this)(Keys[i]);
            for (size_t j=0; j<pair.Keys.Size(); ++j)
            {
            }
        }
    }
    if (SamePairs)
    {
        for (size_t i=0; i<begin()->second.Keys.Size(); ++i)
        {
            SetRowLabelValue (i, begin()->second.Keys[i]);
            SetRowSize       (i, HeiRows);
        }
    }
}


void WxDict::OnCellValueChanging (wxGridEvent & Event)
{
    int      row  = Event.GetRow();
    int      col  = Event.GetCol();
    wxString sval = GetCellValue(row,col);
    double   val;
    if (!sval.ToDouble(&val))
    {
        wxLogMessage ("Value=%s does not correspond to a real number",sval.ToStdString().c_str());
        Event.Veto();
        return;
    }
    Event.Skip();
}

void WxDict::OnCellChange (wxGridEvent & Event)
{
    int      row  = Event.GetRow();
    int      col  = Event.GetCol();
    wxString sval = GetCellValue(row,col);
    double   val;
    if (!sval.ToDouble(&val))
    { 
        wxLogError ("MyGrid::OnDatCellChange: Value=%s does not correspond to a real number",sval.ToStdString().c_str());
        return;
    }
}


#endif // MECHSYS_WXDICT_H
