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

#ifndef MECHSYS_DATAPLT_H
#define MECHSYS_DATAPLT_H

// STL
#include <cmath>   // for ceil and floor
#include <cfloat>  // for DBL_EPSILON
#include <cstring> // for strncpy
#include <iostream>
#include <map>

// wxWidgets
#include <wx/frame.h>
#include <wx/aui/aui.h>
#include <wx/filedlg.h>
#include <wx/checkbox.h>
#include <wx/combobox.h>

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/array.h>
#include <mechsys/util/maps.h>
#include <mechsys/gui/plotxy.h>

namespace GUI
{

class DataPLT : public wxWindow
{
public:
    // Constructor
     DataPLT (wxFrame * Parent);
    ~DataPLT () { Aui.UnInit(); }

    // Methods
    void OnLoad    (wxCommandEvent & Event);
    void OnSelData (wxCommandEvent & Event) { _init_plots (); }
    void OnMultq   (wxCommandEvent & Event) { _init_plots (); }
    void OnPltAll  (wxCommandEvent & Event) { _init_plots (); }

    // Controls
    wxAuiManager   Aui;
    GUI::PlotXY  * qped, * qpev, * eved, * evlp;

    // Data
    wxCheckBox            * mult_q;               // multiply q according to Lode angle ?
    wxCheckBox            * plt_all;              // plot all data at the same time ?
    wxComboBox            * cbx_data;             // control with data filenames
    Array<Array<double> >   ed, ev;               // Octahedral invariants
    Array<Array<double> >   p, q, t, lp, qp, mqp; // Octahedral invariants

    // Events
    DECLARE_EVENT_TABLE()

private:
    void _init_plots ();
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


enum
{
    ID_LOAD = wxID_HIGHEST+1,
    ID_SELDATA,
    ID_MULTQ,
    ID_PLTALL,
};

BEGIN_EVENT_TABLE(DataPLT, wxWindow)
    EVT_BUTTON   (ID_LOAD,    DataPLT::OnLoad)
    EVT_COMBOBOX (ID_SELDATA, DataPLT::OnSelData)
    EVT_CHECKBOX (ID_MULTQ,   DataPLT::OnMultq)
    EVT_CHECKBOX (ID_PLTALL,  DataPLT::OnPltAll)
END_EVENT_TABLE()

inline DataPLT::DataPLT (wxFrame * Parent)
    : wxWindow (Parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE)
{
    // flags for sizers
    long f1 = wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL|wxALL;
    long f2 = wxALIGN_LEFT|wxALL|wxEXPAND;

    // tell wxAuiManager to manage this frame
    Aui.SetManagedWindow (this);

    // plots
    qped = new GUI::PlotXY (this, "Deviatoric stress/strains",             "ed [%]", "q/p");
    qpev = new GUI::PlotXY (this, "Deviatoric stress - volumetric strain", "ev [%]", "q/p");
    eved = new GUI::PlotXY (this, "Volumetric strain - deviatoric strain", "ed [%]", "ev [%]");
    evlp = new GUI::PlotXY (this, "Volumetric strain - natural log (p)",   "ln(p)",  "ev [%]");
    qped->ShowLastY = false;
    qpev->ShowLastY = false;
    eved->ShowLastY = false;
    evlp->ShowLastY = false;

    // control panel
    wxPanel         * pnl = new wxPanel    (this, wxID_ANY);
    wxBoxSizer      * szt = new wxBoxSizer (wxVERTICAL);
    wxFlexGridSizer * sz0 = new wxFlexGridSizer (/*rows*/1,/*cols*/4,/*vgap*/0,/*hgap*/0);
    cbx_data = new wxComboBox (pnl, ID_SELDATA, "Data files",wxDefaultPosition,wxSize(400,10));
    mult_q   = new wxCheckBox (pnl, ID_MULTQ,   "Lode multiply q");
    plt_all  = new wxCheckBox (pnl, ID_PLTALL,  "Plot all");
    sz0->Add ( new wxButton   (pnl, ID_LOAD,    "Load data"), 0,f2,2);
    sz0->Add (cbx_data, 0,f2,2);
    sz0->Add (mult_q,   0,f2,2);
    sz0->Add (plt_all,  0,f2,2);
    szt->Add (sz0, 0,f2,2);
    pnl->SetSizer     (szt);
    szt->Fit          (pnl);
    szt->SetSizeHints (pnl);
    mult_q->SetValue  (true);

    // panes
    Aui.AddPane (pnl,  wxAuiPaneInfo().Name("cpnl").Caption("Data").Top   ().MinSize(wxSize(100,40)) .DestroyOnClose(false).CaptionVisible(true) .CloseButton(false));
    Aui.AddPane (qped, wxAuiPaneInfo().Name("qped").Caption("qped").Centre().Position(0)             .DestroyOnClose(false).CaptionVisible(false).CloseButton(false));
    Aui.AddPane (eved, wxAuiPaneInfo().Name("eved").Caption("eved").Centre().Position(1)             .DestroyOnClose(false).CaptionVisible(false).CloseButton(false));
    Aui.AddPane (evlp, wxAuiPaneInfo().Name("evlp").Caption("evlp").Right ().MinSize(wxSize(300,100)).DestroyOnClose(false).CaptionVisible(false).CloseButton(false));
    Aui.AddPane (qpev, wxAuiPaneInfo().Name("qpev").Caption("qpev").Right ().MinSize(wxSize(300,100)).DestroyOnClose(false).CaptionVisible(false).CloseButton(false));

    // commit all changes to wxAuiManager
    Aui.Update();
}

inline void DataPLT::OnLoad (wxCommandEvent & Event)
{
    wxFileDialog fd(this, "Choose data file", "", "", "*.dat", wxFD_MULTIPLE);
    if (fd.ShowModal()==wxID_OK)
    {
        // get filenames
        wxArrayString paths, fnames;
        fd .GetPaths     (paths);
        fd .GetFilenames (fnames);

        // disconnect plots
        qped->DelCurves ();
        qpev->DelCurves ();
        eved->DelCurves ();
        evlp->DelCurves ();
        cbx_data->Clear ();

        // resize data arrays
        ed .Resize (paths.size());
        ev .Resize (paths.size());
        p  .Resize (paths.size());
        q  .Resize (paths.size());
        t  .Resize (paths.size());
        lp .Resize (paths.size());
        qp .Resize (paths.size());
        mqp.Resize (paths.size());

        // load data
        for (size_t i=0; i<paths.size(); ++i)
        {
            // read tabular data
            Table dat;
            dat.Read (paths[i]);
            bool   old_data_file_xyz = false;
            bool   old_data_file_art = false;
            bool   new_data_file_xyz = true;
            size_t nrow              = 0;
            bool   has_sxy           = false;
            for (size_t j=0; j<dat.Keys.Size(); ++j)
            {
                if      (dat.Keys[j]=="Sx")  { old_data_file_xyz=true;  nrow=dat("Sx").Size(); }
                else if (dat.Keys[j]=="Sa")  { old_data_file_art=true;  nrow=dat("Sa").Size(); }
                else if (dat.Keys[j]=="sx")  { new_data_file_xyz=true;  nrow=dat("sx").Size(); }
                else if (dat.Keys[j]=="Sxy") { has_sxy=true; }
                else if (dat.Keys[j]=="Sar") { has_sxy=true; }
                else if (dat.Keys[j]=="sxy") { has_sxy=true; }
            }
            if (nrow==0) throw new Fatal("DataPLT::OnLoad: Could not find (sx,sy,sz,sxy) columns in file %s",fnames[i].ToStdString().c_str());

            // set pointers to sub arrays of data
            Array<double> const * sx;
            Array<double> const * sy;
            Array<double> const * sz;
            Array<double> const * sxy = NULL;
            Array<double> const * ex;
            Array<double> const * ey;
            Array<double> const * ez;
            Array<double> const * exy = NULL;
            double m = 1.0; // multiplier for sig
            double n = 1.0; // multiplier for eps to convert to %
            if (old_data_file_xyz)
            {
                sx = &dat("Sx");
                sy = &dat("Sy");
                sz = &dat("Sz");
                ex = &dat("Ex");
                ey = &dat("Ey");
                ez = &dat("Ez");
                if (has_sxy) sxy = &dat("Sxy");
                if (has_sxy) exy = &dat("Exy");
                m = -1.0;
            }
            else if (old_data_file_art)
            {
                sz = &dat("Sa");
                sx = &dat("Sr");
                sy = &dat("St");
                ez = &dat("Ea");
                ex = &dat("Er");
                ey = &dat("Et");
                if (has_sxy) sxy = &dat("Srt");
                if (has_sxy) exy = &dat("Ert");
                m = -1.0;
            }
            else
            {
                sx = &dat("sx");
                sy = &dat("sy");
                sz = &dat("sz");
                ex = &dat("ex");
                ey = &dat("ey");
                ez = &dat("ez");
                if (has_sxy) sxy = &dat("sxy");
                if (has_sxy) exy = &dat("exy");
                n = 100.0;
            }

            // resize sub arrays
            ed [i].Resize (nrow);
            ev [i].Resize (nrow);
            p  [i].Resize (nrow);
            q  [i].Resize (nrow);
            t  [i].Resize (nrow);
            lp [i].Resize (nrow);
            qp [i].Resize (nrow);
            mqp[i].Resize (nrow);

            // calculate invariants
            Vec_t sig(4), eps(4);
            for (size_t j=0; j<nrow; ++j)
            {
                sig =   m*(*sx)[j],   m*(*sy)[j],   m*(*sz)[j], (has_sxy ?   m*(*sxy)[j] : 0.0);
                eps = n*m*(*ex)[j], n*m*(*ey)[j], n*m*(*ez)[j], (has_sxy ? n*m*(*exy)[j] : 0.0);
                OctInvs (sig, p[i][j], q[i][j], t[i][j]);
                ev [i][j] = Calc_evoct (eps);
                ed [i][j] = Calc_edoct (eps);
                lp [i][j] = log(p[i][j]);
                qp [i][j] = q[i][j]/p[i][j];
                mqp[i][j] = (t[i][j]<0.0 ? -qp[i][j] : qp[i][j]);
            }
        }

        // refresh cbx => replot
        cbx_data->Set          (fnames);
        cbx_data->SetSelection (0);
    }
}

inline void DataPLT::_init_plots ()
{
    // filenames
    wxArrayString fnames = cbx_data->GetStrings ();

    // disconnect plots
    qped->DelCurves ();
    qpev->DelCurves ();
    eved->DelCurves ();
    evlp->DelCurves ();

    // reconnect plots
    bool   mul = mult_q ->GetValue();
    bool   all = plt_all->GetValue();
    size_t ini = (all ? 0         : cbx_data->GetSelection());
    size_t num = (all ? ed.Size() : ini+1);
    for (size_t i=ini; i<num; ++i)
    {
        GUI::CurveProps & c0 = qped->AddCurve (&ed[i], (mul ? &mqp[i] : &qp[i]), fnames[i].ToStdString().c_str());
        GUI::CurveProps & c1 = qpev->AddCurve (&ev[i], (mul ? &mqp[i] : &qp[i]), fnames[i].ToStdString().c_str());
        GUI::CurveProps & c2 = eved->AddCurve (&ed[i],                  &ev[i],  fnames[i].ToStdString().c_str());
        GUI::CurveProps & c3 = evlp->AddCurve (&lp[i],                  &ev[i],  fnames[i].ToStdString().c_str());
        c0.Typ=GUI::CT_BOTH;  c0.Psz=4;  c0.Pen.Set((all ? GUI::LinClr(i) : "black"), (i>5 ? "dash" : "solid"), 1);
        c1.Typ=GUI::CT_BOTH;  c1.Psz=4;  c1.Pen.Set((all ? GUI::LinClr(i) : "black"), (i>5 ? "dash" : "solid"), 1);
        c2.Typ=GUI::CT_BOTH;  c2.Psz=4;  c2.Pen.Set((all ? GUI::LinClr(i) : "black"), (i>5 ? "dash" : "solid"), 1);
        c3.Typ=GUI::CT_BOTH;  c3.Psz=4;  c3.Pen.Set((all ? GUI::LinClr(i) : "black"), (i>5 ? "dash" : "solid"), 1);
    }

    // redraw
    qped->Redraw ();
    qpev->Redraw ();
    eved->Redraw ();
    evlp->Redraw ();
}

}; // namespace GUI

#endif // MECHSYS_DATAPLT_H
