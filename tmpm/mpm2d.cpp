/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso                                *
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

/* MPM2D - Copyright (C) 2007 Dorival de Moraes Pedroso */

// STL
#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat> // for DBL_EPSILON
#include <ctime>  // for std::clock()

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/util/fatal.h>

// Local
#include <mechsys/mpm/defs.h>
#include <mechsys/mpm/plotxy.h>
#include <mechsys/mpm/drawarea2d.h>
#include <mechsys/mpm/problems.h>
#include <mechsys/mpm/infobox.h>
#include <mechsys/mpm/tiled.h>
#include <mechsys/mpm/output.h>

// FLTK
#include <FL/Fl.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Round_Button.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Double_Window.H>

using std::cout;
using std::endl;
using MPM::Problem2D;
using MPM::DrawArea2D;
using MPM::Problems;
using MPM::DefaultProblem;
using MPM::PlotXY;
using MPM::MkdirP;
using MPM::CT_LINES;
using MPM::Vector3D;
using MPM::FourTiled;

/////////////////////////////////////////////////////////////////////////////////////////// TLframe /////

class TLframe : public Fl_Group
{
public:
    TLframe (int X, int Y, int W, int H, Problem2D const * Prob) : Fl_Group (X,Y,W,H,NULL), _pb(Prob)
    {
        int hp = (H-1)/2;
        _plt_f_u = new PlotXY (X+1         , Y    , (W-1)/2 , hp , "Force x Displacement"    , "Displacement" , "Force");
        _plt_s_e = new PlotXY (X+1+(W-1)/2 , Y    , (W-1)/2 , hp , "Stress x Strain"         , "Strain"       , "Stress");
        _plt_v_u = new PlotXY (X+1         , Y+hp , (W-1)/2 , hp , "Velocity x Displacement" , "Displacement" , "Velocity");
        _plt_f_t = new PlotXY (X+1+(W-1)/2 , Y+hp , (W-1)/2 , hp , "Force x Time"            , "Time"         , "Force");
        end();

        // Set Plot force x displacement
        _plt_f_u->AddCurve ("fx/ux");
        _plt_f_u->AddCurve ("fy/uy");
        _plt_f_u->SetCurve (0).Typ = CT_LINES;
        _plt_f_u->SetCurve (0).Clr = FL_BLUE;
        _plt_f_u->SetCurve (1).Typ = CT_LINES;
        _plt_f_u->SetCurve (1).Clr = FL_GREEN;

        // Set Plot stress x strain
        _plt_s_e->AddCurve ("exx/sxx");
        _plt_s_e->AddCurve ("eyy/syy");
        _plt_s_e->SetCurve (0).Typ = CT_LINES;
        _plt_s_e->SetCurve (0).Clr = FL_BLUE;
        _plt_s_e->SetCurve (1).Typ = CT_LINES;
        _plt_s_e->SetCurve (1).Clr = FL_GREEN;

        // Set Plot velocity x displacement
        _plt_v_u->AddCurve ("vx/ux");
        _plt_v_u->AddCurve ("vy/uy");
        _plt_v_u->SetCurve (0).Typ = CT_LINES;
        _plt_v_u->SetCurve (0).Clr = FL_BLUE;
        _plt_v_u->SetCurve (1).Typ = CT_LINES;
        _plt_v_u->SetCurve (1).Clr = FL_GREEN;

        // Set Plot displacement
        _plt_f_t->AddCurve ("fx");
        _plt_f_t->AddCurve ("fy");
        _plt_f_t->SetCurve (0).Typ = CT_LINES;
        _plt_f_t->SetCurve (0).Clr = FL_BLUE;
        _plt_f_t->SetCurve (1).Typ = CT_LINES;
        _plt_f_t->SetCurve (1).Clr = FL_GREEN;
    }
    void SetCurves ()
    {
        if (_pb->P->nPoints()>0 && _pb->Outp<_pb->P->nPoints())
        {
            _plt_f_u->SetXY (0, &_pb->Out_p_ux_t [_pb->Outp], &_pb->Out_p_fx_t [_pb->Outp]);
            _plt_f_u->SetXY (1, &_pb->Out_p_uy_t [_pb->Outp], &_pb->Out_p_fy_t [_pb->Outp]);
            _plt_s_e->SetXY (0, &_pb->Out_p_exx_t[_pb->Outp], &_pb->Out_p_sxx_t[_pb->Outp]);
            _plt_s_e->SetXY (1, &_pb->Out_p_eyy_t[_pb->Outp], &_pb->Out_p_syy_t[_pb->Outp]);
            _plt_v_u->SetXY (0, &_pb->Out_p_ux_t[_pb->Outp], &_pb->Out_p_vx_t[_pb->Outp]);
            _plt_v_u->SetXY (1, &_pb->Out_p_uy_t[_pb->Outp], &_pb->Out_p_vy_t[_pb->Outp]);
            _plt_f_t->SetXY (0, &_pb->Out_t, &_pb->Out_p_fx_t[_pb->Outp]);
            _plt_f_t->SetXY (1, &_pb->Out_t, &_pb->Out_p_fy_t[_pb->Outp]);
        }
    }
    void draw ()
    {
        _plt_f_u->draw();
        _plt_s_e->draw();
        _plt_v_u->draw();
        _plt_f_t->draw();
        fl_color (0xefebe700);
        fl_rectf (x()       , y()+h()-2 , w() , 3);
        fl_rectf (x()+w()-1 , y()       , 3   , h());
    }
private:
    PlotXY          * _plt_f_u; // plot force x displacement
    PlotXY          * _plt_s_e; // plot stress x strain
    PlotXY          * _plt_v_u; // plot velocity x displacement
    PlotXY          * _plt_f_t; // plot force x time
    Problem2D const * _pb;      // problem
}; // class TLframe

/////////////////////////////////////////////////////////////////////////////////////////// TRframe /////

class TRframe : public Fl_Group
{
public:
    TRframe (int X, int Y, int W, int H, Problem2D const * Prob) : Fl_Group (X,Y,W,H,NULL), _pb(Prob)
    {
        int hp = (H-1)/2;
        _plt_E_t = new PlotXY (X+1 , Y   , W-1 , hp , "Energy x Time"       , "Time" , "Energy");
        _plt_u_t = new PlotXY (X+1 , Y+hp, W-1 , hp , "Displacement x Time" , "Time" , "Displ.");
        end();

        // Set Plot energy
        _plt_E_t->AddCurve ("strain");
        _plt_E_t->AddCurve ("kinetic");
        _plt_E_t->AddCurve ("total");
        _plt_E_t->SetCurve (0).Typ = CT_LINES;
        _plt_E_t->SetCurve (0).Clr = FL_BLUE;
        _plt_E_t->SetCurve (1).Typ = CT_LINES;
        _plt_E_t->SetCurve (1).Clr = FL_GREEN;
        _plt_E_t->SetCurve (2).Typ = CT_LINES;
        _plt_E_t->SetCurve (2).Clr = FL_RED;

        // Set Plot displacement
        _plt_u_t->AddCurve ("ux");
        _plt_u_t->AddCurve ("uy");
        _plt_u_t->SetCurve (0).Typ = CT_LINES;
        _plt_u_t->SetCurve (0).Clr = FL_BLUE;
        _plt_u_t->SetCurve (1).Typ = CT_LINES;
        _plt_u_t->SetCurve (1).Clr = FL_GREEN;
    }
    void SetCurves ()
    {
        if (_pb->P->nPoints()>0 && _pb->Outp<_pb->P->nPoints())
        {
            _plt_E_t->SetXY (0, &_pb->Out_t, &_pb->Out_sE_t);
            _plt_E_t->SetXY (1, &_pb->Out_t, &_pb->Out_kE_t);
            _plt_E_t->SetXY (2, &_pb->Out_t, &_pb->Out_tE_t);
            _plt_u_t->SetXY (0, &_pb->Out_t, &_pb->Out_p_ux_t [_pb->Outp]);
            _plt_u_t->SetXY (1, &_pb->Out_t, &_pb->Out_p_uy_t [_pb->Outp]);
        }
    }
    void draw ()
    {
        _plt_E_t->draw();
        _plt_u_t->draw();
        fl_color (0xefebe700);
        fl_rectf (x(),y()+h()-2,w(),3);
    }
private:
    PlotXY          * _plt_E_t; // plot energy x time
    PlotXY          * _plt_u_t; // plot displacement x time
    Problem2D const * _pb;      // problem
}; // class TRframe

/////////////////////////////////////////////////////////////////////////////////////////// BRframe /////

class BRframe : public Fl_Group
{
public:
    BRframe (int X, int Y, int W, int H, Problem2D const * Prob) : Fl_Group (X,Y,W,H,NULL), _pb(Prob)
    {
        int hp = (H-1)/3;
        _plt_v_t = new PlotXY (X+1, Y+1     , W-1 , hp, "Velocity x Time", "Time", "Velocity");
        _plt_s_t = new PlotXY (X+1, Y+1+  hp, W-1 , hp, "Stress x Time",   "Time", "Stress"  );
        _plt_e_t = new PlotXY (X+1, Y+1+2*hp, W-1 , hp, "Strain x Time",   "Time", "Strain"  );
        end();

        // Set Plot velocity
        _plt_v_t->AddCurve ("vx");
        _plt_v_t->AddCurve ("vy");
        _plt_v_t->AddCurve ("cvx");
        _plt_v_t->SetCurve (0).Typ = CT_LINES;
        _plt_v_t->SetCurve (0).Clr = FL_BLUE;
        _plt_v_t->SetCurve (1).Typ = CT_LINES;
        _plt_v_t->SetCurve (1).Clr = FL_GREEN;
        _plt_v_t->SetCurve (2).Typ = CT_LINES;
        _plt_v_t->SetCurve (2).Clr = FL_RED;

        // Set Plot stress
        _plt_s_t->AddCurve ("sxx");
        _plt_s_t->AddCurve ("syy");
        _plt_s_t->SetCurve (0).Typ = CT_LINES;
        _plt_s_t->SetCurve (0).Clr = FL_BLUE;
        _plt_s_t->SetCurve (1).Typ = CT_LINES;
        _plt_s_t->SetCurve (1).Clr = FL_GREEN;

        // Set Plot strain
        _plt_e_t->AddCurve ("exx");
        _plt_e_t->AddCurve ("eyy");
        _plt_e_t->SetCurve (0).Typ = CT_LINES;
        _plt_e_t->SetCurve (0).Clr = FL_BLUE;
        _plt_e_t->SetCurve (1).Typ = CT_LINES;
        _plt_e_t->SetCurve (1).Clr = FL_GREEN;
    }
    void SetCurves ()
    {
        if (_pb->P->nPoints()>0 && _pb->Outp<_pb->P->nPoints())
        {
            _plt_v_t->SetXY (0, &_pb->Out_t, &_pb->Out_p_vx_t [_pb->Outp]);
            _plt_v_t->SetXY (1, &_pb->Out_t, &_pb->Out_p_vy_t [_pb->Outp]);
            _plt_v_t->SetXY (2, &_pb->Out_t, &_pb->Out_p_cvx_t[_pb->Outp]);
            _plt_s_t->SetXY (0, &_pb->Out_t, &_pb->Out_p_sxx_t[_pb->Outp]);
            _plt_s_t->SetXY (1, &_pb->Out_t, &_pb->Out_p_syy_t[_pb->Outp]);
            _plt_e_t->SetXY (0, &_pb->Out_t, &_pb->Out_p_exx_t[_pb->Outp]);
            _plt_e_t->SetXY (1, &_pb->Out_t, &_pb->Out_p_eyy_t[_pb->Outp]);
        }
    }
    void draw ()
    {
        _plt_v_t->draw();
        _plt_s_t->draw();
        _plt_e_t->draw();
    }
private:
    PlotXY          * _plt_v_t; // plot velocity x time
    PlotXY          * _plt_s_t; // plot stress x time
    PlotXY          * _plt_e_t; // plot strain x time
    Problem2D const * _pb;      // problem

}; // class BRframe

/////////////////////////////////////////////////////////////////////////////////////////// BLframe /////

class BLframe : public Fl_Group
{
public:
    BLframe (int X, int Y, int W, int H, Problem2D * Prob) : Fl_Group (X,Y,W,H,NULL), _pb(Prob),
     _tl(NULL), _tr(NULL), _br(NULL)
    {
        _da = new DrawArea2D (X+1,Y+1,W-2,H-1-3);
        end();

        // Set Draw area
        _da->RecalcSF (true);
        _da->SetTime  (&(Prob->t)); // time
        _da->SetM     (&(Prob->M)); // multiplier for external forces
        _da->callback (_execute_selection_cb, this); // method to be called when a point/node is selected in the DrawArea2D
    }
    void ReSetGridPointsAndOutp ()
    {
        _da->SetGrid   (_pb->G);
        _da->SetPoints (_pb->P);
        _da->ResetSelP (_pb->Outp);
    }
    void SetFrames (Fl_Group * TL, Fl_Group * TR, Fl_Group * BR)
    {
        _tl = static_cast<TLframe*>(TL);
        _tr = static_cast<TRframe*>(TR);
        _br = static_cast<BRframe*>(BR);
    }
    void draw ()
    {
        _da->damage (FL_DAMAGE_EXPOSE);
        _da->draw   ();
        fl_color    (0xefebe700);
        fl_rectf    (x()+w()-1 , y(), 3, h());
    }

    void SavePNG (char const * Filename) { _da->SavePNG (Filename); }

private:
    // Data
    DrawArea2D * _da;      // draw area
    Problem2D  * _pb;      // problem

    // Frames
    TLframe * _tl;
    TRframe * _tr;
    BRframe * _br;

    // Callbacks (must be static)
    static void _execute_selection_cb (Fl_Widget * o, void * v) { ((BLframe*)v)->_execute_selection(); }

    // Private methods
    void _execute_selection ()
    {
        _pb->Outp = _da->SelP(); // selected point
        if (_tl!=NULL) { _tl->SetCurves(); _tl->redraw(); }
        if (_tr!=NULL) { _tr->SetCurves(); _tr->redraw(); }
        if (_br!=NULL) { _br->SetCurves(); _br->redraw(); }
        Fl::wait(0);
    }

}; // class BLframe

///////////////////////////////////////////////////////////////////////////////////////// CallBacks /////

inline Fl_Group * AllocTL (int X, int Y, int W, int H, void * Prob) { return new TLframe (X,Y,W,H,static_cast<Problem2D*>(Prob)); }
inline Fl_Group * AllocTR (int X, int Y, int W, int H, void * Prob) { return new TRframe (X,Y,W,H,static_cast<Problem2D*>(Prob)); }
inline Fl_Group * AllocBL (int X, int Y, int W, int H, void * Prob) { return new BLframe (X,Y,W,H,static_cast<Problem2D*>(Prob)); }
inline Fl_Group * AllocBR (int X, int Y, int W, int H, void * Prob) { return new BRframe (X,Y,W,H,static_cast<Problem2D*>(Prob)); }

/////////////////////////////////////////////////////////////////////////////////////////// Control /////

Fl_Menu_Item NumberPointsPerCell[] =
{
    {"1" , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
    {"4" , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
    {"9" , 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
    {"16", 0, 0, 0, 0, FL_NORMAL_LABEL, 0, 0, 0},
    {0,0,0,0,0,0,0,0,0}
};

class Control : public Fl_Group
{
public:
    // Constructor
    Control (int X, int Y, int W, int H, Problem2D & Prob) : Fl_Group (X,Y,W,H,NULL), _prob(Prob),
       _running(false), _do_stop(false), _out_inc(0)
    {
        // Allocate widgets
        int xa[] = {50, 300, 350, 410, 470, 523, 605, 685};
        _w_prob    = new Fl_Choice      (xa[0], 2 , 200 , 25 , "Prob#");
        _w_npcell  = new Fl_Choice      (xa[1], 2 ,  50 , 25 , "NPcell");
        _w_cpdi    = new Fl_Check_Button(xa[2], 2 ,  30 , 25 , "CPDI");
        _w_mpm     = new Fl_Check_Button(xa[3], 2 ,  30 , 25 , "MPM");
        _w_usf     = new Fl_Check_Button(xa[4], 2 ,  30 , 25 , "USF");
        _w_smallst = new Fl_Check_Button(xa[5], 2 ,  30 , 25 , "SmallST");
        _w_fweuler = new Fl_Check_Button(xa[6], 2 ,  30 , 25 , "FwEuler");
        _w_sav_vtu = new Fl_Check_Button(xa[7], 2 ,  30 , 25 , "Save VTU");
        //_w_prob->box(FL_DOWN_BOX);

        // Widgets: next line
        int xb[] = {15, 105, 240, 304, 374, 444, 514, 584, 685};
        _w_tf      = new Fl_Input       (xb[0], 27 , 65 , 25 , "tf");
        _w_dt      = new Fl_Input       (xb[1], 27 , 80 , 25 , "Dt");
        _w_dtout   = new Fl_Input       (xb[2], 27 , 60 , 25 , "Out Dt");
        _w_setprob = new Fl_Button      (xb[3], 27 , 70 , 25 , "Set &Prob");
        _w_refine  = new Fl_Button      (xb[4], 27 , 70 , 25 , "Re&fine");
        _w_step    = new Fl_Button      (xb[5], 27 , 70 , 25 , "S&tep");
        _w_run     = new Fl_Button      (xb[6], 27 , 70 , 25 , "&Run");
        _w_close   = new Fl_Button      (xb[7], 27 , 70 , 25 , "&Close");
        _w_sav_png = new Fl_Check_Button(xb[8], 27 , 30 , 25 , "Save PNG");
        end();

        // Initialize widgets
        _w_prob  ->menu  (Problems);
        _w_prob  ->value (DefaultProblem-1);
        _w_npcell->menu  (NumberPointsPerCell);
        char buf[256];
        snprintf(buf, 256, "%g", _prob.tf    ); _w_tf     ->value(buf);
        snprintf(buf, 256, "%g", _prob.Dt    ); _w_dt     ->value(buf);
        snprintf(buf, 256, "%g", _prob.Dtout ); _w_dtout  ->value(buf);
        if (_prob.CPDI)                         _w_cpdi   ->set();
        if (_prob.MPM)                          _w_mpm    ->set();
        if (_prob.USF)                          _w_usf    ->set();
        if (_prob.SmallSt)                      _w_smallst->set();
        if (_prob.FwEuler)                      _w_fweuler->set();
        if (_prob.SaveVTU)                      _w_sav_vtu->set();
        if (_prob.SavePNG)                      _w_sav_png->set();

        // Callbacks
        _w_setprob->callback (_setprob_cb, this);
        _w_refine ->callback (_refine_cb , this);
        _w_step   ->callback (_step_cb   , this);
        _w_run    ->callback (_run_cb    , this);
        _w_close  ->callback (_close_cb  , this);

        // Set focus
        _w_run->take_focus();
    }
    
    // Methods
    void SetFrames (Fl_Group * TL, Fl_Group * TR, Fl_Group * BL, Fl_Group * BR)
    {
        _tl = static_cast<TLframe*>(TL);
        _tr = static_cast<TRframe*>(TR);
        _bl = static_cast<BLframe*>(BL);
        _br = static_cast<BRframe*>(BR);
    }

    // Destructor
    ~Control()
    {
        // Close files
        if (_main_file.is_open())
        {
            _main_file << "  </Collection>" << std::endl;
            _main_file << "</VTKFile>" << std::endl;
            _main_file.close ();
        }

        // Deallocate memory
        if (_prob.G!=NULL) delete _prob.G;
        if (_prob.P!=NULL) delete _prob.P;
    }

protected:
    // No resizable
    void resize(int X, int Y, int W, int H) {}

private:
    // Data
    Problem2D       & _prob;
    Fl_Choice       * _w_prob;
    Fl_Choice       * _w_npcell;
    Fl_Check_Button * _w_cpdi;
    Fl_Check_Button * _w_mpm;
    Fl_Check_Button * _w_usf;
    Fl_Check_Button * _w_smallst;
    Fl_Check_Button * _w_fweuler;
    Fl_Check_Button * _w_sav_vtu;
    Fl_Check_Button * _w_sav_png;
    Fl_Input        * _w_tf;
    Fl_Input        * _w_dt;
    Fl_Input        * _w_dtout;
    Fl_Button       * _w_setprob;
    Fl_Button       * _w_refine;
    Fl_Button       * _w_step;
    Fl_Button       * _w_run;
    Fl_Button       * _w_close;
    bool              _running;
    bool              _do_stop;
    char              _str_prob[64];
    char              _str_path[128];
    int               _out_inc; // output increment
    std::ofstream     _main_file;

    // Frames
    TLframe * _tl;
    TRframe * _tr;
    BLframe * _bl;
    BRframe * _br;

    // Callbacks (must be static)
    static void _setprob_cb (Fl_Widget * o, void * v) { ((Control*)v)->_setprob(); }
    static void _refine_cb  (Fl_Widget * o, void * v) { ((Control*)v)->_refine (); }
    static void _step_cb    (Fl_Widget * o, void * v) { ((Control*)v)->_step   (); }
    static void _run_cb     (Fl_Widget * o, void * v) { ((Control*)v)->_run    (); }
    static void _close_cb   (Fl_Widget * o, void * v) { ((Control*)v)->_close  (); }

    // Set problem
    void _setprob()
    {
        // Input
        _prob.Id     = _w_prob->value()+1;
        _prob.NPcell = static_cast<int>(pow(static_cast<double>(_w_npcell->value()+1.0),2.0));

        // Set problem
        SetProblem (_prob);

        std::cout << "Number of material points = " << _prob.P->nPoints() << std::endl;
        
        // Set GUI back
        char buf[256];
        snprintf(buf, 256, "%g", _prob.tf    ); _w_tf    ->value(buf);
        snprintf(buf, 256, "%g", _prob.Dt    ); _w_dt    ->value(buf);
        snprintf(buf, 256, "%g", _prob.Dtout ); _w_dtout ->value(buf);

        // Initial output
        _do_output (0.0);

        // Set draw area
        _bl->ReSetGridPointsAndOutp ();

        // Set plots
        _tl->SetCurves ();
        _tr->SetCurves ();
        _br->SetCurves ();

        // Redraw
        _do_redraw ();

        // set for output
        if (_prob.SavePNG || _prob.SaveVTU)
        {
            // Create directory
            snprintf (_str_prob, 64,  "p%d",        _prob.Id);
            snprintf (_str_path, 128, "results/%s", _str_prob);
            MkdirP   (_str_path);
            cout << "<" << _str_path << "> created." << endl;

            // Create main file
            if (_main_file.is_open())
            {
                _main_file << "  </Collection>" << std::endl;
                _main_file << "</VTKFile>" << std::endl;
                _main_file.close ();
            }
        }
        if (_prob.SaveVTU)
        {
            char str_main_file[256];
            snprintf (str_main_file, 256, "%s/main.pvd", _str_path);
            _main_file.open (str_main_file, std::ios::out);
            _main_file << "<?xml version=\"1.0\" ?>" << std::endl;
            _main_file << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
            _main_file << "  <Collection>" << std::endl;
        }

        // Save files
        _out_inc = 0;
        _save ();
        _out_inc++;
    }

    // Refine
    void _refine()
    {
        if (_prob.G!=NULL)
        {
            // Re-set data and refine
            _prob.t = 0.0;
            _prob.G->Refine   ();
            _prob.P->ReInit   ();
            ReSetOutputArrays (_prob);

            // Redraw
            _tl->SetCurves ();
            _tr->SetCurves ();
            _br->SetCurves ();
            _do_redraw     ();

            // Save files
            _out_inc = 0;
            _save ();
            _out_inc++;
        }
    }

    // Do one step
    bool _step()
    {
        if (_prob.G!=NULL && _prob.P!=NULL)
        {
            // Input
            _prob.tf      = atof(_w_tf   ->value());
            _prob.Dt      = atof(_w_dt   ->value());
            _prob.Dtout   = atof(_w_dtout->value());
            _prob.CPDI    = (_w_cpdi   ->value()==1 ? true : false);
            _prob.MPM     = (_w_mpm    ->value()==1 ? true : false);
            _prob.USF     = (_w_usf    ->value()==1 ? true : false);
            _prob.SmallSt = (_w_smallst->value()==1 ? true : false);
            _prob.FwEuler = (_w_fweuler->value()==1 ? true : false);
            _prob.SaveVTU = (_w_sav_vtu->value()==1 ? true : false);
            _prob.SavePNG = (_w_sav_png->value()==1 ? true : false);

            // Set MPoints
            _prob.P->SetCPDI (_prob.CPDI);
            _prob.P->SetMPM  (_prob.MPM);
            _prob.P->SetUSF  (_prob.USF);

            // Check
            if (fabs(_prob.t-_prob.tf)<DBL_EPSILON*10.0 || _prob.t>_prob.tf) return false;

            // Next t for output
            double to = _prob.t+_prob.Dtout;
            if (to>_prob.tf) to = _prob.tf;

            // Run
            double dse = 0.0;
            bool   ok  = true;
            try
            {
                while (_prob.t<to && ok)
                {
                    if (_prob.t+_prob.Dt>to) _prob.Dt = to-_prob.t;
                    ok = _prob.P->TimeUpdate(_prob.FwEuler, _prob.Dt, _prob.pB, _prob.pM,  _prob.t, dse);
                }
            }
            catch (Fatal * e)
            {
                e->Cout();
                delete e;
                ok = false;
            }

            // Do output
            _do_output (dse);

            // Redraw
            _do_redraw ();

            // Save files
            _save ();

            // Increment output index
            _out_inc++;

            // Next step
            return ok;
        }
        else return false;
    }

    // Run simulation
    void _run ()
    {
        if (_running) _do_stop = true;
        else
        {
            // Set variables
            _running = true;
            _w_run->label("&Stop");

            // Run
            std::clock_t start = std::clock();
            while (!_do_stop && _step());
            std::clock_t total = std::clock() - start;
            std::cout << "Time elapsed = [1;31m" << static_cast<double>(total)/CLOCKS_PER_SEC << "[0m seconds" << std::endl;

            // Reset variables
            _running = false;
            _do_stop = false;
            _w_run->label("&Run");

            // Save point data
            if (_prob.SavePointData)
            {
                std::ostringstream oss;
                oss << Util::_8s << "Time" << Util::_8s << "ux" << Util::_8s << "uy" << std::endl;
                for (size_t i=0; i<_prob.Out_t.Size(); ++i)
                {
                    oss << Util::_8s << _prob.Out_t[i] << Util::_8s << _prob.Out_p_ux_t[_prob.Outp][i] << Util::_8s << _prob.Out_p_uy_t[_prob.Outp][i] << std::endl;
                }
                char buf[256];
                snprintf (buf,  256, "%s/point_%04zd.res", _str_path, _prob.Outp);
                std::ofstream of(buf, std::ios::out);
                of << oss.str();
                of.close();
            }
        }
    }

    // Save
    void _save ()
    {
        char str_vtu[256];
        if (_prob.SaveVTU)
        {
            // Output to string
            size_t np = _prob.P->nPoints();
            std::ostringstream oss;
            oss << "<?xml version=\"1.0\"?>" << std::endl;
            oss << "<!-- Problem Time = " << _prob.t << " -->" << std::endl;
            oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
            oss << "  <UnstructuredGrid>" << std::endl;
            oss << "    <Piece NumberOfPoints=\"" << np << "\" NumberOfCells=\"" << np << "\">" << std::endl;
            oss << "      <Points>" << std::endl;
            oss << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
            int k = 0; oss << "        ";
            for (size_t i=0; i<np; ++i)
            {
                oss << "  " << _prob.P->P(i)(0) << " ";
                oss <<         _prob.P->P(i)(1) << " ";
                oss <<         _prob.P->P(i)(2);
                k++;
                if (k>7) { oss<<std::endl; if (i<np-1) oss<<"        "; k=0; }
                else { if (i==np-1) oss<<std::endl; }
            }
            oss << "        </DataArray>" << std::endl;
            oss << "      </Points>" << std::endl;
            oss << "      <Cells>" << std::endl;
            oss << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
            k = 0; oss << "        ";
            for (size_t i=0; i<np; ++i)
            {
                oss << "  " << i;
                k++;
                if (i==np-1) oss<<std::endl;
            }
            oss << "        </DataArray>" << std::endl;
            oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
            k = 0; oss << "        ";
            int offset = 0;
            for (size_t i=0; i<np; ++i)
            {
                offset += 1;
                oss << " " << offset;
                k++;
                if (i==np-1) oss<<std::endl;
            }
            oss << "        </DataArray>" << std::endl;
            oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
            k = 0; oss << "        ";
            for (size_t i=0; i<np; ++i)
            {
                oss << " " << "1";
                if (i==np-1) oss<<std::endl;
            }
            oss << "        </DataArray>" << std::endl;
            oss << "      </Cells>" << std::endl;
            oss << "      <PointData Scalars=\"TheScalars\" Vectors=\"TheVectors\">" << std::endl;
            MPM::VTU::Out::PointDataVectors    (_prob.P->u(),  "displacements", oss);
            MPM::VTU::Out::PointDataStressComp (_prob.P->Sta(), /*stress_xx*/0, oss);
            MPM::VTU::Out::PointDataStressComp (_prob.P->Sta(), /*stress_yy*/1, oss);
            MPM::VTU::Out::PointDataStressComp (_prob.P->Sta(), /*stress_xy*/3, oss);
            MPM::VTU::Out::PointDataStrainComp (_prob.P->Sta(), /*strain_xx*/0, oss);
            MPM::VTU::Out::PointDataStrainComp (_prob.P->Sta(), /*strain_yy*/1, oss);
            MPM::VTU::Out::PointDataStrainComp (_prob.P->Sta(), /*strain_xy*/3, oss);
            MPM::VTU::Out::PointDataSig        (_prob.P->Sta(), "sig",          oss);
            MPM::VTU::Out::PointDataIvs        (_prob.P->Sta(), "ivs",          oss);
            MPM::VTU::Out::PointDataTensors    (_prob.P->Fs(),  "F",            oss);
            oss << "      </PointData>" << std::endl;
            oss << "      <CellData Scalars=\"TheScalars\">" << std::endl;
            oss << "      </CellData>" << std::endl;
            oss << "    </Piece>" << std::endl;
            oss << "  </UnstructuredGrid>" << std::endl;
            oss << "</VTKFile>" << std::endl;

            // Save file
            snprintf (str_vtu,  256, "%s/res_%09d.vtu", _str_path, _out_inc);
            std::ofstream of(str_vtu, std::ios::out);
            of << oss.str();
            of.close();

            // Add vtu filename to main file
            if (!_main_file.is_open())
            {
                char str_main_file    [256];
                char str_main_file_tmp[256];
                snprintf (str_main_file,     256, "%s/main.pvd", _str_path);
                snprintf (str_main_file_tmp, 256, "%s/main_tmp.pvd", _str_path);
                std::ifstream in(str_main_file,     std::ios::in);
                std::ofstream tm(str_main_file_tmp, std::ios::out);
                std::string line;
                while (std::getline(in,line)) { if (line!="  </Collection>" && line!="</VTKFile>") tm << line << "\n"; }
                in.close();
                tm.close();
                std::remove (str_main_file);
                std::rename (str_main_file_tmp,str_main_file);
                _main_file.open (str_main_file, std::ios::app);
            }
            snprintf (str_vtu,  256, "res_%09d.vtu", _out_inc);
            _main_file << "    <DataSet timestep=\"" << _prob.t << "\" file=\"" << str_vtu << "\" />" << std::endl;
        }

        // Save PNG
        if (_prob.SavePNG)
        {
            snprintf (str_vtu,  256, "%s/fig_%09d.png", _str_path, _out_inc);
            _bl->SavePNG (str_vtu);
        }
    }

    // Close file
    void _close ()
    {
        // Close files
        if (_main_file.is_open())
        {
            _main_file << "  </Collection>" << std::endl;
            _main_file << "</VTKFile>" << std::endl;
            _main_file.close ();
        }
    }

    // Do output (fill arrays)
    void _do_output (double DsE)
    {
        // Check
        if (_prob.P->nPoints()<1) return;

        // Time
        _prob.Out_t.Push (_prob.t);

        // Strain energy
        size_t n = _prob.Out_sE_t.Size(); // current number of elements
        double se = 0.0;                  // strain energy
        if (n>0) se = _prob.Out_sE_t[n-1]+DsE;
        _prob.Out_sE_t.Push (se);

        // Kinetic energy
        double ke = 0.0; // kinetic energy
        for (size_t i=0; i<_prob.P->nPoints(); ++i)
            ke += 0.5*_prob.P->m(i)*blitz::dot(_prob.P->v(i),_prob.P->v(i));
        _prob.Out_kE_t.Push (ke);

        // Total energy
        _prob.Out_tE_t.Push (se+ke);

        // Mat point state values
        for (size_t p=0; p<_prob.P->nPoints(); ++p)
        {
            // Mat point external forces
            _prob.Out_p_fx_t[p].Push (_prob.P->f(p)(0));
            _prob.Out_p_fy_t[p].Push (_prob.P->f(p)(1));

            // Mat point displacement
            _prob.Out_p_ux_t[p].Push (_prob.P->u(p)(0));
            _prob.Out_p_uy_t[p].Push (_prob.P->u(p)(1));

            // Mat point velocity
            _prob.Out_p_vx_t[p].Push (_prob.P->v(p)(0));
            _prob.Out_p_vy_t[p].Push (_prob.P->v(p)(1));
            
            // Mat point correct velocity
            if (_prob.pVeloc!=NULL)
            {
                Vector3D cvp; // correct velocity
                (*_prob.pVeloc) (_prob.t, _prob.P->P(p), cvp);
                _prob.Out_p_cvx_t[p].Push (cvp(0));
            }
            
            // Mat point stress
            _prob.Out_p_sxx_t[p].Push (_prob.P->s(p)(0));
            _prob.Out_p_syy_t[p].Push (_prob.P->s(p)(1));

            // Mat point strain
            _prob.Out_p_exx_t[p].Push (_prob.P->e(p)(0));
            _prob.Out_p_eyy_t[p].Push (_prob.P->e(p)(1));
        }
    }

    // Do redraw
    void _do_redraw ()
    {
        _tl->redraw();
        _tr->redraw();
        _bl->redraw();
        _br->redraw();
        Fl::wait(0);
    }

}; // class Control

////////////////////////////////////////////////////////////////////////////////////////////// main /////

int main(int argc, char **argv) try
{
    // Set the backgroud color of all widgets
    Fl::set_color (FL_BACKGROUND_COLOR, 0xefebe700);

    // Initialize problem
    Problem2D prob;
    prob.Id      = 8;     // default problem
    prob.NPcell  = 4;     // number of points in each cell
    prob.G       = NULL;  // grid
    prob.P       = NULL;  // mat points
    prob.pB      = NULL;  // body mass
    prob.pVeloc  = NULL;  // correct velocity
    prob.pStress = NULL;  // correct stress
    prob.CPDI    = true;  // use CPDI
    prob.MPM     = false; // use original MPM
    prob.USF     = false; // update stress first ?
    prob.SmallSt = false; // small strains ?
    prob.FwEuler = true;  // Use Forward Euler ?
    prob.SaveVTU = false; // Save vtu at each output time ?
    prob.SavePNG = false; // Save png at each output time ?
    prob.t       = 0.0;   // current time
    prob.tf      = 1.0;   // final time
    prob.Dt      = 0.01;  // time increment
    prob.Dtout   = 0.1;   // time increment for output
    prob.Outp    = 0;     // point/particle number to output data
    prob.B       = 0.0;   // body mass

    // Allocate window, control and layout
    int ch = 56; // control heigh: 2x28
    Fl_Double_Window win (1800, 900);
    Control          ctr (0, 0, win.w(), ch, prob);
    FourTiled        lay (0,ch, win.w(), win.h()-ch, &prob, &AllocTL, &AllocTR, &AllocBL, &AllocBR, 0.67, 2./5.);

    // Set frames
    ctr.SetFrames (lay.TL(), lay.TR(), lay.BL(), lay.BR());
    static_cast<BLframe*>(lay.BL())->SetFrames (lay.TL(), lay.TR(), lay.BR());

    // Set and show window
    win.resizable (&lay);
    win.end       ();
    win.show      ();

    // Run
    Fl::run();

    // End
    return 0;
}
catch (Fatal * f)
{
    f->Cout();
    delete f;
    exit (1);
}
catch (char const * m)
{
    std::cout << "Fatal: " << m << std::endl;
    exit (1);
}
catch (...)
{
    std::cout << "Some exception (...) ocurred\n";
}
