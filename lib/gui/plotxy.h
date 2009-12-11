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

#ifndef MECHSYS_PLOTXY_H
#define MECHSYS_PLOTXY_H

// STL
#include <cmath>   // for ceil and floor
#include <cfloat>  // for DBL_EPSILON
#include <cstring> // for strncpy
#include <iostream>
#include <map>

// FLTK
#ifdef USE_FLTK
  #include <FL/Fl.H>
  #include <FL/Fl_Group.H>
  #include <FL/Enumerations.H> // for Fl_Color
#endif

// wxWidgets
#ifdef USE_WXWIDGETS
  #include <wx/frame.h>
  #include <wx/dcbuffer.h>
#endif

// MechSys
#include "util/fatal.h"
#include "util/array.h"
#include "gui/common.h"

namespace GUI
{

enum CurveType { CT_POINTS, CT_LINES, CT_BOTH };

struct CurveProps
{
    CurveType Typ;      ///< Curve type
    Pen_c     Pen;      ///< Pen: color, lty, lwd
    int       Pch;      ///< Point type
    int       Psz;      ///< Point size
    bool      MaxY;     ///< Write Max Y ?
    char      Nam[256]; ///< Name
};

#if defined(USE_FLTK)
class PlotXY : public Fl_Group
#elif defined(USE_WXWIDGETS)
class PlotXY : public wxWindow
#else
#error MechSys:plotxy.h: Either USE_FLTK or USE_WXWIDGETS must be defined
#endif
{
public:
    /* Constructor. */
#if defined(USE_FLTK)
    PlotXY (int xmin, int ymin, int width, int height, char const * Title=NULL, char const * Xlbl=NULL, char const * Ylbl=NULL); // Screen coordinates
#elif defined(USE_WXWIDGETS)
    PlotXY (wxFrame * Parent, char const * Title=NULL, char const * Xlbl=NULL, char const * Ylbl=NULL);
#endif

    /* Destructor. */
    virtual ~PlotXY () {}

    // Methods
    CurveProps & AddCurve   (char const * Name);                                                   ///< Add curve (x-y values)
    CurveProps & AddCurve   (Array<double> const * X, Array<double> const * Y, char const * Name); ///< Add curve (x-y values)
    CurveProps & SetXY      (size_t i, Array<double> const * X, Array<double> const * Y);          ///< Set X and Y pointers
    void         DelCurve   (size_t i);                                                            ///< Remove curve (x-y values)
    void         CalcSF     ();                                                                    ///< Calculate scale factors
    void         DrawRulers ();                                                                    ///< Draw rulers
    void         DrawLegend ();                                                                    ///< Draw legend
    void         Redraw     ();                                                                    ///< Force redraw

    // Set methods
    PlotXY & EqScales     (bool EQScales=true) { _eqsf=EQScales;  return (*this); } ///< Set equal scale factors
    PlotXY & RecalcSF     (bool REcalcSF=true) { _recsf=REcalcSF; return (*this); } ///< Recalculate scale factors during draw ?
    PlotXY & WithFrame    (bool WFrame  =true) { _wframe=WFrame;  return (*this); } ///< Draw an all-around frame ?
    PlotXY & BRuler       (bool Visible =true, int nTicks=10);                      ///< Set bottom ruler
    PlotXY & LRuler       (bool Visible =true, int nTicks=10);                      ///< Set left ruler
    PlotXY & SetBTicksFmt (char const * Fmt="%4.1f", int Size=5);                   ///< Set bottom ruler ticks format, ex.: "%4.1f" => Size==5
    PlotXY & SetLTicksFmt (char const * Fmt="%3.1f", int Size=5);                   ///< Set left ruler ticks format, ex.: "%3.1f" => Size==5

    // Draw methods (internal)
#if defined(USE_FLTK)
    void draw ();
#elif defined(USE_WXWIDGETS)
    void OnPaint           (wxPaintEvent & Event);
    void OnEraseBackground (wxEraseEvent & Event) { /* Empty implementation to prevent flicker */ }
#endif

    // Data
    Array<CurveProps> C; ///< Curve properties

private:
    // Data
    double _sfX;        ///< Scale factor: x=xmin+(X-Xmin)*sf and y=ymin+ymax-(Y-Ymin)*sf, where x and y are screen coordinates
    double _sfY;        ///< Scale factor: x=xmin+(X-Xmin)*sf and y=ymin+ymax-(Y-Ymin)*sf, where x and y are screen coordinates
    bool   _eqsf;       ///< Equal scale factors ?
    bool   _recsf;      ///< Recalculate scale factors during draw ?
    bool   _sfok;       ///< Is scale factor ok (set true after a first time calculation) ?
    bool   _wframe;     ///< With frame? draw an all-around frame ?
    double _Xmin;       ///< Minimum x value (real coordinates)
    double _Ymin;       ///< Minimum y value (real coordinates)
    double _Xmax;       ///< Maximum x value (real coordinates)
    double _Ymax;       ///< Maximum y value (real coordinates)
    int    _lrt;        ///< Left ruler tickness (screen coordinates)
    int    _rrt;        ///< Right ruler tickness (screen coordinates)
    int    _brt;        ///< Bottom ruler tickness (screen coordinates)
    int    _trt;        ///< Top ruler tickness (screen coordinates)
    int    _hb;         ///< Increment for horizontal borders (to the inside) (screen coordinates)
    int    _vb;         ///< Increment for vertical borders (to the inside) (screen coordinates)
    int    _bnt;        ///< Bottom ruler number of ticks
    int    _lnt;        ///< Left ruler number of ticks
    int    _ticfsz;     ///< Ticks font size
    int    _lblfsz;     ///< Labels font size
    int    _titfsz;     ///< Title font size
    char   _btifmt[16]; ///< Bottom ruler ticks number format
    char   _ltifmt[16]; ///< Left ruler ticks number format
    char   _title[256]; ///< The title
    char   _blbl[64];   ///< Bottom ruler label
    char   _llbl[64];   ///< Left ruler lable
    int    _blegh;      ///< Bottom legend height
    int    _rlegw;      ///< Right legend width
    int    _bhpad;      ///< Border horizontal padding
    int    _bvpad;      ///< Border vertical padding
    int    _pahpad;     ///< Plot area horizontal padding
    int    _pavpad;     ///< Plot area vertical padding
    bool   _legatbot;   ///< Legend at bottom ?
    bool   _cptbr;      ///< Compact bottom ruler ?
#if defined(USE_WXWIDGETS)
    wxPen   _bg_pen;
    wxBrush _bg_brush;
    DECLARE_EVENT_TABLE()
#endif

    // Curves
    Array<Array<double> const*> _X; ///< X curves (real coordinates)
    Array<Array<double> const*> _Y; ///< Y curves (real coordinates)

#if defined(USE_WXWIDGETS)
    int x() const { return 0; }
    int y() const { return 0; }
    int w() const { int wi,he; GetClientSize(&wi,&he); return wi; }
    int h() const { int wi,he; GetClientSize(&wi,&he); return he; }
#endif

    // Private Methods
    int  _x      (double X) const { return static_cast<int>((x()+1+_lrt+_hb)              +_sfX*(X-_Xmin)); } ///< X in screen coordinates
    int  _y      (double Y) const { return static_cast<int>((y()+1+_trt+_vb)+(h()-_pavpad)-_sfY*(Y-_Ymin)); } ///< Y in screen coordinates
    int  _l      (double L) const { return static_cast<int>((_sfX>_sfY?_sfY:_sfX)*L); } ///< Length in screen coordinates
    void _pretty (double Lo, double Hi, int nDiv, Array<double> & Vals);                ///< Return a pretty serie of numbers between Lo and HI

}; // class PlotXY


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


#if defined(USE_FLTK)
inline PlotXY::PlotXY (int xmin, int ymin, int width, int height, char const * Title, char const * Xlbl, char const * Ylbl)
    : Fl_Group (xmin,ymin,width,height,0),
#elif defined(USE_WXWIDGETS)
inline PlotXY::PlotXY (wxFrame * Parent, char const * Title, char const * Xlbl, char const * Ylbl)
    : wxWindow (Parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE),
#endif
      _sfX      (1.0),
      _sfY      (1.0),
      _eqsf     (false),
      _recsf    (true),
      _sfok     (false),
      _wframe   (true),
      _Xmin     (0.0),
      _Ymin     (0.0),
      _Xmax     (1.0),
      _Ymax     (1.0),
      _hb       (6),
      _vb       (6),
      _bnt      (5),
      _lnt      (5),
      _ticfsz   (10),
      _lblfsz   (12),
      _titfsz   (14),
      _blegh    (16),
      _rlegw    (0),
      _legatbot (true),
      _cptbr    (true)
{
#ifdef USE_FLTK
    end();
#endif

    // Set dependent constants
    _lrt    = 40;                     // left ruler
    _rrt    = 0+5;                    // +5 to account for the width of the bottom tick text
    _brt    = (_cptbr?20:20+_ticfsz); // +_ticfsz to account for the label
    _trt    = 18;                     // area for the title
    _bhpad  = _lrt+_rrt+_rlegw;       // border horizontal padding
    _bvpad  = _brt+_trt+_blegh;       // border vertical padding
    _pahpad = _bhpad+2+2*_hb;         // plot-area horizontal padding
    _pavpad = _bvpad+2+2*_vb;         // plot-area vertical padding

    // Set tick number formats
    strncpy (_btifmt, "%g", 2);  _btifmt[2]='\0';
    strncpy (_ltifmt, "%g", 2);  _ltifmt[2]='\0';

    // Set title
    if (Title==NULL) { strncpy (_title, "X-Y plot", 8); _title[8]='\0'; }
    else               strncpy (_title, Title, 256);

    // Set x-label
    if (Xlbl==NULL) { strncpy (_blbl, "X", 1); _blbl[1]='\0'; }
    else              strncpy (_blbl, Xlbl, 64);

    // Set y-label
    if (Ylbl==NULL) { strncpy (_llbl, "Y", 1);  _llbl[1]='\0'; }
    else              strncpy (_llbl, Ylbl, 64);
}

inline CurveProps & PlotXY::AddCurve (char const * Name)
{
    // Default properties
    CurveProps cp = { CT_POINTS, Pen_c(), /*Pch*/1, /*Psz*/_hb*2, true };
    snprintf (cp.Nam, 256, "%s", Name);

    // Add curve
    _X.Push(NULL);
    _Y.Push(NULL);
     C.Push(cp);

    return C[C.Size()-1];
}

inline CurveProps & PlotXY::AddCurve (Array<double> const * X, Array<double> const * Y, char const * Name)
{
    // Check
    if (X!=NULL && Y!=NULL)
    {
        if (Y->Size()!=X->Size()) throw new Fatal("PlotXY::AddCurve: X (sz=%d) and Y (sz=%d) arrays must have the same size",X->Size(), Y->Size());
    }
    else throw new Fatal("PlotXY::AddCurve: X and Y must NOT be NULL");

    // Default properties
    CurveProps cp = { CT_POINTS, Pen_c(), /*Pch*/1, /*Psz*/_hb*2, true };
    snprintf (cp.Nam, 256, "%s", Name);

    // Add curve
    _X.Push(X);
    _Y.Push(Y);
     C.Push(cp);

    return C[C.Size()-1];
}

inline CurveProps & PlotXY::SetXY (size_t i, Array<double> const * X, Array<double> const * Y)
{
    // Check
    if (X!=NULL && Y!=NULL)
    {
        if (Y->Size()!=X->Size()) throw new Fatal("PlotXY::SetXY: X (sz=%d) and Y (sz=%d) arrays must have the same size",X->Size(), Y->Size());
    }
    else throw new Fatal("PlotXY::SetXY: X and Y must NOT be NULL");
    if (i<0 || i>=_X.Size()) throw new Fatal("PlotXY::SetXY: There is no Curve # %d added to this object",i);

    // Set curve
    _X[i] = X;
    _Y[i] = Y;

    return C[i];
}

inline void PlotXY::DelCurve (size_t i)
{
    // Check
    if (i<0 || i>=_X.Size()) throw new Fatal("PlotXY::DelCurve: There is no Curve # %d added to this object",i);

    // Remove
    _X.Remove(i);
    _Y.Remove(i);
     C.Remove(i);
}

inline void PlotXY::CalcSF ()
{
    _Xmin = 0.0;
    _Ymin = 0.0;
    _Xmax = 1.0;
    _Ymax = 1.0;
    size_t k = 0;
    while (k<_X.Size())
    {
        // Check input
        if (_X[k]==NULL)     break;
        if (_Y[k]==NULL)     break;
        if (_X[k]->Size()<2) break;
        if (_Y[k]->Size()<2) break;

        // Bounding box
        if (k==0)
        {
            _Xmin = (*_X[k])[0];
            _Ymin = (*_Y[k])[0];
            _Xmax = (*_X[k])[1];
            _Ymax = (*_Y[k])[1];
        }
        for (size_t i=0; i<_X[k]->Size(); ++i)
        {
            if ((*_X[k])[i]<_Xmin) _Xmin = (*_X[k])[i];
            if ((*_Y[k])[i]<_Ymin) _Ymin = (*_Y[k])[i];
            if ((*_X[k])[i]>_Xmax) _Xmax = (*_X[k])[i];
            if ((*_Y[k])[i]>_Ymax) _Ymax = (*_Y[k])[i];
        }

        // Next curve
        k++;
    }

    // Scale factors
    if (fabs(_Xmax-_Xmin)<=DBL_EPSILON)
    {
        Array<double> lim;
        _pretty (_Xmin, _Xmax, 3, lim);
        _Xmin = lim[0];
        _Xmax = lim[lim.Size()-1];
    }
    if (fabs(_Ymax-_Ymin)<=DBL_EPSILON)
    {
        Array<double> lim;
        _pretty (_Ymin, _Ymax, 3, lim);
        _Ymin = lim[0];
        _Ymax = lim[lim.Size()-1];
    }
    _sfX = static_cast<double>((w()-_pahpad)/(_Xmax-_Xmin));
    _sfY = static_cast<double>((h()-_pavpad)/(_Ymax-_Ymin));
    if (_eqsf)
    {
        double sf = (_sfX>_sfY ? _sfY : _sfX);
        _sfX = sf;
        _sfY = sf;
    }
}

inline void PlotXY::DrawRulers ()
{
#ifdef USE_FLTK
    if (_brt>0) // bottom ruler
    {
        // variables
        int yi = y()+h()-_brt-_blegh;  // initial y-position (screen coordinates)
        int X  = x()+_lrt;             // position
        int W  = w()-_rlegw-_lrt-_rrt; // width

        // background
        fl_color (FL_WHITE);
        fl_rectf (x(), yi, w(), _brt);

        // ticks
        const int     len = 9;                     // tick length
        Array<double> ticks;                       // ticks position
        char          buf[256];                    // buffer for ticks text
        _pretty       (_Xmin, _Xmax, _bnt, ticks); // find ticks position
        fl_color      (FL_BLACK);                  // set color
        fl_line_style (FL_SOLID, 1);               // set line style
        fl_font       (0,_ticfsz);                 // set font for ticks
        for (size_t i=0; i<ticks.Size(); ++i)
        {
            int xi = _x(ticks[i]);                  // tick initial x-position (screen coordinates)
            if (xi>=X && xi<=X+W)
            {
                snprintf (buf, 256, _btifmt, ticks[i]); // format text
                fl_line  (xi, yi, xi, yi+len);          // draw tick mark
                fl_draw  (buf, xi, yi+len+_ticfsz/2+1, 0, 0, FL_ALIGN_CENTER); // draw tick text
            }
        }

        // label 
        if (!_cptbr) fl_draw (_blbl, X+W/2, yi+len+3*_ticfsz/2+1, 0, 0, FL_ALIGN_CENTER); // draw label
    }
    if (_lrt>0) // left ruler
    {
        // variables
        int Y = y()+_trt;             // position
        int H = h()-_trt-_brt-_blegh; // height

        // background
        fl_color (FL_WHITE);
        fl_rectf (x(), Y, _lrt, H);

        // ticks
        const int     len = 9;                     // tick length
        Array<double> ticks;                       // ticks position
        char          buf[256];                    // buffer for ticks text
        _pretty       (_Ymin, _Ymax, _lnt, ticks); // find ticks position
        fl_color      (FL_BLACK);                  // set color
        fl_line_style (FL_SOLID, 1);               // set line style
        fl_font       (0,_ticfsz);                 // set font for ticks
        for (size_t i=0; i<ticks.Size(); ++i)
        {
            int xi = x()+_lrt-len;                  // tick initial x-position (screen coordinates)
            int yi = _y(ticks[i]);                  // tick initial y-position (screen coordinates)
            if (yi>=Y && yi<=Y+H)
            {
                snprintf (buf, 256, _ltifmt, ticks[i]); // format text
                fl_line  (xi, yi, xi+len, yi);          // draw tick mark
                fl_draw  (buf, x()+_lrt-len-2, yi, 0, 0, FL_ALIGN_RIGHT); // draw tick text
            }
        }
    }
    if (_trt>0) // top ruler
    {
        // background
        fl_color (FL_WHITE);
        fl_rectf (x(), y(), w(), _trt); // (x()+_lrt, y(), w()-_rrt-_lrt, _trt)

        // left ruler label and title
        fl_color (FL_BLACK);                 // set color
        fl_font  (0,_lblfsz);                // set font for label
        fl_draw  (_llbl, x()+2, y()+_trt-2); // draw left label
        fl_font  (0,_titfsz);                // set font for title
        fl_draw  (_title, x()+w()/2, y()+_trt-_titfsz/2-2, 0, 0, FL_ALIGN_CENTER); // draw title
    }
    if (_rrt>0) // right ruler
    {
        // background
        fl_color (FL_WHITE);
        fl_rectf (x()+w()-_rrt-_rlegw, y()+_trt, _rrt, h()-_trt-_brt-_blegh);
    }
#endif
}

inline void PlotXY::DrawLegend ()
{
#ifdef USE_FLTK
    // Variables
    char buf[256];
    int  xi   = x()+5;
    int  yi   = y()+h()-_blegh/2;
    int  len  = 30; // line length
    int  tlen = 30; // text width
    int  nc   = _X.Size(); // number of curves
    int  ilen = 5+len+2+tlen; // icon length
    int  dx   = (_cptbr ? 20 : static_cast<int>(static_cast<double>(w()-ilen*nc)/static_cast<double>(nc-1)));
    int  xf   = xi+30;

    // Draw legend
    fl_font (0,_lblfsz); // set font for labels
    if (_legatbot) // bottom legend
    {
        // background
        fl_color (FL_WHITE);
        fl_rectf (x(), y()+h()-_blegh, w(), _blegh);

        // Draw legend
        if (_X.Size()>=1 && _Y.Size()>=1)
        {
            size_t k = 0;
            while (k<_X.Size())
            {
                // Check
                if (_X[k]==NULL || _Y[k]==NULL) return;
                if (_X[k]->Size()!=_Y[k]->Size()) throw new Fatal("PlotXY::DrawLegend(): (%s) X(%s)[%d] and Y(%s)[%d] must have the same size",_title,_blbl,_X[k]->Size(),_llbl,_Y[k]->Size());

                // Draw points
                if (C[k].Typ==CT_POINTS || C[k].Typ==CT_BOTH)
                {
                    fl_color      (C[k].Pen.clr);
                    fl_line_style (FL_SOLID, 1);
                    if (C[k].Pch==1) // Open circle
                    {
                        fl_circle (xi+len/2, yi, C[k].Psz/2);
                    }
                    else if (C[k].Pch==16) // Filled circle
                    {
                        int r = C[k].Psz/2;
                        for (size_t i=0; i<_X[k]->Size(); ++i)
                            fl_pie (xi+len/2-r, yi-r, 2*r+1, 2*r+1, 0,360);
                    }
                }

                // Draw lines
                if (C[k].Typ==CT_LINES || C[k].Typ==CT_BOTH)
                {
                    if (_X[k]->Size()>1)
                    {
                        fl_color      (C[k].Pen.clr);
                        fl_line_style (C[k].Pen.lty, C[k].Pen.lwd);
                        fl_line       (xi, yi, xf, yi);
                    }
                }

                // Draw names
                snprintf (buf, 256, "%s", C[k].Nam); // format text
                fl_color (FL_BLACK);
                fl_draw  (buf, xf+2, yi, 0, 0, FL_ALIGN_LEFT); // draw name

                // Next curve
                xi += ilen+dx;
                xf  = xi+len;
                k++;
            }
        }
        if (_cptbr) fl_draw (_blbl, x()+w()-3, y()+h()-_ticfsz/2-3, 0, 0, FL_ALIGN_RIGHT); // draw label
    }
    else
    {
    }
#endif
}

inline void PlotXY::Redraw ()
{
#if defined(USE_FLTK)
    redraw();
    Fl::wait(0);
#elif defined(USE_WXWIDGETS)
    Refresh();
    Update();
#endif
}

inline PlotXY & PlotXY::BRuler (bool Visible, int nTicks)
{
    _brt = (Visible ? 27+12 : 0); // +12 to account for the label
    _bnt = nTicks;
    return (*this);
}

inline PlotXY & PlotXY::LRuler (bool Visible, int nTicks)
{
    _lrt = (Visible ? 27 : 0);
    _lnt = nTicks;
    return (*this);
}

inline PlotXY & PlotXY::SetBTicksFmt (char const * Fmt, int Size)
{
    strncpy (_btifmt, Fmt, Size);  _btifmt[Size]='\0';
    return (*this);
}

inline PlotXY & PlotXY::SetLTicksFmt (char const * Fmt, int Size)
{
    strncpy (_ltifmt, Fmt, Size);  _ltifmt[Size]='\0';
    return (*this);
}

#if defined(USE_FLTK)
inline void PlotXY::draw ()
#elif defined(USE_WXWIDGETS)
inline void PlotXY::OnPaint (wxPaintEvent & Event)
#endif
{
    // Clear the background
#if defined(USE_FLTK)
    fl_color (FL_WHITE);
    fl_rectf (x()+_lrt, y()+_trt, w()-_bhpad, h()-_bvpad);
#elif defined(USE_WXWIDGETS)
    wxBufferedPaintDC dc(this); // Device context
    PrepareDC (dc);
    dc.SetBrush      (_bg_brush);
    dc.SetPen        (_bg_pen);
    wxRect win_rect  (wxPoint(0,0), GetClientSize());
    dc.DrawRectangle (win_rect);
#endif

    // Calculate scale factors and draw rulers
    if (_recsf || !_sfok)
    {
        CalcSF     ();
        DrawRulers ();
        DrawLegend ();
        _sfok = true;
    }

    // Draw Curves
    if (_X.Size()>=1 && _Y.Size()>=1)
    {
        size_t k = 0;
        while (k<_X.Size())
        {
            // Check
            if (_X[k]==NULL || _Y[k]==NULL) break;

            // Draw points
            if (C[k].Typ==CT_POINTS || C[k].Typ==CT_BOTH)
            {
#if defined(USE_FLTK)
                C[k].Pen.Activate (/*ForceSolid*/true);
#elif defined(USE_WXWIDGETS)
                C[k].Pen.Activate (dc, /*ForceSolid*/true);
#endif
                if (C[k].Pch==1) // Open circle
                {
                    for (size_t i=0; i<_X[k]->Size(); ++i)
                        GUI_DRAW_CIRCLE (_x((*_X[k])[i]), _y((*_Y[k])[i]), C[k].Psz/2)
                }
                else if (C[k].Pch==16) // Filled circle
                {
                    for (size_t i=0; i<_X[k]->Size(); ++i)
                        GUI_DRAW_FCIRCLE (_x((*_X[k])[i]), _y((*_Y[k])[i]), C[k].Psz/2)
                }
            }

            // Draw lines
            if (C[k].Typ==CT_LINES || C[k].Typ==CT_BOTH)
            {
                if (_X[k]->Size()>1)
                {
#if defined(USE_FLTK)
                    C[k].Pen.Activate ();
#elif defined(USE_WXWIDGETS)
                    C[k].Pen.Activate (dc);
#endif
                    for (size_t i=0; i<_X[k]->Size()-1; ++i)
                        GUI_DRAW_LINE (_x((*_X[k])[i]), _y((*_Y[k])[i]), _x((*_X[k])[i+1]), _y((*_Y[k])[i+1]))
                }
            }

            // Draw text
            if (C[k].MaxY)
            {
                if (_X[k]->Size()>1)
                {
                    char buf[256];
                    //snprintf (buf, 250, "%g, %g", (*_X[k])[_X[k]->Size()-1], (*_Y[k])[_X[k]->Size()-1]);
                    snprintf (buf, 250, "y=%g", (*_Y[k])[_X[k]->Size()-1]);
#if defined(USE_FLTK)
                    fl_color (FL_BLACK);
                    fl_font  (0, 10);
#elif defined(USE_WXWIDGETS)
#endif
                    GUI_DRAW_TEXT (buf, _x((*_X[k])[_X[k]->Size()-1]), _y((*_Y[k])[_X[k]->Size()-1]), true);
                }
            }

            // Next curve
            k++;
        }
    }

#if defined(USE_FLTK)
    // Draw border
    fl_color      (FL_BLACK);
    fl_line_style (FL_SOLID, 1);
    fl_rect       (x()+_lrt, y()+_trt, w()-_bhpad, h()-_bvpad);
    // Draw an all-around frame
    if (_wframe) fl_rect (x(), y(), w(), h());
#endif
}

#if defined(USE_WXWIDGETS)
BEGIN_EVENT_TABLE(PlotXY, wxWindow)
    EVT_PAINT            (PlotXY::OnPaint)
    EVT_ERASE_BACKGROUND (PlotXY::OnEraseBackground)
END_EVENT_TABLE()
#endif

inline void PlotXY::_pretty (double Lo, double Hi, int nDiv, Array<double> & Vals)
{
    // Constants
    const double rounding_eps   = 1.0e-7;
    const double eps_correction = 0.0;
    const double shrink_sml     = 0.75;
    const double h              = 1.5;
    const double h5             = 0.5+1.5*h;

    // Local variables
    int    min_n = static_cast<int>(nDiv)/static_cast<int>(3);
    double lo    = Lo;
    double hi    = Hi;
    double dx    = hi-lo;
    double cell  = 1;    // cell := "scale" here
    double ub    = 0;    // upper bound on cell/unit
    bool   isml  = true; // is small ?

    // Check range
    if (!(dx==0 && hi==0)) // hi=lo=0
    {
        cell = (fabs(lo)>fabs(hi) ? fabs(lo) : fabs(hi));
        ub   = (1+(h5>=1.5*h+.5) ? 1/(1+h) : 1.5/(1+h5));
        isml = (dx<cell*ub*(nDiv>1 ? nDiv : 1)*DBL_EPSILON*3); // added times 3, as several calculations here
    }

    // Set cell
    if (isml)
    {
        if (cell>10) cell = 9 + cell/10;
        cell *= shrink_sml;
        if (min_n>1) cell /= min_n;
    }
    else 
    {
        cell = dx;
        if (nDiv>1) cell /= nDiv;
    }
    if      (cell<20*DBL_MIN) cell = 20*DBL_MIN; // very small range.. corrected
    else if (cell*10>DBL_MAX) cell = .1*DBL_MAX; // very large range.. corrected

    // Find base and unit
    double base = pow(10., floor(log10(cell))); // base <= cell < 10*base
    double unit = base;
    if ((ub = 2*base)-cell <  h*(cell-unit)) { unit = ub;
    if ((ub = 5*base)-cell < h5*(cell-unit)) { unit = ub;
    if ((ub =10*base)-cell <  h*(cell-unit))   unit = ub; }}

    // Find number of 
    double ns = floor (lo/unit+rounding_eps);
    double nu = ceil  (hi/unit-rounding_eps);
    if (eps_correction && (eps_correction>1 || !isml))
    {
        if (lo) lo *= (1-DBL_EPSILON); else lo = -DBL_MIN;
        if (hi) hi *= (1+DBL_EPSILON); else hi = +DBL_MIN;
    }
    while (ns*unit>lo+rounding_eps*unit) ns--;
    while (nu*unit<hi-rounding_eps*unit) nu++;

    // Find number of divisions
    int ndiv = static_cast<int>(.5+nu-ns);
    if (ndiv<min_n)
    {
        int k = min_n-ndiv;
        if (ns>=0.0)
        {
            nu += k/2;
            ns -= k/2 + k%2;
        } 
        else
        {
            ns -= k/2;
            nu += k/2 + k%2;
        }
        ndiv = min_n;
    }
    ndiv++;

    // Ensure that result covers original range
    if (ns*unit<lo) lo = ns * unit;
    if (nu*unit>hi) hi = nu * unit;

    // Fill array
    const double MINZERO = sqrt(DBL_EPSILON); // Minimum value to be replaced to 0.0 in _pretty method
    Vals.Resize(ndiv);
    Vals[0] = lo;
    for (int i=1; i<ndiv; ++i)
    {
        Vals[i] = Vals[i-1]+unit;
        if (fabs(Vals[i])<MINZERO) Vals[i] = 0.0;
    }
}

}; // namespace GUI

#endif // MECHSYS_PLOTXY_H
