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
    PenType   Pen;      ///< Pen: color, lty, lwd
    int       Pch;      ///< Point type
    int       Psz;      ///< Point size
    bool      LstY;     ///< Write Last Y ?
    char      Nam[256]; ///< Name
};

#if defined(USE_FLTK)
class PlotXY : public Fl_Group
#elif defined(USE_WXWIDGETS)
class PlotXY : public wxWindow
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
    void         DelCurves  ();                                                                    ///< Remove curves
    void         Redraw     ();                                                                    ///< Force redraw

    // Set methods
    PlotXY & SetBTicksFmt (char const * Fmt="%4.1f") { strcpy (_btifmt, Fmt); return (*this); } ///< Set bottom ruler ticks format
    PlotXY & SetLTicksFmt (char const * Fmt="%3.1f") { strcpy (_ltifmt, Fmt); return (*this); } ///< Set left ruler ticks format

    // Draw methods (internal)
    void CalcSF     ();                   ///< Calculate scale factors
    void DrawRulers (DeviceContext & DC); ///< Draw rulers
    void DrawLegend (DeviceContext & DC); ///< Draw legend
    void DrawCurves (DeviceContext & DC); ///< Draw curves
#if defined(USE_FLTK)
    void draw ();                         ///< Window redraw
#elif defined(USE_WXWIDGETS)
    void OnPaint           (wxPaintEvent & Event); ///< Window redraw
    void OnEraseBackground (wxEraseEvent & Event) { /* Empty implementation to prevent flicker */ }
#endif

    // Settings
    bool EqSF;      ///< Equal scale factors ?
    bool Grid;      ///< With primary grid ?
    bool RecSF;     ///< Recalculate scale factors during draw ?
    bool WFrame;    ///< With frame? draw an all-around frame ?
    int  TicFsz;    ///< Ticks font size
    int  LblFsz;    ///< Labels font size
    int  TitFsz;    ///< Title font size
    bool LegAtBot;  ///< Legend at bottom ?
    int  LegHei;    ///< Bottom legend height
    int  LegWid;    ///< Right legend width
    bool CompactBR; ///< Compact bottom ruler ?
    int  DefPsz;    ///< Default point size
    int  LRth;      ///< Left ruler tickness (screen coordinates)
    int  RRth;      ///< Right ruler tickness (screen coordinates)
    int  BRth;      ///< Bottom ruler tickness (screen coordinates)
    int  TRth;      ///< Top ruler tickness (screen coordinates)
    int  DelHB;     ///< Increment for horizontal borders (to the inside) (screen coordinates)
    int  DelVB;     ///< Increment for vertical borders (to the inside) (screen coordinates)
    int  TicLen;    ///< Length of tick line
    int  LinLen;    ///< Line length in legend
    int  LinDx;     ///< Spacing between legend items

    // Data
    Array<CurveProps> C;        ///< Curve properties
    FontType          LblFnt;   ///< Font for labels
    FontType          TitFnt;   ///< Font for title
    FontType          TicFnt;   ///< Font for ticks
    PenType           PLPen;    ///< Plot area pen
    PenType           BRPen;    ///< Bottom ruler pen
    PenType           TRPen;    ///< Top ruler pen
    PenType           LRPen;    ///< Left ruler pen
    PenType           RRPen;    ///< Righ ruler pen
    PenType           LEPen;    ///< Legend pen
    PenType           GrdPen;   ///< Grid pen
    Array<double>     BTicks;   ///< Bottom ticks
    Array<double>     LTicks;   ///< Left ticks

private:
    // Data
    double _sfX;        ///< Scale factor: x=xmin+(X-Xmin)*sf and y=ymin+ymax-(Y-Ymin)*sf, where x and y are screen coordinates
    double _sfY;        ///< Scale factor: x=xmin+(X-Xmin)*sf and y=ymin+ymax-(Y-Ymin)*sf, where x and y are screen coordinates
    bool   _sfok;       ///< Is scale factor ok (set true after a first time calculation) ?
    double _Xmin;       ///< Minimum x value (real coordinates)
    double _Ymin;       ///< Minimum y value (real coordinates)
    double _Xmax;       ///< Maximum x value (real coordinates)
    double _Ymax;       ///< Maximum y value (real coordinates)
    int    _bnt;        ///< Bottom ruler number of ticks
    int    _lnt;        ///< Left ruler number of ticks
    char   _btifmt[16]; ///< Bottom ruler ticks number format
    char   _ltifmt[16]; ///< Left ruler ticks number format
    char   _title[256]; ///< The title
    char   _blbl[64];   ///< Bottom ruler label
    char   _llbl[64];   ///< Left ruler lable
    int    _bhpad;      ///< Border horizontal padding
    int    _bvpad;      ///< Border vertical padding
    int    _pahpad;     ///< Plot area horizontal padding
    int    _pavpad;     ///< Plot area vertical padding

#if defined(USE_WXWIDGETS)
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
    int  _x      (double X) const { return static_cast<int>((x()+1+LRth+DelHB)              +_sfX*(X-_Xmin)); } ///< X in screen coordinates
    int  _y      (double Y) const { return static_cast<int>((y()+1+TRth+DelVB)+(h()-_pavpad)-_sfY*(Y-_Ymin)); } ///< Y in screen coordinates
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
      EqSF      (false),
      Grid      (true),
      RecSF     (true),
      WFrame    (true),
      TicFsz    (10),
      LblFsz    (12),
      TitFsz    (14),
      LegAtBot  (true),
      LegHei    (16),
      LegWid    (0),
      CompactBR (true),
      DefPsz    (8),
      LRth      (40),
      RRth      (5),
      BRth      (CompactBR ? 20 : 20+TicFsz),
      TRth      (18),
      DelHB     (6),
      DelVB     (6),
      TicLen    (9),
      LinLen    (40),
      LinDx     (20),
      _sfX      (1.0),
      _sfY      (1.0),
      _sfok     (false),
      _Xmin     (0.0),
      _Ymin     (0.0),
      _Xmax     (1.0),
      _Ymax     (1.0),
      _bnt      (10),
      _lnt      (10)
{
#ifdef USE_FLTK
    end();
#endif

    // Set dependent constants
    _bhpad  = LRth+RRth+(LegAtBot ? 0 : LegWid);
    _bvpad  = BRth+TRth+(LegAtBot ? LegHei : 0);
    _pahpad = _bhpad+2+2*DelHB;
    _pavpad = _bvpad+2+2*DelVB;

    // Set tick number formats
    strcpy (_btifmt, "%g");
    strcpy (_ltifmt, "%g");

    // Set title
    if (Title==NULL) { strncpy (_title, "X-Y plot", 8); _title[8]='\0'; }
    else               strncpy (_title, Title, 256);

    // Set x-label
    if (Xlbl==NULL) { strncpy (_blbl, "X", 1); _blbl[1]='\0'; }
    else              strncpy (_blbl, Xlbl, 64);

    // Set y-label
    if (Ylbl==NULL) { strncpy (_llbl, "Y", 1);  _llbl[1]='\0'; }
    else              strncpy (_llbl, Ylbl, 64);

    // Settings
    LblFnt.Set ("arial", LblFsz, false, false);
    TitFnt.Set ("arial", TitFsz, true,  false);
    TicFnt.Set ("arial", TicFsz, false, false);
    PLPen .Set (255,255,255);
    LEPen .Set (205,217,222);
    BRPen .Set (205,217,222);
    TRPen .Set (205,217,222);
    LRPen .Set (205,217,222);
    RRPen .Set (205,217,222);
    GrdPen.Set (140,140,140, "dash");
}

inline CurveProps & PlotXY::AddCurve (char const * Name)
{
    // Default properties
    CurveProps cp = { CT_POINTS, PenType(), /*Pch*/1, DefPsz, true };
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
    CurveProps cp = { CT_POINTS, PenType(), /*Pch*/1, DefPsz, true };
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

inline void PlotXY::DelCurves ()
{
    _X.Resize(0);
    _Y.Resize(0);
     C.Resize(0);
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
    if (EqSF)
    {
        double sf = (_sfX>_sfY ? _sfY : _sfX);
        _sfX = sf;
        _sfY = sf;
    }

    // ticks
    _pretty (_Xmin, _Xmax, _bnt, BTicks);
    _pretty (_Ymin, _Ymax, _lnt, LTicks);
}

inline void PlotXY::DrawRulers (DeviceContext & DC)
{
    if (BRth>0) // bottom ruler
    {
        // variables
        int yi = y()+h()-BRth-LegHei;  // initial y-position (screen coordinates)
        int X  = x()+LRth;             // position
        int W  = w()-LegWid-LRth-RRth; // width

        // background
        BRPen.Activate (DC);
        GUI_DRAW_FRECT (DC, x(), yi, w(), BRth);

        // ticks
        char buf[256]; // buffer for ticks text
        Black .Activate (DC);
        TicFnt.Activate (DC);
        for (size_t i=0; i<BTicks.Size(); ++i)
        {
            int xi = _x(BTicks[i]); // tick initial x-position (screen coordinates)
            if (xi>=X && xi<=X+W)
            {
                snprintf (buf, 256, _btifmt, BTicks[i]);                        // format text
                GUI_DRAW_LINE (DC, xi, yi, xi, yi+TicLen);                      // draw tick mark
                GUI_DRAW_TEXT (DC, buf, xi, yi+TicLen+TicFsz/2+1, /*center*/0); // draw tick text
            }
        }

        // label 
        if (!CompactBR) 
        {
            LblFnt.Activate (DC);
            GUI_DRAW_TEXT (DC, _blbl, X+W/2, yi+TicLen+3*TicFsz/2+1, /*center*/0); // draw label
        }
    }
    if (LRth>0) // left ruler
    {
        // variables
        int Y = y()+TRth;             // position
        int H = h()-TRth-BRth-LegHei; // height

        // background
        LRPen.Activate (DC);
        GUI_DRAW_FRECT (DC, x(), Y, LRth, H);

        // ticks
        char buf[256]; // buffer for ticks text
        Black .Activate (DC);
        TicFnt.Activate (DC);
        for (size_t i=0; i<LTicks.Size(); ++i)
        {
            int xi = x()+LRth-TicLen; // tick initial x-position (screen coordinates)
            int yi = _y(LTicks[i]);   // tick initial y-position (screen coordinates)
            if (yi>=Y && yi<=Y+H)
            {
                snprintf (buf, 256, _ltifmt, LTicks[i]); // format text
                GUI_DRAW_LINE (DC, xi, yi, xi+TicLen, yi); // draw tick mark
                GUI_DRAW_TEXT (DC, buf, x()+LRth-TicLen-2, yi, /*right*/1); // draw tick text
            }
        }
    }
    if (TRth>0) // top ruler
    {
        // background
        TRPen.Activate (DC);
        GUI_DRAW_FRECT (DC, x(), y(), w(), TRth); // (x()+LRth, y(), w()-RRth-LRth, TRth)

        // draw left label
        Black .Activate (DC);
        LblFnt.Activate (DC);
        GUI_DRAW_TEXT (DC, _llbl, x()+2, y()+TRth-6, /*left*/-1);

        // left ruler label and title
        TitFnt.Activate (DC);
        GUI_DRAW_TEXT (DC, _title, x()+w()/2, y()+TRth-TitFsz/2-2, /*center*/0); // draw title
    }
    if (RRth>0) // right ruler
    {
        // background
        RRPen.Activate (DC);
        GUI_DRAW_FRECT (DC, x()+w()-RRth-LegWid, y()+TRth, RRth, h()-TRth-BRth-LegHei);
    }
}

inline void PlotXY::DrawLegend (DeviceContext & DC)
{
    // Variables
    char buf[256];
    int  xi   = x()+5;
    int  yi   = y()+h()-LegHei/2;
    int  ilen = 5+LinLen+2; // icon length
    int  xf   = xi+30;

    // Draw legend
    if (LegAtBot) // bottom legend
    {
        // clear background
        LEPen.Activate (DC);
        GUI_DRAW_FRECT (DC, x(), y()+h()-LegHei, w(), LegHei);

        // Draw legend
        if (_X.Size()>=1 && _Y.Size()>=1)
        {
            for (size_t k=0; k<_X.Size(); ++k)
            {
                // Check
                if (_X[k]==NULL || _Y[k]==NULL) return;

                // Draw points
                if (C[k].Typ==CT_POINTS || C[k].Typ==CT_BOTH)
                {
                    C[k].Pen.Activate (DC, /*ForceSolid*/true);
                    if (C[k].Pch==1) // open circle
                    {
                        GUI_DRAW_CIRCLE (DC, xi+LinLen/2, yi, C[k].Psz/2)
                    }
                    else if (C[k].Pch==2) // filled circle
                    {
                        GUI_DRAW_FCIRCLE (DC, xi+LinLen/2, yi, C[k].Psz/2)
                    }
                    else if (C[k].Pch==3) // open square
                    {
                        GUI_DRAW_SQUARE (DC, xi+LinLen/2, yi, C[k].Psz)
                    }
                    else if (C[k].Pch==4) // filled square
                    {
                        GUI_DRAW_FSQUARE (DC, xi+LinLen/2, yi, C[k].Psz)
                    }
                }

                // Draw lines
                if (C[k].Typ==CT_LINES || C[k].Typ==CT_BOTH)
                {
                    if (_X[k]->Size()>1)
                    {
                        C[k].Pen.Activate (DC);
                        GUI_DRAW_LINE (DC, xi, yi, xf, yi);
                    }
                }

                // Draw names
                snprintf (buf, 256, "%s", C[k].Nam); // format text
                Black .Activate (DC);
                LblFnt.Activate (DC);
                GUI_DRAW_TEXT (DC, buf, xf+2, yi, /*left*/-1);

                // Next curve
                GUI_DRAW_ADD_WIDTH (DC,buf,xi)
                xi = xi+ilen+LinDx;
                xf = xi+LinLen;
            }
        }
        if (CompactBR)
        {
            LblFnt.Activate (DC);
            GUI_DRAW_TEXT (DC, _blbl, x()+w()-3, y()+h()-TicFsz/2-3, /*right*/1); // draw label
        }
    }
    else
    {
        // TODO
    }
}

inline void PlotXY::DrawCurves (DeviceContext & DC)
{
    // clear background
    PLPen.Activate (DC);
    GUI_DRAW_FRECT (DC, x()+LRth, y()+TRth, w()-_bhpad, h()-_bvpad);

    // draw grid
    if (Grid)
    {
        GrdPen.Activate (DC);
        for (size_t i=0; i<BTicks.Size(); ++i)
        {
            int xi = _x(BTicks[i]); // tick initial x-position (screen coordinates)
            GUI_DRAW_LINE (DC, xi, y()+TRth, xi, y()+h()-BRth-LegHei)
        }
        for (size_t i=0; i<LTicks.Size(); ++i)
        {
            int yi = _y(LTicks[i]); // tick initial y-position (screen coordinates)
            GUI_DRAW_LINE (DC, x()+LRth, yi, x()+w()-RRth, yi)
        }
    }

    // draw curves
    if (_X.Size()>=1 && _Y.Size()>=1)
    {
        for (size_t k=0; k<_X.Size(); ++k)
        {
            // Check
            if (_X[k]==NULL || _Y[k]==NULL) break;

            // Draw points
            if (C[k].Typ==CT_POINTS || C[k].Typ==CT_BOTH)
            {
                C[k].Pen.Activate (DC, /*ForceSolid*/true);
                if (C[k].Pch==1) // open circle
                {
                    for (size_t i=0; i<_X[k]->Size(); ++i)
                        GUI_DRAW_CIRCLE (DC, _x((*_X[k])[i]), _y((*_Y[k])[i]), C[k].Psz/2)
                }
                else if (C[k].Pch==2) // filled circle
                {
                    for (size_t i=0; i<_X[k]->Size(); ++i)
                        GUI_DRAW_FCIRCLE (DC, _x((*_X[k])[i]), _y((*_Y[k])[i]), C[k].Psz/2)
                }
                else if (C[k].Pch==3) // open square
                {
                    for (size_t i=0; i<_X[k]->Size(); ++i)
                        GUI_DRAW_SQUARE (DC, _x((*_X[k])[i]), _y((*_Y[k])[i]), C[k].Psz)
                }
                else if (C[k].Pch==4) // filled square
                {
                    for (size_t i=0; i<_X[k]->Size(); ++i)
                        GUI_DRAW_FSQUARE (DC, _x((*_X[k])[i]), _y((*_Y[k])[i]), C[k].Psz)
                }
            }

            // Draw lines
            if (C[k].Typ==CT_LINES || C[k].Typ==CT_BOTH)
            {
                if (_X[k]->Size()>1)
                {
                    C[k].Pen.Activate (DC);
                    for (size_t i=0; i<_X[k]->Size()-1; ++i)
                        GUI_DRAW_LINE (DC, _x((*_X[k])[i]), _y((*_Y[k])[i]), _x((*_X[k])[i+1]), _y((*_Y[k])[i+1]))
                }
            }

            // Draw text
            if (C[k].LstY)
            {
                if (_X[k]->Size()>1)
                {
                    char buf[256];
                    //snprintf (buf, 250, "%g, %g", (*_X[k])[_X[k]->Size()-1], (*_Y[k])[_X[k]->Size()-1]);
                    snprintf (buf, 250, "y=%g", (*_Y[k])[_X[k]->Size()-1]);
                    Black .Activate (DC);
                    TicFnt.Activate (DC);
                    GUI_DRAW_TEXT (DC, buf, _x((*_X[k])[_X[k]->Size()-1]), _y((*_Y[k])[_X[k]->Size()-1]), /*right*/1)
                }
            }
        }
    }
}

#if defined(USE_FLTK)
inline void PlotXY::draw ()
#elif defined(USE_WXWIDGETS)
inline void PlotXY::OnPaint (wxPaintEvent & Event)
#endif
{
    // get device context
#if defined(USE_FLTK)
    DeviceContext dc;
#elif defined(USE_WXWIDGETS)
    wxBufferedPaintDC dc(this);
    PrepareDC (dc);
#endif

    // Calculate scale factors and draw rulers
    if (RecSF || !_sfok)
    {
        CalcSF     ();
        DrawRulers (static_cast<DeviceContext&>(dc));
        DrawLegend (static_cast<DeviceContext&>(dc));
        _sfok = true;
    }

    // Draw Curves
    DrawCurves (static_cast<DeviceContext&>(dc));

    // Draw inner border
    Black.Activate (static_cast<DeviceContext&>(dc));
    GUI_DRAW_RECT  (dc, x()+LRth, y()+TRth, w()-_bhpad, h()-_bvpad);

    // Draw an all-around frame
    if (WFrame) GUI_DRAW_RECT (dc, x(), y(), w(), h());
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
    const double rounding_eps   = sqrt(DBL_EPSILON);
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
