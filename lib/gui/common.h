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

#ifndef MECHSYS_GUI_COMMON_H
#define MECHSYS_GUI_COMMON_H

// STL
#include <map>

#if defined(USE_FLTK)
  // FLTK
  #include <FL/Enumerations.H> // for Fl_Color
  #include <FL/fl_draw.H>
  #define GUI_DRAW_CIRCLE(x,y,r)     fl_circle (x, y, r);
  #define GUI_DRAW_FCIRCLE(x,y,r)    fl_pie    (x-r, y-r, 2*r+1, 2*r+1, 0, 360);
  #define GUI_DRAW_LINE(x0,y0,x1,y1) fl_line   (x0,y0, x1,y1);
  #define GUI_DRAW_TEXT(t,x,y,ralign)               \
    if (ralign) fl_draw (t,x,y,0,0,FL_ALIGN_RIGHT); \
    else        fl_draw (t,x,y);
#elif defined(USE_WXWIDGETS)
  // wxWidgets
  #include <wx/pen.h>
  #include <wx/gdicmn.h>
  #define GUI_DRAW_CIRCLE(x,y,r)     dc.DrawCircle (x, y, r);
  #define GUI_DRAW_FCIRCLE(x,y,r)    dc.DrawCircle (x, y, r);
  #define GUI_DRAW_LINE(x0,y0,x1,y1) dc.DrawLine   (x0,y0, x1,y1);
  #define GUI_DRAW_TEXT(t,x,y,ralign)    \
    wxCoord w,h;                         \
    dc.GetTextExtent (t,&w,&h);          \
    if (ralign) dc.DrawText (t,x-w,y-h); \
    else        dc.DrawText (t,x,  y-h);
#endif

namespace GUI
{

#if defined(USE_FLTK)
typedef std::map<String,Fl_Color> Clr_t;
typedef std::map<String,int>      Lty_t;
#elif defined(USE_WXWIDGETS)
typedef std::map<String,wxColour> Clr_t;
typedef std::map<String,int>      Lty_t;
#endif

#if defined(USE_WXWIDGETS)
wxColourDatabase ClrDb;
#endif

class Clr_c : public Clr_t
{
public:
    // Constructor
    Clr_c()
    {
#if defined(USE_FLTK)
        (*this)["black"]    = FL_BLACK;
        (*this)["blue"]     = FL_BLUE;
        (*this)["cyan"]     = FL_CYAN;
        (*this)["green"]    = FL_GREEN;
        (*this)["magenta"]  = FL_MAGENTA;
        (*this)["red"]      = FL_RED;
        (*this)["white"]    = FL_WHITE;
        (*this)["yellow"]   = FL_YELLOW;
        (*this)["dblue"]    = FL_DARK_BLUE;
        (*this)["dcyan"]    = FL_DARK_CYAN;
        (*this)["dgreen"]   = FL_DARK_GREEN;
        (*this)["dmagenta"] = FL_DARK_MAGENTA;
        (*this)["dred"]     = FL_DARK_RED;
        (*this)["dyellow"]  = FL_DARK_YELLOW;
#elif defined(USE_WXWIDGETS)
        if (ClrDb.Find("BLACK")      .Ok()) (*this)["black"]      = ClrDb.Find("BLACK");
        if (ClrDb.Find("BLUE")       .Ok()) (*this)["blue"]       = ClrDb.Find("BLUE");
        if (ClrDb.Find("CYAN")       .Ok()) (*this)["cyan"]       = ClrDb.Find("CYAN");
        if (ClrDb.Find("GREEN")      .Ok()) (*this)["green"]      = ClrDb.Find("GREEN");
        if (ClrDb.Find("MAGENTA")    .Ok()) (*this)["magenta"]    = ClrDb.Find("MAGENTA");
        if (ClrDb.Find("RED")        .Ok()) (*this)["red"]        = ClrDb.Find("RED");
        if (ClrDb.Find("WHITE")      .Ok()) (*this)["white"]      = ClrDb.Find("WHITE");
        if (ClrDb.Find("YELLOW")     .Ok()) (*this)["yellow"]     = ClrDb.Find("YELLOW");
        if (ClrDb.Find("BLUE VIOLET").Ok()) (*this)["blueviolet"] = ClrDb.Find("BLUE VIOLET");
        if (ClrDb.Find("DARK GREEN") .Ok()) (*this)["dgreen"]     = ClrDb.Find("DARK GREEN");
        if (ClrDb.Find("MAROON")     .Ok()) (*this)["maroon"]     = ClrDb.Find("MAROON");
        if (ClrDb.Find("FIREBRICK")  .Ok()) (*this)["firebrick"]  = ClrDb.Find("FIREBRICK");
        if (ClrDb.Find("ORANGE")     .Ok()) (*this)["orange"]     = ClrDb.Find("ORANGE");
#endif
    }

    // Operators
#if defined(USE_FLTK)
    Fl_Color operator() (char const * Name)
#elif defined(USE_WXWIDGETS)
    wxColour operator() (char const * Name)
#endif
    {
        Clr_t::const_iterator p = this->find(Name);
        if (p==this->end()) throw new Fatal("GUI::Clr(): Colour '%s' is not available",Name);
        return p->second;
    }
};

class Lty_c : public Lty_t
{
public:
    // Constructor
    Lty_c()
    {
#if defined(USE_FLTK)
        (*this)["solid"]   = FL_SOLID;
        (*this)["dash"]    = FL_DASH;
        (*this)["dot"]     = FL_DOT;
        (*this)["dashdot"] = FL_DASHDOT;
#elif defined(USE_WXWIDGETS)
        (*this)["solid"]   = wxSOLID;
        (*this)["dash"]    = wxSHORT_DASH;
        (*this)["dot"]     = wxDOT;
        (*this)["dashdot"] = wxDOT_DASH;
#endif
    }

    // Operators
#if defined(USE_FLTK)
    int operator() (char const * Name)
#elif defined(USE_WXWIDGETS)
    int operator() (char const * Name)
#endif
    {
        Lty_t::const_iterator p = this->find(Name);
        if (p==this->end()) throw new Fatal("GUI::Lty(): Linetype '%s' is not available",Name);
        return p->second;
    }
};

Clr_c Clr;
Lty_c Lty;

class Pen_c
{
public:
    // Constructor
    Pen_c (char const * ClrName="black", char const * LtyName="solid", int Lwd=1) { Set(ClrName,LtyName,Lwd); }

    // Methods
    void Set (char const * ClrName, char const * LtyName=NULL, int Lwd=1)
    {
#if defined(USE_FLTK)
        clr = Clr(ClrName);
        lty = Lty(LtyName);
        lwd = Lwd;
#elif defined(USE_WXWIDGETS)
        pen .SetColour (Clr(ClrName));
        pen .SetStyle  (Lty(LtyName));
        pen .SetWidth  (Lwd);
        spen.SetColour (Clr(ClrName));
        spen.SetStyle  (Lty("solid"));
        spen.SetWidth  (Lwd);
#endif
    }

#if defined(USE_FLTK)
    void Activate (bool ForceSolid=false)
    {
        fl_color (clr);
        if (ForceSolid) fl_line_style (FL_SOLID, 1);
    }
#elif defined(USE_WXWIDGETS)
    void Activate (wxDC & DC, bool ForceSolid=false)
    {
        if (ForceSolid) DC.SetPen (spen);
        else            DC.SetPen (pen);
    }
#endif

    // Data
#if defined(USE_FLTK)
    Fl_Color clr;
    int      lty;
    int      lwd;
#elif defined(USE_WXWIDGETS)
    wxPen    pen;
    wxPen    spen; // solid pen with the same color and width
#endif
};

}; // namespace GUI

#endif // MECHSYS_GUI_COMMON_H
