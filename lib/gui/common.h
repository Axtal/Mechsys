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

// FLTK
#ifdef USE_FLTK
  #include <FL/Enumerations.H> // for Fl_Color
#endif

// wxWidgets
#ifdef USE_WXWIDGETS
  #include <wx/pen.h>
  #include <wx/gdicmn.h>
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
        (*this)["black"]    = (*wxBLACK);
        (*this)["blue"]     = (*wxBLUE);
        (*this)["cyan"]     = (*wxCYAN);
        (*this)["green"]    = (*wxGREEN);
        (*this)["red"]      = (*wxRED);
        (*this)["white"]    = (*wxWHITE);
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
        pen.SetColour (Clr(ClrName));
        pen.SetStyle  (Lty(LtyName));
        pen.SetWidth  (Lwd);
#endif
    }

    // Data
#if defined(USE_FLTK)
    Fl_Color clr;
    int      lty;
    int      lwd;
#elif defined(USE_WXWIDGETS)
    wxPen    pen;
#endif
};

}; // namespace GUI

#endif // MECHSYS_GUI_COMMON_H
