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
#include <sstream>
#include <cstring> // for strcmp

// FLTK or wxWidgets
#if defined(USE_FLTK)
  #include <FL/Enumerations.H>
  #include <FL/fl_draw.H>
#elif defined(USE_WXWIDGETS)
  #include <wx/pen.h>
  #include <wx/font.h>
  #include <wx/gdicmn.h>
  #include <wx/dcbuffer.h>
#else
  #error MechSys:gui/common.h: Either USE_FLTK or USE_WXWIDGETS must be defined
#endif

// MechSys
#include <mechsys/util/string.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/array.h>

#ifdef USE_WXWIDGETS

  #define CREATE_WXPANEL(PNL, SZR, M, N)                                      \
        wxScrolled<wxPanel> * PNL = new wxScrolled<wxPanel> (this, wxID_ANY); \
        wxFlexGridSizer * SZR = new wxFlexGridSizer (M,N,0,0);                \
        PNL->SetSizer (new wxBoxSizer(wxVERTICAL));                           \
        PNL->SetScrollbars(10,10,0,0);                                        \
        PNL->GetSizer()->Add (SZR,0,wxALIGN_LEFT|wxALL|wxEXPAND);             \
        PNL->GetSizer()->Fit (PNL);                                           \
        PNL->GetSizer()->SetSizeHints (PNL);

  #define ADD_WXREALNUMINPUT(PNL, SZR, ID, VAR, CTRL, LBL)                                             \
        CTRL = new WxRealNumInput (PNL, ID, VAR);                                                      \
        SZR->Add (new wxStaticText(PNL,wxID_ANY,LBL), 0,wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL|wxALL,2); \
        SZR->Add (CTRL, 0,wxALIGN_LEFT|wxALL|wxEXPAND,2);

  #define ADD_WXCHECKBOX(PNL, SZR, ID, VAR, CTRL, LBL)                                                 \
        CTRL = new wxCheckBox (PNL, ID, wxEmptyString);                                                \
        CTRL->SetValue (static_cast<bool>(VAR));                                                       \
        SZR->Add (new wxStaticText(PNL,wxID_ANY,LBL), 0,wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL|wxALL,2); \
        SZR->Add (CTRL, 0,wxALIGN_LEFT|wxALL|wxEXPAND,2);

  #define ADD_WXTEXTCTRL(PNL, SZR, ID, VAR, CTRL, LBL)                                                 \
        CTRL = new wxTextCtrl (PNL, ID, VAR);                                                          \
        SZR->Add (new wxStaticText(PNL,wxID_ANY,LBL), 0,wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL|wxALL,2); \
        SZR->Add (CTRL, 0,wxALIGN_LEFT|wxALL|wxEXPAND,2);

#endif

namespace GUI
{


//////////////////////////////////////////////////////////////////////////////////// Device Context /////


#if defined(USE_FLTK)
class DeviceContext
#elif defined(USE_WXWIDGETS)
class DeviceContext : public wxBufferedPaintDC
#endif
{};


/////////////////////////////////////////////////////////////////////////////////////////// Colours /////


#if defined(USE_FLTK)
typedef std::map<String,Fl_Color> Colours_t;
#elif defined(USE_WXWIDGETS)
typedef std::map<String,wxColour> Colours_t;
#endif

#if defined(USE_WXWIDGETS)
wxColourDatabase ClrDb;
#endif

class Colours : public Colours_t
{
public:
    // Constructor
    Colours()
    {
        // basic colours
#if defined(USE_FLTK)
        (*this)["black"]   = FL_BLACK;
        (*this)["blue"]    = FL_BLUE;
        (*this)["cyan"]    = FL_CYAN;
        (*this)["green"]   = FL_GREEN;
        (*this)["magenta"] = FL_MAGENTA;
        (*this)["red"]     = FL_RED;
        (*this)["white"]   = FL_WHITE;
        (*this)["yellow"]  = FL_YELLOW;
#elif defined(USE_WXWIDGETS)
        if (ClrDb.Find("BLACK")  .Ok()) (*this)["black"]   = ClrDb.Find("BLACK");
        if (ClrDb.Find("BLUE")   .Ok()) (*this)["blue"]    = ClrDb.Find("BLUE");
        if (ClrDb.Find("CYAN")   .Ok()) (*this)["cyan"]    = ClrDb.Find("CYAN");
        if (ClrDb.Find("GREEN")  .Ok()) (*this)["green"]   = ClrDb.Find("GREEN");
        if (ClrDb.Find("MAGENTA").Ok()) (*this)["magenta"] = ClrDb.Find("MAGENTA");
        if (ClrDb.Find("RED")    .Ok()) (*this)["red"]     = ClrDb.Find("RED");
        if (ClrDb.Find("WHITE")  .Ok()) (*this)["white"]   = ClrDb.Find("WHITE");
        if (ClrDb.Find("YELLOW") .Ok()) (*this)["yellow"]  = ClrDb.Find("YELLOW");
#endif
        // additional colours
        (*this)["pink"]    = (*this)(250, 204, 228);
        (*this)["lblue"]   = (*this)(217, 228, 255);
        (*this)["lgreen"]  = (*this)(100, 241, 193);
        (*this)["dblue"]   = (*this)( 45,   0, 160);
        (*this)["orange"]  = (*this)(241, 125,   0);
        (*this)["lyellow"] = (*this)(234, 228, 179);
        (*this)["dgreen"]  = (*this)(  0, 111,   5);
    }

    // Operators
#if defined(USE_FLTK)
    Fl_Color operator() (char const * Name)
#elif defined(USE_WXWIDGETS)
    wxColour operator() (char const * Name)
#endif
    {
        if (Name==NULL) throw new Fatal("GUI::Colours(): Name must not be NULL");
        Colours_t::const_iterator p = this->find(Name);
        if (p==this->end()) throw new Fatal("GUI::Colours(): Colour '%s' is not available",Name);
        return p->second;
    }

#if defined(USE_FLTK)
    Fl_Color operator() (unsigned char R, unsigned char G, unsigned char B, char const * Name=NULL)
    {
        if (Name==NULL) return fl_rgb_color (R,G,B);
        else
        {
            Colours_t::const_iterator p = this->find(Name);
            if (p==this->end()) // add new
            {
                Fl_Color clr = fl_rgb_color (R,G,B);
                (*this)[Name] = clr;
                return clr;
            }
            else return p->second;
        }
    }
#elif defined(USE_WXWIDGETS)
    wxColour operator() (unsigned char R, unsigned char G, unsigned char B, char const * Name=NULL)
    {
        if (Name==NULL) return wxColour(R,G,B);
        else
        {
            Colours_t::const_iterator p = this->find(Name);
            if (p==this->end()) // add new
            {
                wxColour clr(R,G,B);
                (*this)[Name] = clr;
                return clr;
            }
            else return p->second;
        }
    }
#endif

};

// Allocate colours
Colours Clr;

// Line colors
Array<String> __lin_clr;
int __init_lin_clr ()
{
    __lin_clr.Resize (10);
    __lin_clr = "red", "blue", "dgreen", "magenta", "dblue", "green", "orange", "cyan", "pink", "yellow";
    return 0;
}
int __dummy_init_lin_clr = __init_lin_clr();

inline char const * LinClr (int Idx=0) { return __lin_clr[Idx % __lin_clr.Size()].CStr(); }


//////////////////////////////////////////////////////////////////////////////////////// Line types /////


typedef std::map<String,int> LineTypes_t;

class LineTypes : public LineTypes_t
{
public:
    // Constructor
    LineTypes()
    {
#if defined(USE_FLTK)
        (*this)["solid"]   = FL_SOLID;
        (*this)["dash"]    = FL_DASH;
        (*this)["dot"]     = FL_DOT;
        (*this)["dashdot"] = FL_DASHDOT;
#elif defined(USE_WXWIDGETS)
        (*this)["solid"]   = wxSOLID;
        (*this)["dash"]    = wxLONG_DASH;
        (*this)["dot"]     = wxDOT;
        (*this)["dashdot"] = wxDOT_DASH;
#endif
    }

    // Operators
    int operator() (char const * Name)
    {
        if (Name==NULL) throw new Fatal("GUI::LineTypes(): Name must not be NULL");
        LineTypes_t::const_iterator p = this->find(Name);
        if (p==this->end()) throw new Fatal("GUI::LineTypes(): Linetype '%s' is not available",Name);
        return p->second;
    }
};

// Allocate linetypes
LineTypes Lty;


/////////////////////////////////////////////////////////////////////////////////////////// PenType /////


class PenType
{
public:
    // Constructors
    PenType (char const * ClrName="black", char const * LtyName="solid", int Lwd=1) { Set(ClrName,LtyName,Lwd); }
    PenType (unsigned char R, unsigned char G, unsigned char B, char const * LtyName="solid", int Lwd=1) { Set(R,G,B,  LtyName,Lwd); }

    // Methods
    void Set (char const * ClrName, char const * LtyName="solid", int Lwd=1)
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
        spen.SetWidth  (1);
#endif
    }

    void Set (unsigned char R, unsigned char G, unsigned char B, char const * LtyName="solid", int Lwd=1)
    {
#if defined(USE_FLTK)
        clr = Clr(R,G,B);
        lty = Lty(LtyName);
        lwd = Lwd;
#elif defined(USE_WXWIDGETS)
        pen .SetColour (Clr(R,G,B));
        pen .SetStyle  (Lty(LtyName));
        pen .SetWidth  (Lwd);
        spen.SetColour (Clr(R,G,B));
        spen.SetStyle  (Lty("solid"));
        spen.SetWidth  (1);
#endif
    }

    void Activate (DeviceContext & DC, bool ForceSolid=false)
    {
#if defined(USE_FLTK)
        fl_color (clr);
        if (ForceSolid) fl_line_style (FL_SOLID, 1);
        else            fl_line_style (lty, lwd);
#elif defined(USE_WXWIDGETS)
        if (ForceSolid) DC.SetPen (spen);
        else            DC.SetPen (pen);
#endif
    }

    // Data
#if defined(USE_FLTK)
    Fl_Color clr;
    int      lty;
    int      lwd;
#elif defined(USE_WXWIDGETS)
    wxPen    pen;
    wxPen    spen; // solid pen with the same color and width=1
#endif
};

PenType Black("black", "solid", 1);


////////////////////////////////////////////////////////////////////////////////////////////// Pens /////


typedef std::map<String,PenType> Pens_t;

class Pens : public Pens_t
{
public:
    // Constructor
    Pens()
    {
        (*this)["black_solid_1"] = PenType("black","solid",1);
        (*this)["red_solid_1"]   = PenType("red",  "solid",1);
        (*this)["green_solid_1"] = PenType("green","solid",1);
        (*this)["blue_solid_1"]  = PenType("blue", "solid",1);
    }

    // Operators
    PenType operator() (char const * Name)
    {
        Pens_t::const_iterator p = this->find(Name);

        // add new pen
        if (p==this->end())
        {
            // replace "_" with spaces
            String name(Name);
            size_t pos = name.find_first_of("_");
            while (pos!=String::npos)
            {
                name[pos] = ' ';
                pos = name.find_first_of("_",pos+1);
            }

            // parse clr, lty, and lwd
            std::istringstream iss(name);
            String clr, lty;
            int    lwd;
            iss >> clr >> lty >> lwd;

            // add new pen
            PenType pen(clr.CStr(),lty.CStr(),lwd);
            (*this)[Name] = pen;
            return pen;
        }
        else return p->second;
    }
};

// Allocate pens
Pens Pen;


/////////////////////////////////////////////////////////////////////////////////////////// FontType /////


class FontType
{
public:
    // Constructor
    FontType (char const * FntName="arial", int Size=10, bool Bold=false, bool Italic=false) { Set(FntName,Size,Bold,Italic); }

    // Methods
    void Set (char const * FntName, int Size, bool Bold, bool Italic)
    {
#if defined(USE_FLTK)
        if (strcmp(FntName,"arial")==0)
        {
            if (Bold && Italic) type = FL_HELVETICA_BOLD_ITALIC;
            else if (Bold)      type = FL_HELVETICA_BOLD;
            else if (Italic)    type = FL_HELVETICA_ITALIC;
            else                type = FL_HELVETICA;
        }
        else if (strcmp(FntName,"courier")==0)
        {
            if (Bold && Italic) type = FL_COURIER_BOLD_ITALIC;
            else if (Bold)      type = FL_COURIER_BOLD;
            else if (Italic)    type = FL_COURIER_ITALIC;
            else                type = FL_COURIER;
        }
        else throw new Fatal("GUI::FontType::FontType: Font named '%s' is not available",FntName);
        size = Size;
#elif defined(USE_WXWIDGETS)
        if      (strcmp(FntName,"arial")  ==0) font.SetFamily (wxFONTFAMILY_ROMAN);
        else if (strcmp(FntName,"courier")==0) font.SetFamily (wxFONTFAMILY_TELETYPE);
        else throw new Fatal("GUI::FontType::FontType: Font named '%s' is not available",FntName);
        if (Bold)   font.SetWeight (wxFONTWEIGHT_BOLD);
        if (Italic) font.SetStyle  (wxFONTSTYLE_ITALIC);
        font.SetPointSize (Size);
#endif
    }

    void Activate (DeviceContext & DC)
    {
#if defined(USE_FLTK)
        fl_font (type, size);
#elif defined(USE_WXWIDGETS)
        DC.SetFont (font);
#endif
    }

    // Data
#if defined(USE_FLTK)
    int type;
    int size;
#elif defined(USE_WXWIDGETS)
    wxFont font;
#endif
};


///////////////////////////////////////////////////////////////////////////////////////////// Fonts /////


typedef std::map<String,FontType> Fonts_t;

class Fonts : public Fonts_t
{
public:
    // Constructor
    Fonts () : _initialized(false) {}

    // Initialize
    void Initialize ()
    {
        (*this)["arial_10_n"]    = FontType("arial",   10, false, false);
        (*this)["arial_10_b"]    = FontType("arial",   10, true,  false);
        (*this)["arial_10_i"]    = FontType("arial",   10, false, true);
        (*this)["arial_10_bi"]   = FontType("arial",   10, true,  true);
        (*this)["courier_10_n"]  = FontType("courier", 10, false, false);
        (*this)["courier_10_b"]  = FontType("courier", 10, true,  false);
        (*this)["courier_10_i"]  = FontType("courier", 10, false, true);
        (*this)["courier_10_bi"] = FontType("courier", 10, true,  true);
        _initialized = true;
    }

    // Operators
    FontType operator() (char const * Name)
    {
        if (_initialized) Initialize();

        Fonts_t::const_iterator p = this->find(Name);

        // add new font
        if (p==this->end())
        {
            // replace "_" with spaces
            String name(Name);
            size_t pos = name.find_first_of("_");
            while (pos!=String::npos)
            {
                name[pos] = ' ';
                pos = name.find_first_of("_",pos+1);
            }

            // parse fontname, size, bold, italic
            std::istringstream iss(name);
            String fntname, key;
            int    size;
            iss >> fntname >> size >> key;
            bool bold   = false;
            bool italic = false;
            if (key.size()==1)
            {
                if (key[0]=='b') bold   = true;
                if (key[0]=='i') italic = true;
            }
            else if (key.size()==2)
            {
                if (key[0]=='b') bold   = true;
                if (key[1]=='i') italic = true;
            }
            else throw new Fatal("GUI::Fonts(): key for bold/italic == %d is invalid (must be: n, b, i, bi)",key.CStr());

            // add new font
            FontType font(fntname.CStr(),size,bold,italic);
            (*this)[Name] = font;
            return font;
        }
        else return p->second;
    }
private:
    bool _initialized;
};

// Allocate fonts
Fonts Fnt;


////////////////////////////////////////////////////////////////////////////////// Drawing functions /////


#if defined(USE_FLTK)

    // Circle
    #define GUI_DRAW_CIRCLE(dc,x,y,r) fl_circle (x, y, r);

    // Filled circle
    #define GUI_DRAW_FCIRCLE(dc,x,y,r) fl_pie (x-r, y-r, 2*r+1, 2*r+1, 0, 360);

    // Line
    #define GUI_DRAW_LINE(dc,x0,y0,x1,y1) fl_line (x0,y0, x1,y1);

    // Text
    #define GUI_DRAW_TEXT(dc,t,x,y,align) {                    \
        if      (align<0) fl_draw (t,x,y,0,0,FL_ALIGN_LEFT);   \
        else if (align>0) fl_draw (t,x,y,0,0,FL_ALIGN_RIGHT);  \
        else              fl_draw (t,x,y,0,0,FL_ALIGN_CENTER); }

    // Square (centered at x,y)
    #define GUI_DRAW_SQUARE(dc,x,y,l) fl_rect (x-l/2, y-l/2, l, l);

    // Filled square (centered at x,y)
    #define GUI_DRAW_FSQUARE(dc,x,y,l) fl_rectf (x-l/2, y-l/2, l, l);

    // Rectangle
    #define GUI_DRAW_RECT(dc,x,y,w,h) fl_rect (x, y, w, h);

    // Filled rectangle
    #define GUI_DRAW_FRECT(dc,x,y,w,h) fl_rectf (x, y, w, h);

    // Add width of text to p
    #define GUI_DRAW_ADD_WIDTH(dc,t,p) p += fl_width(t);

#elif defined(USE_WXWIDGETS)

    // Circle
    #define GUI_DRAW_CIRCLE(dc,x,y,r) dc.DrawCircle (x, y, r);

    // Filled circle
    #define GUI_DRAW_FCIRCLE(dc,x,y,r) {                \
        dc.SetBrush (wxBrush(dc.GetPen().GetColour())); \
        dc.DrawCircle (x, y, r);                        \
        dc.SetBrush ((*wxTRANSPARENT_BRUSH));           }

    // Line
    #define GUI_DRAW_LINE(dc,x0,y0,x1,y1) dc.DrawLine (x0,y0, x1,y1);

    // Text
    #define GUI_DRAW_TEXT(dc,t,x,y,align) {                \
        wxCoord wid,hei;                                   \
        dc.GetTextExtent (t,&wid,&hei);                    \
        if      (align<0) dc.DrawText (t,x,      y-hei/2); \
        else if (align>0) dc.DrawText (t,x-wid,  y-hei/2); \
        else              dc.DrawText (t,x-wid/2,y-hei/2); }

    // Square (centered at x,y)
    #define GUI_DRAW_SQUARE(dc,x,y,l) dc.DrawRectangle (x-l/2, y-l/2, l, l);

    // Filled square (centered at x,y)
    #define GUI_DRAW_FSQUARE(dc,x,y,l) {                \
        dc.SetBrush (wxBrush(dc.GetPen().GetColour())); \
        dc.DrawRectangle (x-l/2, y-l/2, l, l);          \
        dc.SetBrush ((*wxTRANSPARENT_BRUSH));           }

    // Rectangle
    #define GUI_DRAW_RECT(dc,x,y,w,h) dc.DrawRectangle (x, y, w, h);

    // Filled rectangle
    #define GUI_DRAW_FRECT(dc,x,y,w,h) {                \
        dc.SetBrush (wxBrush(dc.GetPen().GetColour())); \
        dc.DrawRectangle (x, y, w, h);                  \
        dc.SetBrush ((*wxTRANSPARENT_BRUSH));           }

    // Add width of text to p
    #define GUI_DRAW_ADD_WIDTH(dc,t,p) { \
        wxCoord wid,hei;                 \
        dc.GetTextExtent (t,&wid,&hei);  \
        p += wid;                        }

#endif

}; // namespace GUI

#endif // MECHSYS_GUI_COMMON_H
