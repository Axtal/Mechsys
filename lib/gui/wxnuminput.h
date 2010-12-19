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

#ifndef MECHSYS_WXNUMINPUT_H
#define MECHSYS_WXNUMINPUT_H

// wxWidgets
#include <wx/textctrl.h>
#include <wx/event.h>
#include <wx/validate.h>


///////////////////////////////////////////////////////////////////////////////////// WxNumInput_Event /////


class WxNumInput_Event : public wxNotifyEvent
{
public:
    // Constructor
    WxNumInput_Event(wxEventType CmdType=wxEVT_NULL, int Id=0)
        : wxNotifyEvent(CmdType, Id)
    {}

    // Copy constructor
    WxNumInput_Event(WxNumInput_Event const & Event)
        : wxNotifyEvent(Event)
    {}

    // Clone function
    virtual wxEvent * Clone() const { return new WxNumInput_Event(*this); }

    // Declare dynamic class
    DECLARE_DYNAMIC_CLASS(WxNumInput_Event);

    // Methods
    void   SetValue (double Val) { _value = Val;  }
    double GetValue ()     const { return _value; }

private:
    // Data
    double _value;
};

typedef void (wxEvtHandler::*WxNumInput_Event_Fun)(WxNumInput_Event&);

BEGIN_DECLARE_EVENT_TYPES()
    DECLARE_EVENT_TYPE(EVT_REALNUM_CHANGED, 801)
END_DECLARE_EVENT_TYPES()

// Macros for handling events
#define EVT_REALNUM_CHANGED(id, fn)                                           \
    DECLARE_EVENT_TABLE_ENTRY                                                 \
    (                                                                         \
        EVT_REALNUM_CHANGED                                                 , \
        id                                                                  , \
        -1                                                                  , \
        (wxObjectEventFunction)(wxEventFunction)(WxNumInput_Event_Fun) & fn , \
        (wxObject*) NULL                                                      \
    ),

DEFINE_EVENT_TYPE       (EVT_REALNUM_CHANGED)
IMPLEMENT_DYNAMIC_CLASS (WxNumInput_Event, wxNotifyEvent)


/////////////////////////////////////////////////////////////////////////////////////////// WxNumInput /////


class WxNumInput : public wxTextCtrl
{
public:
    // Constructor
    WxNumInput(wxWindow * Parent, wxWindowID Id, double Val, const wxPoint & Pos=wxDefaultPosition, const wxSize & Size=wxDefaultSize, int Style=0, wxValidator const & Validator=wxDefaultValidator)
        : wxTextCtrl(Parent, Id, wxEmptyString, Pos, Size, Style|wxTE_PROCESS_ENTER, Validator)
    {
        wxString str;  str.Printf("%g",Val);
        SetValue (str);
    }

    // Events
    void OnChar      (wxKeyEvent     & Event);
    void OnSetFocus  (wxFocusEvent   & Event);
    void OnKillFocus (wxFocusEvent   & Event);
    void OnTextEnter (wxCommandEvent & Event); /// TODO: correct this: It seems to send changed event twice!

    // Methods
    double GetVal () const { double val; GetValue().ToDouble(&val); return val; }
    void   SetVal (double Val) { wxString str; str.Printf("%g",Val); SetValue(str); }

private:
    // Event table
    DECLARE_EVENT_TABLE()

    // Data
    wxString _buffer;

    // Methods
    bool _validate_value             (double & TmpValue);
    void _send_realnum_changed_event (double Val);
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void WxNumInput::OnChar(wxKeyEvent & Event)
{
    if (!isalpha(Event.GetKeyCode())) Event.Skip(); // Input is ok, so skip to main loop
    else
    {
        if (Event.GetKeyCode()==69)  Event.Skip(); // "E"
        if (Event.GetKeyCode()==101) Event.Skip(); // "e"
    }
}

inline void WxNumInput::OnSetFocus(wxFocusEvent & Event)
{
    _buffer = GetValue(); // save a buffered value
    Event.Skip();
}

inline void WxNumInput::OnKillFocus(wxFocusEvent & Event)
{
    double temp;
    if (_validate_value(temp)) _send_realnum_changed_event(temp);
}

inline void WxNumInput::OnTextEnter(wxCommandEvent & Event)
{
    double temp;
    if (_validate_value(temp)) _send_realnum_changed_event(temp);
}

inline bool WxNumInput::_validate_value(double & OutValue)
{
    bool ok_and_changed = false;
    if (!GetValue().ToDouble(&OutValue)) SetValue(_buffer); // not valid, so restore buffered value
    else
    {
        if (GetValue()!=_buffer) // really changed ?
        {
            double buffer_temp;
            if (_buffer.ToDouble(&buffer_temp))
            {
                if (buffer_temp==OutValue) SetValue(_buffer); // did not change, just the format is different, ex.: 0 != 000
                else ok_and_changed=true;
            }
        }
    }
    return ok_and_changed;
}

inline void WxNumInput::_send_realnum_changed_event(double Val)
{
    WxNumInput_Event event(EVT_REALNUM_CHANGED, GetId());
    event.SetEventObject (this);
    event.SetValue       (Val);
    GetEventHandler()->ProcessEvent(event);
}

BEGIN_EVENT_TABLE(WxNumInput, wxTextCtrl)
    EVT_CHAR       (          WxNumInput::OnChar     )
    EVT_SET_FOCUS  (          WxNumInput::OnSetFocus )
    EVT_KILL_FOCUS (          WxNumInput::OnKillFocus)
    EVT_TEXT_ENTER (wxID_ANY, WxNumInput::OnTextEnter)
END_EVENT_TABLE()


/////////////////////////////////////////////////////////////////////////////////////////// Validator /////


class WxNumValidator : public wxValidator
{
public:
    WxNumValidator (double * ptVal=NULL) : wxValidator(), pVal(ptVal) {}
    WxNumValidator (WxNumValidator const & Other) { pVal = Other.pVal; }
    virtual ~WxNumValidator() {}

    virtual wxObject * Clone              () const { return new WxNumValidator(*this); }
    virtual bool       TransferFromWindow ();
    virtual bool       TransferToWindow   ();
    virtual bool       Validate           (wxWindow * parent) { return true; }

    WxNumValidator & operator= (WxNumValidator const & Other) { pVal = Other.pVal; return (*this); }

private:
    double * pVal;
    DECLARE_DYNAMIC_CLASS(WxNumValidator)
    DECLARE_EVENT_TABLE()
};

IMPLEMENT_DYNAMIC_CLASS(WxNumValidator, wxValidator)

BEGIN_EVENT_TABLE(WxNumValidator, wxValidator)
END_EVENT_TABLE()

bool WxNumValidator::TransferFromWindow()
{
    (*pVal) = static_cast<WxNumInput*>(m_validatorWindow)->GetVal();
    return true;
}

bool WxNumValidator::TransferToWindow()
{
    static_cast<WxNumInput*>(m_validatorWindow)->SetVal((*pVal));
    return true;
}

#endif // MECHSYS_WXNUMINPUT_H
