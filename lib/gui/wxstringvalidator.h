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

#ifndef MECHSYS_WXSTRINGVALIDATOR_H
#define MECHSYS_WXSTRINGVALIDATOR_H

// wxWidgets
#include <wx/textctrl.h>
#include <wx/validate.h>

// MechSys
#include <mechsys/util/string.h>

class WxStringValidator : public wxValidator
{
public:
    WxStringValidator (String * ptStr=NULL) : wxValidator(), pStr(ptStr) {}
    WxStringValidator (WxStringValidator const & Other) { pStr=Other.pStr; }
    virtual ~WxStringValidator() {}

    virtual wxObject * Clone              () const { return new WxStringValidator(*this); }
    virtual bool       TransferFromWindow ();
    virtual bool       TransferToWindow   ();
    virtual bool       Validate           (wxWindow * parent) { return true; }

    WxStringValidator & operator= (WxStringValidator const & Other) { pStr=Other.pStr; return (*this); }

private:
    String * pStr;
    DECLARE_DYNAMIC_CLASS(WxStringValidator)
    DECLARE_EVENT_TABLE()
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation


IMPLEMENT_DYNAMIC_CLASS(WxStringValidator, wxValidator)

BEGIN_EVENT_TABLE(WxStringValidator, wxValidator)
END_EVENT_TABLE()

bool WxStringValidator::TransferFromWindow()
{
    (*pStr) = static_cast<wxTextCtrl*>(m_validatorWindow)->GetValue().ToStdString();
    return true;
}

bool WxStringValidator::TransferToWindow()
{
    static_cast<wxTextCtrl*>(m_validatorWindow)->SetValue((*pStr));
    return true;
}

#endif
