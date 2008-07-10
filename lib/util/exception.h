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

#ifndef MECHSYS_EXCEPTION_H
#define MECHSYS_EXCEPTION_H

#include <iostream> // for cout
#include <cstdarg>  // for va_list, va_start, va_end

#include "util/string.h"

class Exception
{
public:
	// Destructor
	virtual ~Exception() {}
	// Methods
	virtual void Cout     () const =0;
	virtual bool IsFatal  () const =0;
	virtual bool IsWarning() const =0;
	String       Msg      () const { return _msg; }
protected:
	// Data
	String _msg;
};

class Message : public Exception
{
public:
	// Constructor
	Message(String const & Fmt, ...);
	// Methods
	void Cout     () const { std::cout << "[1;32m" << "Message: " << _msg.GetSTL() << "[0m" << std::endl; }
	bool IsFatal  () const { return false; }
	bool IsWarning() const { return false; }
}; // class Message

class Warning : public Exception
{
public:
	// Constructor
	Warning(String const & Fmt, ...);
	// Methods
	void Cout     () const { std::cout << "[34m" << "Warning: " << _msg.GetSTL() << "[0m" << std::endl; }
	bool IsFatal  () const { return false; }
	bool IsWarning() const { return true;  }
}; // class Warning

class Fatal : public Exception
{
public:
	// Constructor
	Fatal(String const & Fmt, ...);
	// Methods
	void Cout     () const { std::cout << "[1;31m" << "Fatal: " << _msg.GetSTL() << "[0m" << std::endl; }
	bool IsFatal  () const { return true;  }
	bool IsWarning() const { return false; }
}; // class Fatal


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Message::Message(String const & Fmt, ...)
{
	va_list       arg_list;
	va_start     (arg_list, Fmt);
	_msg.PrintfV (Fmt, arg_list);
	va_end       (arg_list);
}

inline Warning::Warning(String const & Fmt, ...)
{
	va_list       arg_list;
	va_start     (arg_list, Fmt);
	_msg.PrintfV (Fmt, arg_list);
	va_end       (arg_list);
}

inline Fatal::Fatal(String const & Fmt, ...)
{
	va_list       arg_list;
	va_start     (arg_list, Fmt);
	_msg.PrintfV (Fmt, arg_list);
	va_end       (arg_list);
}


#ifdef USE_BOOST_PYTHON
// {

void PyExceptTranslator (Exception * e)
{
	String msg;
	msg.Printf("[1;31m%s[0m",e->Msg().GetSTL().c_str());
	PyErr_SetString(PyExc_UserWarning, msg.GetSTL().c_str());
	//if (e->IsFatal()) {delete e; exit(1);}
	delete e;
}

// }
#endif // USE_BOOST_PYTHON


#endif // MECHSYS_EXCEPTION_H
