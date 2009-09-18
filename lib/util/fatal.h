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

#ifndef MECHSYS_FATAL_H
#define MECHSYS_FATAL_H

// Std lib
#include <iostream> // for cout
#include <cstdarg>  // for va_list, va_start, va_end
#include <exception>

// Boost::Python
#ifdef USE_BOOST_PYTHON
  #include <boost/python.hpp> // this includes everything
  namespace BPy = boost::python;
#endif

// MechSys
#include "util/string.h"


#ifdef USE_BOOST_PYTHON
  #define MECHSYS_CATCH catch (Fatal      * e)         { e->Cout();  delete e;  exit(1); }                                             \
                        catch (char const * m)         { std::cout<<"[1;31mFatal: "<<m<<"[0m\n";  exit(1); }                       \
                        catch (BPy::error_already_set) { std::cout<<"[1;31mFatal: "; PyErr_Print(); std::cout<<"[0m\n"; exit(1); } \
                        catch (std::exception & e)     { std::cout<<"[1;31mFatal: "<<e.what()                <<"[0m\n"; exit(1); } \
                        catch (...)                    { std::cout<<"[1;31mFatal: Some exception (...) ocurred[0m\n"; exit(1); }
#else
  #define MECHSYS_CATCH catch (Fatal      * e)     { e->Cout();  delete e;  exit(1); }                                           \
                        catch (char const * m)     { std::cout<<"[1;31mFatal: "<<m<<"[0m\n";  exit(1); }                     \
                        catch (std::exception & e) { std::cout<<"[1;31mFatal: "<<e.what()              <<"[0m\n"; exit(1); } \
                        catch (...)                { std::cout<<"[1;31mFatal: Some exception (...) ocurred[0m\n"; exit(1); }
#endif
 

class Fatal
{
public:
	// Constructor
	Fatal (String const & Fmt, ...);

	// Methods
	void   Cout () const { std::cout << "[1;31m" << "Fatal: " << _msg.CStr() << "[0m" << std::endl; }
	String Msg  () const { return _msg; }

private:
	String _msg;
};



/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Fatal::Fatal(String const & Fmt, ...)
{
	va_list       arg_list;
	va_start     (arg_list, Fmt);
	_msg.PrintfV (Fmt, arg_list);
	va_end       (arg_list);
}


#ifdef USE_BOOST_PYTHON

void PyExceptTranslator (Fatal * e)
{
	PyErr_SetString(PyExc_UserWarning, e->Msg().CStr());
	delete e;
}

#endif // USE_BOOST_PYTHON


#endif // MECHSYS_FATAL_H
