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

// STL
#include <iostream>

// Boost-Python
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python.hpp>

// MechSys
#include "mesh/structured.h"
#include "linalg/matrix.h"
#include "util/exception.h"

// Python string to char const *
#define S2C(py_str) extract<char const *>(py_str)

using std::cout;
using std::endl;

using namespace boost::python;

inline int Max2(int A, int B)        { return (A>B ? A       : B      ); }
inline int Max3(int A, int B, int C) { return (A>B ? A>C?A:C : B>C?B:C); }

////////////////////////////////////////////////////////////////////////////////////////// Wrapper classes


//////////////////////////////////////////////////////////////////////////////////////// Wrapper functions


/////////////////////////////////////////////////////////////////////////////////// Extra Python functions

void except_translator (Exception * e)
{
	String msg;
	msg.Printf("[1;31m%s[0m",e->Msg().GetSTL().c_str());
	PyErr_SetString(PyExc_UserWarning, msg.GetSTL().c_str());
	//if (e->IsFatal()) {delete e; exit(1);}
	delete e;
}

//////////////////////////////////////////////////////////////////////////////////////////// Define Module

BOOST_PYTHON_MODULE (msysmesh)
{
	// Global classes


	// Global functions

	// Exceptions
	register_exception_translator<Exception *>(&except_translator);
}
