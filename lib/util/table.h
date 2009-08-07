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

#ifndef MECHSYS_TABLE_H
#define MECHSYS_TABLE_H

// Std lib
#include <iostream> // for cout
#include <sstream> // for cout
#include <cstdarg>  // for va_list, va_start, va_end

// MechSys
#include "util/string.h"
#include "util/numstreams.h"

typedef std::map< String,Array<double> > Table_t;

class Table : public Table_t
{
public:
	// Constructor
	Table (const char * StrKeys, size_t NumRows, ...);

	// Data
	Array<String> Keys;
	size_t        NRows;

}; // class Table


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Table::Table(const char * StrKeys, size_t NumRows, ...)
	: NRows (NumRows)
{
	// retrieve keys and initialize table
	String             key;
	std::istringstream iss(StrKeys);
	while (iss>>key)
	{
		Keys.Push(key);
		(*this)[key].Resize(NumRows);
	}

	// read values
	va_list   arg_list;
	va_start (arg_list, NumRows);
	for (size_t i=0; i<NumRows; ++i)
	{
		for (size_t j=0; j<Keys.Size(); ++j)
		{
			(*this)[Keys[j]][i] = va_arg(arg_list,double);
		}
	}
	//(*this)[key][i] = 1.0;//va_arg(arg_list,double);
	va_end (arg_list);
}

std::ostream & operator<< (std::ostream & os, Table const & T)
{
	// keys
	for (size_t i=0; i<T.Keys.Size(); ++i) os << Util::_8s << T.Keys[i];
	os << "\n";

	// values
	for (size_t i=0; i<T.NRows; ++i)
	{
		for (size_t j=0; j<T.Keys.Size(); ++j)
		{
			Table_t::const_iterator p = T.find(T.Keys[j]);
			os << Util::_8s << p->second[i];
		}
		os << "\n";
	}

	return os;
}


#endif // MECHSYS_TABLE_H
