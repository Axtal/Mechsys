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
#include "fem/data.h"
#include "fem/element.h"
#include "fem/solver.h"
#include "util/exception.h"

// MechSys -- Solvers
#include "fem/solvers/autome.h"
#include "fem/solvers/forwardeuler.h"

// MechSys -- Elements
#include "fem/elems/elasticrod.h"
#include "fem/elems/tri6equilib.h"
#include "fem/elems/hex8equilib.h"

// MechSys -- Models
#include "models/equilibs/linelastic.h"

using FEM::Nodes;
using FEM::Elems;
using std::cout;
using std::endl;

using namespace boost::python;

// Wrapper functions
void addnode_2d (double X, double Y)           { FEM::AddNode (X,Y);   }
void addnode_3d (double X, double Y, double Z) { FEM::AddNode (X,Y,Z); }
void addelem    (str Type, bool IsActive)      { FEM::AddElem (extract<char const *>(Type),IsActive); }

// Exceptions
void except_translator (Exception * e)
{
	String msg;
	msg.Printf("[1;31m%s[0m",e->Msg().GetSTL().c_str());
	PyErr_SetString(PyExc_UserWarning, msg.GetSTL().c_str());
	//if (e->IsFatal()) {delete e; exit(1);}
	delete e;
}

// Extra Python functions
void printelems ()
{
	cout << "[1;34mMechSys:[0m Elements available: " << endl;
	FEM::ElementFactory_t::const_iterator it;
	for (it=FEM::ElementFactory.begin(); it!=FEM::ElementFactory.end(); it++)
	{
		cout << "\t" << it->first << endl;
	}
}

BOOST_PYTHON_MODULE (mechsys)
{
	// Global functions
    def ("addnode",    addnode_2d);
    def ("addnode",    addnode_3d);
    def ("addelem",    addelem   );
    def ("printelems", printelems);

	// Exceptions
	register_exception_translator<Exception *>(&except_translator);
}
