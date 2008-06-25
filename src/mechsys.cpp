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

#define USE_BOOST_PYTHON

// STL
#include <iostream>

// Boost-Python
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python.hpp> // this includes everything

// MechSys -- mesh
#include "mesh/structured.h"

// MechSys -- fem -- basic
#include "fem/data.h"
#include "fem/element.h"
#include "fem/solver.h"
#include "models/model.h"
#include "util/exception.h"

// MechSys -- fem -- Solvers
#include "fem/solvers/autome.h"
#include "fem/solvers/forwardeuler.h"

// MechSys -- fem -- Elements
#include "fem/elems/rod.h"
#include "fem/elems/tri6pstrain.h"
#include "fem/elems/tri6pstress.h"
#include "fem/elems/quad4pstrain.h"
#include "fem/elems/quad4pstress.h"
#include "fem/elems/quad8pstrain.h"
#include "fem/elems/quad8pstress.h"
#include "fem/elems/hex8equilib.h"

// MechSys -- Models
#include "models/equilibs/linelastic.h"

// Python string to char const *
#define S2C(py_str) extract<char const *>(py_str)

using std::cout;
using std::endl;

using namespace boost::python;

////////////////////////////////////////////////////////////////////////////////////////// Wrapper classes

class PyNode
{
public:
	PyNode (int iNode) : _id(iNode) {}
	PyNode & bry (str Name, double Value) { FEM::Nodes[_id]->Bry(S2C(Name), Value);  return (*this); }
	double   val (str Name)               { return FEM::Nodes[_id]->Val(S2C(Name)); }
private:
	int _id;
}; // class PyNode

class PyElement
{
public:
	PyElement (int iElem) : _id(iElem) {}
	PyElement & set_node  (int iNodeLocal, int iNodeGlobal) { FEM::Elems[_id]->SetNode(iNodeLocal, iNodeGlobal);        return (*this); }
	PyElement & set_model (str Name, str Prms, str Inis)    { FEM::Elems[_id]->SetModel(S2C(Name),S2C(Prms),S2C(Inis)); return (*this); }
	double      val       (int iNodeLocal, str Name)        { return FEM::Elems[_id]->Val(iNodeLocal,S2C(Name)); }
private:
	int _id;
}; // class PyElement

class PySolver
{
public:
	PySolver (str Name) { sol = FEM::AllocSolver(S2C(Name)); }
	void       solve          ()                 { sol->Solve(); }
	PySolver & set_lin_sol    (str Key)          { sol->SetLinSol   (S2C(Key));  return (*this); }
	PySolver & set_num_div    (int Numdiv)       { sol->SetNumDiv   (Numdiv);    return (*this); }
	PySolver & set_delta_time (double DeltaTime) { sol->SetDeltaTime(DeltaTime); return (*this); }
private:
	FEM::Solver * sol;
}; // class PySolver

//////////////////////////////////////////////////////////////////////////////////////// Wrapper functions

void      add_node_2d   (double X, double Y)           { FEM::AddNode (X,Y); }
void      add_node_3d   (double X, double Y, double Z) { FEM::AddNode (X,Y,Z); }
void      add_elem      (str Type, bool IsActive)      { FEM::AddElem (S2C(Type),IsActive); }
void      dim           (int  nDim)                    { FEM::Dim = nDim; }
PyNode    nodes         (int  iNode)                   { PyNode    tmp(iNode); return tmp; }
PyElement elems         (int  iElem)                   { PyElement tmp(iElem); return tmp; }

/////////////////////////////////////////////////////////////////////////////////// Extra Python functions

void except_translator (Exception * e)
{
	String msg;
	msg.Printf("[1;31m%s[0m",e->Msg().GetSTL().c_str());
	PyErr_SetString(PyExc_UserWarning, msg.GetSTL().c_str());
	//if (e->IsFatal()) {delete e; exit(1);}
	delete e;
}

void print_elems ()
{
	cout << "[1;34mMechSys:[0m Elements available: " << endl;
	FEM::ElementFactory_t::const_iterator it;
	for (it=FEM::ElementFactory.begin(); it!=FEM::ElementFactory.end(); it++)
	{
		cout << "\t" << it->first << endl;
	}
}

void print_models ()
{
	cout << "[1;34mMechSys:[0m Models available: " << endl;
	ModelFactory_t::const_iterator it;
	for (it=ModelFactory.begin(); it!=ModelFactory.end(); it++)
	{
		cout << "\t" << it->first << endl;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////// Define Module

BOOST_PYTHON_MODULE (mechsys)
{
	// Global classes

	// --------------------------------------------------------------------------- Mesh

	class_<PyMeshBlock>("mesh_block")
	    .def("set",       PMBSet1)
	    .def("set",       PMBSet2)
	    .def("set_etags", &PyMeshBlock::SetETags)
	    .def("set_ftags", &PyMeshBlock::SetFTags)
	    ;

	class_<PyMeshStruct>("mesh_struct")
	    .def(init<double>())
	    .def("generate",  &PyMeshStruct::Generate)
	    .def("write_vtu", &PyMeshStruct::WriteVTU)
	    .def("get_verts", &PyMeshStruct::GetVerts)
	    .def("get_elems", &PyMeshStruct::GetElems)
	    .def("get_etags", &PyMeshStruct::GetETags)
	    ;

	// ---------------------------------------------------------------------------- FEM
	
	class_<PyNode>("node", init<int>())
	    .def("bry", &PyNode::bry, return_internal_reference<>())
	    .def("val", &PyNode::val)
	    ;

	class_<PyElement>("element", init<int>())
	    .def("set_node",  &PyElement::set_node , return_internal_reference<>())
	    .def("set_model", &PyElement::set_model, return_internal_reference<>())
	    .def("val",       &PyElement::val)
	    ;

	class_<PySolver>("solver", init<str>())
	    .def("solve",          &PySolver::solve)
	    .def("set_lin_sol",    &PySolver::set_lin_sol,    return_internal_reference<>())
	    .def("set_num_div",    &PySolver::set_num_div,    return_internal_reference<>())
	    .def("set_delta_time", &PySolver::set_delta_time, return_internal_reference<>())
	    ;

	// Global functions
	def ("add_node",      add_node_2d  );
	def ("add_node",      add_node_3d  );
	def ("add_elem",      add_elem     );
	def ("dim",           dim          );
	def ("nodes",         nodes        );
	def ("elems",         elems        );
	def ("print_elems",   print_elems  );
	def ("print_models",  print_models );

	// Exceptions
	register_exception_translator<Exception *>(&except_translator);
}
