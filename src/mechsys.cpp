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
#include <cfloat>    // for DBL_EPSILON

// Boost-Python
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/module_init.hpp>
#include <boost/python/def.hpp>
#include <boost/python/call_method.hpp>
#include <boost/ref.hpp>
#include <boost/utility.hpp>
#include <boost/python.hpp> // this includes everything

// MechSys -- mesh
#include "mesh/mesh.h"
#include "mesh/structured.h"

// MechSys -- fem -- basic
#include "fem/geometry.h"
#include "fem/functions.h"
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

using namespace boost::python;

BOOST_PYTHON_MODULE (mechsys)
{
	// --------------------------------------------------------------------------- Mesh

	class_<Mesh::Generic>("mesh_generic")
	    .def(init<double>())
	    .def("write_vtu",   &Mesh::Generic::PyWriteVTU)
	    .def("get_verts",   &Mesh::Generic::PyGetVerts)
	    .def("get_elems",   &Mesh::Generic::PyGetElems)
	    .def("get_eletags", &Mesh::Generic::PyGetEleTags)
	    .def("get_etags",   &Mesh::Generic::PyGetETags)
	    .def("set_nverts",  &Mesh::Generic::PySetNVerts)
	    .def("set_nelems",  &Mesh::Generic::PySetNElems)
	    .def("set_vert",    &Mesh::Generic::PySetVert2D)
	    .def("set_vert",    &Mesh::Generic::PySetVert3D)
	    .def("set_elem",    &Mesh::Generic::PySetElem)
	    .def(self_ns::str(self))
	    ;

	class_<Mesh::Block>("mesh_block")
	    .def("set2d",     &Mesh::Block::PySet2D)
	    .def("set3d",     &Mesh::Block::PySet3D)
	    .def("set_etags", &Mesh::Block::PySetETags)
	    .def("set_ftags", &Mesh::Block::PySetFTags)
	    ;

	class_<Mesh::Structured>("mesh_structured")
	    .def(init<double>())
	    .def("write_vtu",   &Mesh::Structured::PyWriteVTU)
	    .def("get_verts",   &Mesh::Structured::PyGetVerts)
	    .def("get_elems",   &Mesh::Structured::PyGetElems)
	    .def("get_eletags", &Mesh::Structured::PyGetEleTags)
	    .def("get_etags",   &Mesh::Structured::PyGetETags)
	    .def("generate",    &Mesh::Structured::PyGenerate)
	    .def(self_ns::str(self))
	    ;

	// ---------------------------------------------------------------------------- FEM
	
	class_<FEM::Node>("node")
	    .def("get_id",  &FEM::Node::GetID)
	    .def("x",       &FEM::Node::X)
	    .def("y",       &FEM::Node::Y)
	    .def("z",       &FEM::Node::Z)
	    .def("bry",     &FEM::Node::PyBry, return_internal_reference<>())
	    .def("val",     &FEM::Node::PyVal)
	    .def(self_ns::str(self))
	    ;

	class_<PyElem>("elem", init<FEM::Element *>())
	    .def("get_id",    &PyElem::GetID)
	    .def("nnodes",    &PyElem::nNodes)
	    .def("nod",       &PyElem::Nod,      return_internal_reference<>())
	    .def("connect",   &PyElem::Connect,  return_internal_reference<>())
	    .def("set_model", &PyElem::SetModel, return_internal_reference<>())
	    .def("bry",       &PyElem::Bry,      return_internal_reference<>())
	    .def("val",       &PyElem::Val1)
	    .def("val",       &PyElem::Val2)
	    .def(self_ns::str(self))
	    ;

	class_<FEM::Geom>("geom", init<int>())
	    .def("set_nnodes", &FEM::Geom::SetNNodes)
	    .def("set_nelems", &FEM::Geom::SetNElems)
	    .def("nnodes",     &FEM::Geom::nNodes)
	    .def("nelems",     &FEM::Geom::nElems)
	    .def("set_node",   &FEM::Geom::PySetNode2D, return_internal_reference<>())
	    .def("set_node",   &FEM::Geom::PySetNode3D, return_internal_reference<>())
	    .def("set_elem",   &FEM::Geom::PySetElem1)
	    .def("set_elem",   &FEM::Geom::PySetElem2)
	    .def("nod",        &FEM::Geom::PyNod,       return_internal_reference<>())
	    .def("ele",        &FEM::Geom::PyEle)
	    .def(self_ns::str(self))
	    ;

	class_<PySolver>("solver", init<str const &>())
	    .def("set_geom",       &PySolver::SetGeom,      return_internal_reference<>())
	    .def("set_lin_sol",    &PySolver::SetLinSol,    return_internal_reference<>())
	    .def("set_num_div",    &PySolver::SetNumDiv,    return_internal_reference<>())
	    .def("set_delta_time", &PySolver::SetDeltaTime, return_internal_reference<>())
	    .def("solve",          &PySolver::Solve)
	    ;

	// ----------------------------------------------------------------------- functions
	
	// Global functions
	def ("write_vtu_equilib", PyWriteVTUEquilib  );
	def ("write_vtk",         PyWriteVTK         );
	def ("set_geom",          PySetGeom          );
	def ("set_geom",          PySetGeomStructured);

	// ---------------------------------------------------------------------- Exceptions
	
	register_exception_translator<Exception *>(&PyExceptTranslator);
}
