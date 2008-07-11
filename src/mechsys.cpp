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
#include <boost/python/def.hpp>
#include <boost/python.hpp> // this includes everything

// MechSys -- mesh
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
	
	class_<PyNode>("node", init<FEM::Node *>())
	    .def("bry",     &PyNode::Bry, return_internal_reference<>())
	    .def("val",     &PyNode::Val)
	    .def("x",       &PyNode::X)
	    .def("y",       &PyNode::Y)
	    .def("z",       &PyNode::Z)
	    .def("get_id",  &PyNode::GetID)
	    ;

	class_<PyElem>("elem", init<FEM::Element *>())
	    .def("connect",   &PyElem::Connect,  return_internal_reference<>())
	    .def("set_model", &PyElem::SetModel, return_internal_reference<>())
	    .def("bry",       &PyElem::Bry,      return_internal_reference<>())
	    .def("val",       PElemVal1)
	    .def("val",       PElemVal2)
	    .def("nnodes",    &PyElem::nNodes)
	    .def("nod",       &PyElem::Nod)
	    .def("get_id",    &PyElem::GetID)
	    ;

	class_<PyGeom>("geom", init<int>())
	    .def("set_nnodes", &PyGeom::SetNNodes)
	    .def("set_nelems", &PyGeom::SetNElems)
	    .def("set_node",   PGSetNode1)
	    .def("set_node",   PGSetNode2)
	    .def("set_elem",   PGSetElem1)
	    .def("set_elem",   PGSetElem2)
	    .def("nnodes",     &PyGeom::nNodes)
	    .def("nelems",     &PyGeom::nElems)
	    .def("nod",        &PyGeom::Nod)
	    .def("ele",        &PyGeom::Ele)
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
	def ("add_nodes_elems",   PyAddNodesElems  );
	def ("set_node_brys",     PySetNodeBrys    );
	def ("set_face_brys",     PySetFaceBrys    );
	def ("write_vtu_equilib", PyWriteVTUEquilib);
	def ("write_vtk",         PyWriteVTK       );
	def ("set_geom",          PySetGeom        );

	// ---------------------------------------------------------------------- Exceptions
	
	register_exception_translator<Exception *>(&PyExceptTranslator);
}
