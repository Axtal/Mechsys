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

// Boost-Python
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/module_init.hpp>
#include <boost/python/def.hpp>
#include <boost/python/call_method.hpp>
#include <boost/ref.hpp>
#include <boost/utility.hpp>
#include <boost/python/list.hpp>
#include <boost/python/dict.hpp>

#define USE_BOOST_PYTHON

namespace BPy = boost::python;

// MechSys
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/mesh/mesh.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/mesh/unstructured.h>
#include <mechsys/fem/element.h>
#include <mechsys/fem/rod.h>
#include <mechsys/fem/beam.h>
#include <mechsys/fem/equilibelem.h>
#include <mechsys/fem/geomelem.h>
#include <mechsys/fem/elems/tri3.h>
#include <mechsys/fem/elems/tri6.h>
#include <mechsys/fem/elems/quad4.h>
#include <mechsys/fem/elems/quad8.h>
#include <mechsys/fem/elems/hex8.h>
#include <mechsys/fem/elems/hex20.h>
#include <mechsys/fem/elems/tet10.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solver.h>
#include <mechsys/models/model.h>
#include <mechsys/models/linelastic.h>
#include <mechsys/models/elastoplastic.h>
#include <mechsys/models/camclay.h>
#include <mechsys/linalg/matvec.h>

// MechSys -- DEM
#include <mechsys/dem/domain.h>

// overloadings
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MG_SetVert,      SetVert,      4, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MG_WriteVTU,     WriteVTU,     1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MG_WriteMPY,     WriteMPY,     1, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MS_Generate,     PyGenerate,   1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MS_GenBox,       GenBox,       0, 7)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MS_GenQRing,     GenQRing,     0, 9)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MU_Generate,     Generate,     0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MU_GenBox,       GenBox,       0, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MU_WritePLY,     WritePLY,     1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (DO_PrintResults, PrintResults, 0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (DO_WriteMPY,     WriteMPY,     1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (SO_Solve,        Solve,        0, 1)

// module
BOOST_PYTHON_MODULE (mechsys)
{

//////////////////////////////////////////////////////////////////////////////////// util /////

// SDPair
BPy::class_<SDPair>("SDPair")
    .def("Set", &SDPair::PySet)
    .def(BPy::self_ns::str(BPy::self))
    ;

// Dict
BPy::class_<Dict>("Dict")
    .def("Set", &Dict::PySet)
    .def(BPy::self_ns::str(BPy::self))
    ;

// Table
BPy::class_<Table>("Table")
    .def("Set", &Table::PySet)
    .def(BPy::self_ns::str(BPy::self))
    ;

// Fatal
BPy::register_exception_translator<Fatal *>(&PyExceptTranslator);

/////////////////////////////////////////////////////////////////////////////////// linalg ////

BPy::def("pqt2L", Pypqt2L);

//////////////////////////////////////////////////////////////////////////////////// mesh /////

// Block
BPy::class_<Mesh::Block>("Block")
    .def("Set", &Mesh::Block::PySet)
    .def(BPy::self_ns::str(BPy::self))
    ;

// Generic
BPy::class_<Mesh::Generic>("Generic","generic mesh", BPy::init<int>())
    .def("SetSize",       &Mesh::Generic::SetSize)
    .def("SetVert",       &Mesh::Generic::SetVert,  MG_SetVert())
    .def("SetCell",       &Mesh::Generic::PySetCell)
    .def("SetBryTag",     &Mesh::Generic::SetBryTag)
    .def("AddLinCells",   &Mesh::Generic::PyAddLinCells)
    .def("WriteVTU",      &Mesh::Generic::WriteVTU, MG_WriteVTU())
    .def("WriteMPY",      &Mesh::Generic::WriteMPY, MG_WriteMPY())
    .def("GetVertsEdges", &Mesh::Generic::PyGetVertsEdges)
    .def(BPy::self_ns::str(BPy::self))
    ;

// Structured
BPy::class_<Mesh::Structured, BPy::bases<Mesh::Generic> >("Structured","structured mesh", BPy::init<int>())
    .def("Generate", &Mesh::Structured::PyGenerate, MS_Generate())
    .def("GenBox",   &Mesh::Structured::GenBox,     MS_GenBox())
    .def("GenQRing", &Mesh::Structured::GenQRing,   MS_GenQRing())
    .def(BPy::self_ns::str(BPy::self))
    ;

// Unstructured
BPy::class_<Mesh::Unstructured, BPy::bases<Mesh::Generic> >("Unstructured","Unstructured mesh", BPy::init<int>())
    .def("Set",      &Mesh::Unstructured::PySet)
    .def("Generate", &Mesh::Unstructured::Generate, MU_Generate())
    .def("GenBox",   &Mesh::Unstructured::GenBox,   MU_GenBox())
    .def("WritePLY", &Mesh::Unstructured::WritePLY, MU_WritePLY())
    .def(BPy::self_ns::str(BPy::self))
    ;

///////////////////////////////////////////////////////////////////////////////////// fem /////

// PROB, GEOM, and MODEL
BPy::def("PROB",  PyPROB);
BPy::def("GEOM",  PyGEOM);
BPy::def("MODEL", PyMODEL);

// Domain
BPy::class_<FEM::Domain>("FEM_Domain", "FEM domain", BPy::init<Mesh::Generic const &, Dict const &, Dict const &, Dict const &>())
    .def("SetBCs",       &FEM::Domain::SetBCs)
    .def("SetOutNods",   &FEM::Domain::PySetOutNods)
    .def("SetOutEles",   &FEM::Domain::PySetOutEles)
    .def("PrintResults", &FEM::Domain::PrintResults, DO_PrintResults())
    .def("WriteMPY",     &FEM::Domain::WriteMPY,     DO_WriteMPY())
    .def("WriteVTU",     &FEM::Domain::WriteVTU)
    .def("CheckError",   &FEM::Domain::CheckError)
    .def(BPy::self_ns::str(BPy::self))
    ;

// Solver
BPy::class_<FEM::Solver>("FEM_Solver", "FEM solver", BPy::init<FEM::Domain const &>())
    .def("Solve", &FEM::Solver::Solve, SO_Solve())
    ;

///////////////////////////////////////////////////////////////////////////////////// dem /////

BPy::class_<DEM::Domain>("DEM_Domain")
    .def("GenSpheres",  &DEM::Domain::GenSpheres)
    .def("GenBox",      &DEM::Domain::GenBox)
    .def("AddVoroPack", &DEM::Domain::AddVoroPack)
    .def("AddSphere",   &DEM::Domain::PyAddSphere)
    .def("AddCube",     &DEM::Domain::PyAddCube)
    .def("AddTetra",    &DEM::Domain::PyAddTetra)
    .def("AddRice",     &DEM::Domain::PyAddRice)
    .def("AddPlane",    &DEM::Domain::PyAddPlane)
    .def("Solve",       &DEM::Domain::Solve)
    .def("WritePOV",    &DEM::Domain::WritePOV)
    .def("WriteBPY",    &DEM::Domain::WriteBPY)
    .def("SetCamPos",   &DEM::Domain::PySetCamPos)
    ;

} // BOOST_PYTHON_MODULE
