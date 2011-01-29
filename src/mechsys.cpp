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
#include <boost/python.hpp> // this includes everything

#define USE_BOOST_PYTHON

namespace BPy = boost::python;

// MechSys
#include <mechsys/matfile.h>
#include <mechsys/inpfile.h>
#include <mechsys/linalg/jacobirot.h>
#include <mechsys/numerical/odesolver.h>
#include <mechsys/models/smpinvs.h>

// MechSys -- FEM
#include <mechsys/fem/fem.h>

// functions overloadings
BOOST_PYTHON_FUNCTION_OVERLOADS (FUN_PHI2M,        Phi2M,        1, 2)
BOOST_PYTHON_FUNCTION_OVERLOADS (FUN_M2PHI,        M2Phi,        1, 2)
BOOST_PYTHON_FUNCTION_OVERLOADS (FUN_JACOBIROT,    PyJacobiRot,  3, 4)

// member overloadings
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (OD_Init,         Init,         3, 8)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (SI_Calc,         PyCalc,       1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (IN_SetPrmsInis,  SetPrmsInis,  1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MG_ReadMesh,     ReadMesh,     1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MG_SetVert,      SetVert,      4, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MG_WriteVTU,     WriteVTU,     1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MG_WriteMPY,     WriteMPY,     1, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MS_Generate,     PyGenerate,   1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MS_GenBox,       GenBox,       0, 7)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MU_Generate,     Generate,     0, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MU_GenBox,       GenBox,       0, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MU_WritePLY,     WritePLY,     1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (DO_PrintResults, PrintResults, 0, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (DO_WriteVTU,     WriteVTU,     1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (SO_Solve,        Solve,        0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (SO_DynSolve,     DynSolve,     3, 4)

// module
BOOST_PYTHON_MODULE (mechsys)
{

//////////////////////////////////////////////////////////////////////////////////// util /////

// String
BPy::class_<String>("String")
    .def("PyStr", &String::PyStr)
    .def(BPy::self_ns::str(BPy::self))
    ;

// SDPair
BPy::class_<SDPair>("SDPair")
    .def("Set", &SDPair::PySet)
    .def(BPy::self_ns::str(BPy::self))
    ;

// Dict
BPy::class_<Dict>("Dict")
    .def("Set", &Dict::PySet)
    .def("Get", &Dict::PyGet)
    .def(BPy::self_ns::str(BPy::self))
    ;

// Table
BPy::class_<Table>("Table")
    .def("Set", &Table::PySet)
    .def(BPy::self_ns::str(BPy::self))
    ;

// Fatal
BPy::register_exception_translator<Fatal *>(&PyExceptTranslator);

/////////////////////////////////////////////////////////////////////////////// Numerical /////

BPy::class_<Numerical::PyODESolver>("ODESolver")
    .def("Init",        &Numerical::PyODESolver::Init, OD_Init())
    .def("Evolve",      &Numerical::PyODESolver::Evolve)
    .def_readwrite("Y", &Numerical::PyODESolver::y)
    ;

/////////////////////////////////////////////////////////////////////////////////// linalg ////

BPy::def("pqth2L",    Pypqth2L);
BPy::def("Phi2M",     Phi2M,       FUN_PHI2M());
BPy::def("M2Phi",     M2Phi,       FUN_M2PHI());
BPy::def("JacobiRot", PyJacobiRot, FUN_JACOBIROT());
BPy::def("EigenProjAnalytic", PyEigenProjAnalytic);

/////////////////////////////////////////////////////////////////////////////// Invariants ////

BPy::class_<SMPInvs>("SMPInvs")
    .def("Calc", &SMPInvs::PyCalc, SI_Calc())
    .def_readwrite("b", &SMPInvs::b)
    .def(BPy::self_ns::str(BPy::self))
    ;

///////////////////////////////////////////////////////////////////////////////// MatFile /////

BPy::class_<MatFile>("MatFile")
    .def("Read",      &MatFile::Read)
    .def("GetPrmIni", &MatFile::PyGetPrmIni)
    .def(BPy::self_ns::str(BPy::self))
    ;

///////////////////////////////////////////////////////////////////////////////// InpFile /////

BPy::class_<InpFile>("InpFile")
    .def("Read",        &InpFile::Read)
    .def("SetPrmsInis", &InpFile::SetPrmsInis, IN_SetPrmsInis())
    .def("GetPrmIni",   &InpFile::PyGetPrmIni)
    .def(BPy::self_ns::str(BPy::self))
    .def_readwrite("matid"      , &InpFile::matid     ) //   1
    .def_readwrite("flwid"      , &InpFile::flwid     ) //   2
    .def_readwrite("ninc"       , &InpFile::ninc      ) //   3
    .def_readwrite("cdrift"     , &InpFile::cdrift    ) //   4
    .def_readwrite("stol"       , &InpFile::stol      ) //   5
    .def_readwrite("ssout"      , &InpFile::ssout     ) //   6
    .def_readwrite("ctetg"      , &InpFile::ctetg     ) //   7
    .def_readwrite("fem"        , &InpFile::fem       ) //   8
    .def_readwrite("dyn"        , &InpFile::dyn       ) //   9
    .def_readwrite("hm"         , &InpFile::hm        ) //  10
    .def_readwrite("tf"         , &InpFile::tf        ) //  11
    .def_readwrite("dt"         , &InpFile::dt        ) //  12
    .def_readwrite("dtout"      , &InpFile::dtout     ) //  13
    .def_readwrite("tsw"        , &InpFile::tsw       ) //  14
    .def_readwrite("ndiv"       , &InpFile::ndiv      ) //  15
    .def_readwrite("nip"        , &InpFile::nip       ) //  16
    .def_readwrite("o2"         , &InpFile::o2        ) //  17
    .def_readwrite("ray"        , &InpFile::ray       ) //  18
    .def_readwrite("am"         , &InpFile::am        ) //  19
    .def_readwrite("ak"         , &InpFile::ak        ) //  20
    .def_readwrite("rk"         , &InpFile::rk        ) //  21
    .def_readwrite("rkscheme"   , &InpFile::rkscheme  ) //  22
    .def_readwrite("rkstol"     , &InpFile::rkstol    ) //  23
    .def_readwrite("refdat"     , &InpFile::refdat    ) //  24
    .def_readwrite("refsim"     , &InpFile::refsim    ) //  25
    .def_readwrite("refana"     , &InpFile::refana    ) //  26
    .def_readwrite("idxvert1"   , &InpFile::idxvert1  ) //  27
    .def_readwrite("idxvert2"   , &InpFile::idxvert2  ) //  28
    .def_readwrite("idxvert3"   , &InpFile::idxvert3  ) //  29
    .def_readwrite("optdbl1"    , &InpFile::optdbl1   ) //  30
    .def_readwrite("optdbl2"    , &InpFile::optdbl2   ) //  31
    .def_readwrite("optdbl3"    , &InpFile::optdbl3   ) //  32
    .def_readwrite("hasoptdbl1" , &InpFile::hasoptdbl1) //  30b
    .def_readwrite("hasoptdbl2" , &InpFile::hasoptdbl2) //  31b
    .def_readwrite("hasoptdbl3" , &InpFile::hasoptdbl3) //  32b
    .def_readwrite("nldt_nsml"  , &InpFile::nldt_nsml ) //  33
    .def_readwrite("nldt_nn"    , &InpFile::nldt_nn   ) //  34
    .def_readwrite("nldt_n"     , &InpFile::nldt_n    ) //  35
    .def_readwrite("nldt_ll"    , &InpFile::nldt_ll   ) //  36
    .def_readwrite("nldt_sch"   , &InpFile::nldt_sch  ) //  37
    .def_readwrite("nldt_m"     , &InpFile::nldt_m    ) //  38
    .def_readwrite("maxit"      , &InpFile::maxit     ) //  39
    .def_readwrite("tolr"       , &InpFile::tolr      ) //  40
    .def_readwrite("fnkey"      , &InpFile::fnkey     ) //  41
    .def_readwrite("pcam0"      , &InpFile::pcam0     ) //  42
    .def_readwrite("haspcam0"   , &InpFile::haspcam0  ) //  42b
    .def_readwrite("scheme"     , &InpFile::scheme    ) //  43
    .def_readwrite("vtufile"    , &InpFile::vtufile   ) //  44
    .def_readwrite("suscheme"   , &InpFile::suscheme  ) //  45
    .def_readwrite("sustol"     , &InpFile::sustol    ) //  46
    .def_readwrite("surkscheme" , &InpFile::surkscheme) //  47
    .def_readwrite("dcmaxit"    , &InpFile::dcmaxit   ) //  48
    .def_readwrite("dcftol"     , &InpFile::dcftol    ) //  49
    ;

//////////////////////////////////////////////////////////////////////////////////// mesh /////

// Block
BPy::class_<Mesh::Block>("Block")
    .def("Set", &Mesh::Block::PySet)
    .def(BPy::self_ns::str(BPy::self))
    ;

// Generic
BPy::class_<Mesh::Generic>("Generic","generic mesh", BPy::init<int>())
    .def("ReadMesh",      &Mesh::Generic::ReadMesh, MG_ReadMesh())
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
    .def("PrintResults", &FEM::Domain::PrintResults, DO_PrintResults())
    .def("WriteVTU",     &FEM::Domain::WriteVTU,     DO_WriteVTU())
    .def(BPy::self_ns::str(BPy::self))
    ;

// Solver
BPy::class_<FEM::Solver>("FEM_Solver", "FEM solver", BPy::init<FEM::Domain &>())
    .def("Solve",     &FEM::Solver::Solve,    SO_Solve())
    .def("DynSolve",  &FEM::Solver::DynSolve, SO_DynSolve())
    .def("SetScheme", &FEM::Solver::SetScheme)
    .def_readwrite("MaxIt", &FEM::Solver::MaxIt)
    .def_readwrite("TolR",  &FEM::Solver::TolR)
    ;

} // BOOST_PYTHON_MODULE
