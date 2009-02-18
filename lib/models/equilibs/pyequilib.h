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

#ifndef MECHSYS_PYEQUILIB_H
#define MECHSYS_PYEQUILIB_H

// Std Lib
#include <iostream>

// Boost::Python
#include <boost/python.hpp> // this includes everything
namespace BPy = boost::python;

// MechSys
#include "models/equilibmodel.h"
#include "tensors/tensors.h"
#include "util/string.h"
#include "util/util.h"
#include "util/lineparser.h"
#include "util/fileparser.h"

class PyEquilib : public EquilibModel
{
public:
	// Constants
	static const char PYEQUILIB_PN[36][8];

	// Destructor
	virtual ~PyEquilib () {}

	// Derived methods
	int         NPrms   () const { return 36;           }
	PrmName_t * Prms    () const { return PYEQUILIB_PN; }
	Str_t       Name    () const { return "PyEquilib";  }
	void        InitIVS (Ini_t const & Ini, Tensor2 const & Sig, Tensor2 const & Eps, IntVals & Ivs) const;

	/* Set script file name for models that can be extended using Python. */
	void SetPyName (Str_t ScriptFileName) { _py_fn = ScriptFileName; }

private:
	// Data
	String      _py_fn;       ///< Python: script filename
	BPy::object _py_stiff;    ///< Python: implementation of "stiff" function
	BPy::object _py_init_ivs; ///< Python: implementation of "init_ivs" function
	BPy::dict   _py_prms;     ///< Python: parameters structure

	// Private methods
	void _initialize ();
	void _stiff      (Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B, bool First) const;

}; // class PyEquilib

const char PyEquilib::PYEQUILIB_PN[36][8] = { "a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10", "a11", "a12", "a13", "a14", "a15", "a16", "a17", "a18", "a19", "a20", "a21", "a22", "a23", "a24", "a25", "a26", "a27", "a28", "a29", "a30", "a31", "a32", "a33", "a34", "a35" };


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void PyEquilib::InitIVS(Ini_t const & Ini, Tensor2 const & Sig, Tensor2 const & Eps, IntVals & Ivs) const
{
	// Fill python structures
	BPy::tuple sig (BPy::make_tuple( Sig(0), Sig(1), Sig(2), Sig(3), Sig(4), Sig(5)));
	BPy::tuple eps (BPy::make_tuple( Eps(0), Eps(1), Eps(2), Eps(3), Eps(4), Eps(5)));
	BPy::dict  ini;
	for (Ini_t::const_iterator it=Ini.begin(); it!=Ini.end(); ++it) ini[it->first.CStr()] = it->second;

	// Call _py_init_ivs
	BPy::list const & res = BPy::extract<BPy::list>( _py_init_ivs(_py_prms,ini,sig,eps) )();
	
	// Extract return vals
	size_t nivs = BPy::len(res);
	Ivs.Resize (nivs);
	for (size_t i=0; i<nivs; ++i) Ivs[i] = BPy::extract<double>(res[i])();
}

inline void PyEquilib::_initialize()
{
	if (FileParser::CheckForFile(_py_fn)==false) throw new Fatal("PyEquilib::_initialize: Could not find < %s > file (Python script with function < [1;33mD,B = stiff(prms,deps,sig,eps,ivs)[0m > )",_py_fn.CStr());

	// Run script and get _py_stiff
	// TODO: call this only if "stiff" wasn't defined yet
	try
	{
		BPy::object main_module((BPy::handle<>(BPy::borrowed(PyImport_AddModule("__main__")))));
		BPy::object main_namespace = main_module.attr("__dict__");
		BPy::exec_file (_py_fn.CStr(), main_namespace, main_namespace);
		_py_stiff    = BPy::extract<BPy::object>(main_namespace["stiff"   ])();
		_py_init_ivs = BPy::extract<BPy::object>(main_namespace["init_ivs"])();
	}
	catch (BPy::error_already_set const & Err)
	{
		std::cout << "[1;31mFatal: PyEquilib: The following error happend when running script < " << _py_fn.CStr() << " > [0m\n"; 
		PyErr_Print();
		std::cout << "[1;31mFunction < [1;33mz = init_ivs(prms,ini,sig,eps)  and  D,B = stiff(prms,deps,sig,eps,ivs)[0m[1;31m > must be defined.[0m\n";
		throw new Fatal("An error happend when running script < %s >. Please, check console for further information.",_py_fn.CStr());
	}

	// Fill parameters structure
	for (int i=0; i<NPrms(); ++i) _py_prms[PYEQUILIB_PN[i]] = Prm(PYEQUILIB_PN[i]);
}

inline void PyEquilib::_stiff(Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B, bool First) const
{
	// Fill python structures
	BPy::tuple deps (BPy::make_tuple(DEps(0),DEps(1),DEps(2),DEps(3),DEps(4),DEps(5)));
	BPy::tuple  sig (BPy::make_tuple( Sig(0), Sig(1), Sig(2), Sig(3), Sig(4), Sig(5)));
	BPy::tuple  eps (BPy::make_tuple( Eps(0), Eps(1), Eps(2), Eps(3), Eps(4), Eps(5)));
	BPy::list   ivs;
	for (size_t k=0; k<Ivs.Size(); ++k) ivs.append (Ivs[k]);

	// Call _py_stiff
	BPy::tuple const & res = BPy::extract<BPy::tuple>( _py_stiff(_py_prms,deps,sig,eps,ivs) )();

	// Extract return vals
	BPy::tuple const & pyD = BPy::extract<BPy::tuple>(res[0])();
	BPy::list  const & pyB = BPy::extract<BPy::list> (res[1])();

	// Fill C++ structures: pyD => D
	if (BPy::len(pyD)!=36) throw new Fatal("PyEquilib::_stiff: Function < [1;33mD,B = stiff(prms,deps,sig,eps,ivs)[0m > in script < %s > must return a tuple with two values: D and B, in which D is a tuple with 36 values (4th order tensor in Mandel's basis)",_py_fn.CStr());
	size_t m = 0;
	for (size_t i=0; i<6; ++i)
	for (size_t j=0; j<6; ++j)
	{
		D(i,j) = BPy::extract<double>(pyD[m])();
		m++;
	}

	// Fill C++ structures: pyB => B
	for (int k=0; k<BPy::len(pyB); ++k)
	{
		BPy::tuple const & tup = BPy::extract<BPy::tuple>(pyB[k])();
		if (BPy::len(tup)!=6) throw new Fatal("PyEquilib::_stiff: Function < [1;33mD,B = stiff(prms,deps,sig,eps,ivs)[0m > in script < %s > must return a tuple with two values: D and B, in which B is a list of tuples (list of 2nd order tensors -- with 6 values -- in Mandel's basis)",_py_fn.CStr());
		Tensor2 tmp;
		for (size_t i=0; i<6; ++i) tmp(i) = BPy::extract<double>(tup[i])();
		B.Push (tmp);
	}
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new model
Model * PyEquilibMaker() { return new PyEquilib(); }

// Register model
int PyEquilibRegister() { ModelFactory["PyEquilib"]=PyEquilibMaker;  return 0; }

// Call register
int __PyEquilib_dummy_int = PyEquilibRegister();


#endif // MECHSYS_PYEQUILIB_H
