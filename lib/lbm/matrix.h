/*************************************************************************************
 * MechSys - A C++ library to simulate (Continuum) Mechanical Systems                *
 * Copyright (C) 2005 Dorival de Moraes Pedroso <dorival.pedroso at gmail.com>       *
 * Copyright (C) 2005 Raul Dario Durand Farfan  <raul.durand at gmail.com>           *
 *                                                                                   *
 * This file is part of MechSys.                                                     *
 *                                                                                   *
 * MechSys is free software; you can redistribute it and/or modify it under the      *
 * terms of the GNU General Public License as published by the Free Software         *
 * Foundation; either version 2 of the License, or (at your option) any later        *
 * version.                                                                          *
 *                                                                                   *
 * MechSys is distributed in the hope that it will be useful, but WITHOUT ANY        *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A   *
 * PARTICULAR PURPOSE. See the GNU General Public License for more details.          *
 *                                                                                   *
 * You should have received a copy of the GNU General Public License along with      *
 * MechSys; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, *
 * Fifth Floor, Boston, MA 02110-1301, USA                                           *
 *************************************************************************************/

// -- IMPORTANT --
// OBS.:
//   Internally, the values are stored as Col-Major, thus a pointer
//   can be passed to FORTRAN routines very easily.

#ifndef MECHSYS_LINALG_MATRIX_H
#define MECHSYS_LINALG_MATRIX_H

#include <sstream>
#include <cassert>
#include <cmath>

template<typename Type>
class Matrix
{
public:
	// Default constructor
	Matrix() : _rows(0), _cols(0), _values(NULL) {}
	// Alternative Constructor
	Matrix(int Rows, int Cols)                              : _values(NULL) { _set(Rows,Cols); }
	// Destructor
	~Matrix()
	{
		if (_values!=NULL) delete [] _values;
	}
	// Copy constructor
	Matrix(Matrix<Type> const & Other)
	{
		  _rows = Other.Rows();
		  _cols = Other.Cols();
		_values = new Type [_rows*_cols];
		for (int i=0; i<_rows*_cols; ++i)
			_values[i] = Other._values[i];
	}
	// Get
	       int   Rows()   const { assert(_values!=NULL); return _rows;   }
	       int   Cols()   const { assert(_values!=NULL); return _cols;   }
	      Type * GetPtr()       { assert(_values!=NULL); return _values; }
	const Type * GetPtr() const { assert(_values!=NULL); return _values; }
	// Set
	void Set(int Rows, int Cols)                              { _set(Rows,Cols); }
	void Resize(int Rows, int Cols) { _resize(Rows,Cols); }
	// operators
	void operator= (const Matrix<Type> & R) // {{{
	{
		assert(&R!=this);
		if (R.Rows() != _rows || R.Cols() != _cols)
		{
			_rows = R.Rows();
			_cols = R.Cols();
			if (_values != NULL) delete [] _values;
			_values = new Type [_rows*_cols];
		}
		
		int _components=_rows*_cols;
		for (int i=0; i<_components; ++i)
			_values[i] = R._values[i];
	} // }}}

	Type & operator() (int i) // {{{
	{
		return _values[i]; // Col major
	} // }}}

	Type & operator() (int i, int j) // {{{
	{
		assert(_values!=NULL);
		assert(i>=0 && i<_rows);
		assert(j>=0 && j<_cols);
		return _values[j*_rows+i]; // Col major
	} // }}}

	const Type & operator() (int i, int j) const // {{{
	{
		assert(_values!=NULL);
		assert(i>=0 && i<_rows);
		assert(j>=0 && j<_cols);
		return _values[j*_rows+i]; // Col major
	} // }}}
private:
	int    _rows;
	int    _cols;
	Type * _values;
	void _set(int Rows, int Cols) // {{{
	{
		assert(_values==NULL);
		assert(Rows>0);
		assert(Cols>0);
		  _rows = Rows;
		  _cols = Cols;
		_values = new Type [_rows*_cols];
		//for (int i=0; i<_rows*_cols; ++i)
		//	_values[i] = static_cast<Type>(0);
	} // }}}
	void _resize(int Rows, int Cols) // {{{
	{
		assert(Rows>=0);
		assert(Cols>=0);
	    if (Rows==_rows && Cols==_cols && _values!=NULL) return;	
		if (_values!=NULL) delete [] _values;
		_rows = Rows;
		_cols = Cols;
		if (Rows!=0 && Cols!=0)
			_values = new Type [_rows*_cols];
		else
			_values = NULL;
	} // }}}

}; // class Matrix<Type>

template<typename Type>
std::ostream & operator<< (std::ostream & os, const Matrix<Type> & M) // {{{
{
	os << std::endl;
	for (int i=0; i<M.Rows(); ++i)
	{
		for (int j=0; j<M.Cols(); ++j)
		os << std::endl;
	}
	return os;
} // }}}

#endif // MECHSYS_LINALG_MATRIX_H

// vim:fdm=marker
