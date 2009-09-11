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

#ifndef MECHSYS_ARRAY_H
#define MECHSYS_ARRAY_H

// STL
#include <algorithm> // for std::find, std::min_element, and std::max_element

// MechSys
#include "util/fatal.h"
#include "util/numstreams.h"

template<typename Value_T>
class Array
{
public:
	// Constructors
	Array ()            : _values(NULL), _ns(&Util::_8s) { Resize(0, 2); } ///< Constructor
	Array (size_t Size) : _values(NULL), _ns(&Util::_8s) { Resize(Size); } ///< Alternative constructor
	Array (Array<Value_T> const & Other);                                  ///< Copy constructor (needed when using Array< Array<...> >)

	/** Destructor. */
	~Array() { if (_values!=NULL) delete [] _values; }

	// Methods
	void                    SetNS     (Util::NumStream & NS)         { _ns = &NS; } ///< Set the NumStream, a structure to aid format output of numbers
	Util::NumStream const & NS        () const                   { return (*_ns); } ///< Return the NumStream, a structure to aid format output of numbers
	size_t                  Size      () const                    { return _size; } ///< Returns the size
	Value_T *               GetPtr    ()                        { return _values; } ///< Returns a pointer to the values
	void                    Resize    (size_t  Size,  double SzFactor=1.2);         ///< Resize the array
	void                    Push      (Value_T const & Value, double SzFactor=1.2); ///< Add a new entry increasing the size if necessary
	void                    PushN     (Value_T const & Value, size_t Num, double SzFactor=1.2); ///< Add a new entry increasing the size if necessary
	void                    Remove    (size_t i, size_t Length=1);                  ///< Remove item i from the array
	long                    Find      (Value_T const & Value) const;                ///< Find a value: returns -1 if not found, otherwise, returns the index of the element found
	long                    Min       () const;                                     ///< Find the minimum value: returns the index of the minimum element
	long                    Max       () const;                                     ///< Find the maximum value: returns the index of the maximum element
	Value_T                 Mean      () const;                                     ///< Calculate the mean value (Value_T must have addition operators)
	Value_T                 Norm      () const;                                     ///< Calculate the norm value (Value_T must have addition operators)
	void                    SetValues (Value_T const & V);                          ///< Set all values to be equal to V
	void                    Clear     () { Resize(0); }                             ///< Clear array

	// Operators
	Value_T       & operator[] (size_t i);                 ///< Access operator (write)
	Value_T const & operator[] (size_t i) const;           ///< Access operator (read)
	void            operator=  (Array<Value_T> const & R); ///< Assignment operator (needed when using Array< Array<...> >)
	void            operator+= (Array<Value_T> const & R); ///< Plus-assignment operator
	void            operator-= (Array<Value_T> const & R); ///< Minus-assignment operator

	// Assign values separated by commas
	class CommaAssign
	{
	public:
		CommaAssign (Array<Value_T> * ptArray, Value_T const & FirstValue) : _ptarray(ptArray), _i(0)
		{
			if (_ptarray->Size()>_i) (*_ptarray)[_i] = FirstValue;
			else throw new Fatal("Array::CommaAssign: The array must be resized with a size greater than or equal to %d before calling this method.",_i+1);
		}
		CommaAssign & operator, (Value_T const & Value)
		{
			_i++;
			if (_ptarray->Size()>_i) (*_ptarray)[_i] = Value;
			else throw new Fatal("Array::CommaAssign: The array must be resized with a size greater than or equal to %d before calling this method.",_i+1);
			return *this;
		}
	private:
		Array<Value_T> * _ptarray;
		size_t           _i;
	};
	CommaAssign operator= (Value_T const & Value) { return CommaAssign(this, Value); }

private:
	// Variables
	size_t            _size;   ///< Current number of components
	size_t            _space;  ///< Available space
	Value_T         * _values; ///< Space to hold all values
	Util::NumStream * _ns;     ///< Structure to aid format output of numbers

}; // class Array


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Constructors

template<typename Value_T>
inline Array<Value_T>::Array(Array<Value_T> const & Other)
	: _values(NULL)
{
	Resize(Other.Size());
	for (size_t i=0; i<_size; ++i)
		_values[i] = Other[i];
}


// Methods

template<typename Value_T>
inline void Array<Value_T>::Resize(size_t Size, double SzFactor)
{
	// Check
	if (SzFactor<1.0) throw new Fatal("Array::Resize: SzFactor==%f must be greater than 1.0", SzFactor);

	// Clear previous memory
	if (_values!=NULL) delete [] _values;

	// Allocate new memory
	_size   = Size;
	_space  = static_cast<size_t>((_size+1)*SzFactor+1);
	_values = new Value_T [_space];
}

template<typename Value_T>
inline void Array<Value_T>::Push(Value_T const & Value, double SzFactor)
{
	if (_size==_space)
	{
		size_t oldsz = _size;
		Value_T * tmp = new Value_T [oldsz];
		for (size_t i=0; i<oldsz; ++i) tmp[i] = _values[i];
		Resize(oldsz+1, SzFactor);
		for (size_t i=0; i<oldsz; ++i) _values[i] = tmp[i];
		delete [] tmp;
	}
	else _size++;
	_values[_size-1] = Value;
}

template<typename Value_T>
inline void Array<Value_T>::PushN(Value_T const & Value, size_t Num, double SzFactor)
{
	if (_size+Num>=_space)
	{
		size_t oldsz = _size;
		Value_T * tmp = new Value_T [oldsz];
		for (size_t i=0; i<oldsz; ++i) tmp[i] = _values[i];
		Resize(oldsz+Num, SzFactor);
		for (size_t i=0; i<oldsz; ++i) _values[i] = tmp[i];
		delete [] tmp;
	}
	else _size += Num;
    for (size_t i=0; i<Num; ++i) _values[_size-Num+i] = Value;
}

template<typename Value_T>
inline void Array<Value_T>::Remove(size_t i, size_t Length)
{
	size_t oldsz = _size;
	Value_T * tmp = new Value_T [oldsz];
	for (size_t j=0; j<oldsz; ++j) tmp[j] = _values[j];
	Resize(_size-Length);
	size_t k = 0;
	for (size_t j=0;        j<i;      ++j) { _values[k]=tmp[j]; k++; }
	for (size_t j=i+Length; j<oldsz;  ++j) { _values[k]=tmp[j]; k++; }
	delete [] tmp;
}

template<typename Value_T>
inline long Array<Value_T>::Find(Value_T const & Value) const
{
	Value_T * res = std::find(_values, _values+_size, Value);
	if (res==_values+_size) return -1;
	else return res-_values;
}

template<typename Value_T>
inline long Array<Value_T>::Min() const
{
	Value_T * res = std::min_element(_values, _values+_size);
	return res-_values;
}

template<typename Value_T>
inline long Array<Value_T>::Max() const
{
	Value_T * res = std::max_element(_values, _values+_size);
	return res-_values;
}

template<typename Value_T>
inline Value_T Array<Value_T>::Mean() const
{
	Value_T sum = 0.0;
	for (size_t i=0; i<_size; ++i) sum += _values[i];
	return sum/_size;
}

template<typename Value_T>
inline Value_T Array<Value_T>::Norm() const
{
	Value_T sum = 0.0;
	for (size_t i=0; i<_size; ++i) sum += _values[i]*_values[i];
	return sqrt(sum);
}

template<typename Value_T>
inline void Array<Value_T>::SetValues(Value_T const & V)
{
	// Set all values to be equal to V
	for (size_t i=0; i<_size; ++i) _values[i] = V;
}


// Operators

template<typename Value_T>
inline Value_T & Array<Value_T>::operator[] (size_t i)
{
#ifndef NDEBUG
	if (i<0 || i>=_size) throw new Fatal("Array::operator[] (write) Subscript==%d (size==%d) is out of range.", i, _size);
#endif
	return _values[i];
}

template<typename Value_T>
inline Value_T const & Array<Value_T>::operator[] (size_t i) const
{
#ifndef NDEBUG
	if (i<0 || i>=_size) throw new Fatal("Array::operator[] (read) Subscript==%d (size==%d) is out of range.", i, _size); 
#endif
	return _values[i];
}

template<typename Value_T>
inline void Array<Value_T>::operator= (Array<Value_T> const & R)
{
#ifndef DNDEBUG
	if (&R==this) throw new Fatal("Array::operator= The right-hand-size of this operation (LHS = RHS) must not be equal to the LHS.");
#endif
	// Reallocate if they are different (LHS != RHS)
	if (_size!=R.Size()) Resize(R.Size());
	
	// Copy values
	for (size_t i=0; i<_size; ++i) _values[i] = R[i];
}

template<typename Value_T>
inline void Array<Value_T>::operator+= (Array<Value_T> const & R)
{
#ifndef DNDEBUG
	if (_size!=R.Size()) throw new Fatal("Array::operator+= The number of components of the LHS (%d) must be equal to the number of components of the RHS (%d).",_size,R.Size());
#endif
	// Add values
	for (int i=0; i<_size; ++i) _values[i] += R[i];
}

template<typename Value_T>
inline void Array<Value_T>::operator-= (Array<Value_T> const & R)
{
#ifndef DNDEBUG
	if (_size!=R.Size()) throw new Fatal("Array::operator-= The number of components of the LHS (%d) must be equal to the number of components of the RHS (%d).",_size,R.Size());
#endif
	// Subtract values
	for (int i=0; i<_size; ++i) _values[i] -= R[i];
}


/** Outputs an array. */
template<typename Value_T>
std::ostream & operator<< (std::ostream & os, const Array<Value_T> & V)
{
	for (size_t i=0; i<V.Size(); ++i)
    {
		os << V[i];
        if (i!=V.Size()-1) os << ", ";
    }
	return os;
}

#endif // MECHSYS_ARRAY_H
