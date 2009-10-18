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

template<typename Value_T>
class Array
{
public:
    // Constructors
    Array ()            : _values(NULL) { Resize(0); }    ///< Constructor
    Array (size_t Size) : _values(NULL) { Resize(Size); } ///< Alternative constructor
    Array (Array<Value_T> const & Other);                 ///< Copy constructor (needed when using Array< Array<...> >)

#ifdef USE_BOOST_PYTHON
    Array (BPy::list const & Dat);
#endif

    // Alternative constructors
    Array (Value_T const & v0, bool JustOne);
    Array (Value_T const & v0, Value_T const & v1);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15, Value_T const & v16);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15, Value_T const & v16, Value_T const & v17);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15, Value_T const & v16, Value_T const & v17, Value_T const & v18);
    Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15, Value_T const & v16, Value_T const & v17, Value_T const & v18, Value_T const & v19);

    /** Destructor. */
    ~Array() { if (_values!=NULL) delete [] _values; }

    // Methods
    size_t    Size      () const { return _size; }                    ///< Returns the size
    Value_T * GetPtr    ()       { return _values; }                  ///< Returns a pointer to the values
    void      Resize    (size_t  Size,  double SzFactor=1.2);         ///< Resize the array
    void      Push      (Value_T const & Value, double SzFactor=1.2); ///< Add a new entry increasing the size if necessary
    void      PushN     (Value_T const & Value, size_t Num,
                         double SzFactor=1.2);                        ///< Add a new entry increasing the size if necessary
    void      Remove    (size_t i, size_t Length=1);                  ///< Remove item i from the array
    long      Find      (Value_T const & Value) const;                ///< Find a value: returns -1 if not found, otherwise, returns the index of the element found
    long      Min       () const;                                     ///< Find the minimum value: returns the index of the minimum element
    long      Max       () const;                                     ///< Find the maximum value: returns the index of the maximum element
    Value_T   Mean      () const;                                     ///< Calculate the mean value (Value_T must have addition operators)
    Value_T   Norm      () const;                                     ///< Calculate the norm value (Value_T must have addition operators)
    void      SetValues (Value_T const & V);                          ///< Set all values to be equal to V
    void      Clear     () { Resize(0); }                             ///< Clear array

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
    size_t    _size;   ///< Current number of components
    size_t    _space;  ///< Available space
    Value_T * _values; ///< Space to hold all values
};


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

#ifdef USE_BOOST_PYTHON

template<typename Value_T>
inline Array<Value_T>::Array (BPy::list const & Dat)
    : _values(NULL)
{
    size_t size = BPy::len(Dat);
    Resize (size);
    for (size_t i=0; i<size; ++i)
        _values[i] = BPy::extract<Value_T>(Dat[i])();
}

#endif

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, bool JustOne)
    : _values(NULL)
{
    Resize(1);
    _values[0] = v0;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1)
    : _values(NULL)
{
    Resize(2);
    _values[0] = v0;
    _values[1] = v1;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2)
    : _values(NULL)
{
    Resize(3);
    _values[0] = v0;
    _values[1] = v1;
    _values[2] = v2;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3)
    : _values(NULL)
{
    Resize(4);
    _values[0] = v0;
    _values[1] = v1;
    _values[2] = v2;
    _values[3] = v3;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4)
    : _values(NULL)
{
    Resize(5);
    _values[0] = v0;
    _values[1] = v1;
    _values[2] = v2;
    _values[3] = v3;
    _values[4] = v4;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5)
    : _values(NULL)
{
    Resize(6);
    _values[0] = v0;
    _values[1] = v1;
    _values[2] = v2;
    _values[3] = v3;
    _values[4] = v4;
    _values[5] = v5;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6)
    : _values(NULL)
{
    Resize(7);
    _values[0] = v0;
    _values[1] = v1;
    _values[2] = v2;
    _values[3] = v3;
    _values[4] = v4;
    _values[5] = v5;
    _values[6] = v6;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7)
    : _values(NULL)
{
    Resize(8);
    _values[0] = v0;
    _values[1] = v1;
    _values[2] = v2;
    _values[3] = v3;
    _values[4] = v4;
    _values[5] = v5;
    _values[6] = v6;
    _values[7] = v7;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8)
    : _values(NULL)
{
    Resize(9);
    _values[0] = v0;
    _values[1] = v1;
    _values[2] = v2;
    _values[3] = v3;
    _values[4] = v4;
    _values[5] = v5;
    _values[6] = v6;
    _values[7] = v7;
    _values[8] = v8;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9)
    : _values(NULL)
{
    Resize(10);
    _values[0] = v0;
    _values[1] = v1;
    _values[2] = v2;
    _values[3] = v3;
    _values[4] = v4;
    _values[5] = v5;
    _values[6] = v6;
    _values[7] = v7;
    _values[8] = v8;
    _values[9] = v9;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10)
    : _values(NULL)
{
    Resize(11);
    _values[ 0] = v0;
    _values[ 1] = v1;
    _values[ 2] = v2;
    _values[ 3] = v3;
    _values[ 4] = v4;
    _values[ 5] = v5;
    _values[ 6] = v6;
    _values[ 7] = v7;
    _values[ 8] = v8;
    _values[ 9] = v9;
    _values[10] = v10;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11)
    : _values(NULL)
{
    Resize(12);
    _values[ 0] = v0;
    _values[ 1] = v1;
    _values[ 2] = v2;
    _values[ 3] = v3;
    _values[ 4] = v4;
    _values[ 5] = v5;
    _values[ 6] = v6;
    _values[ 7] = v7;
    _values[ 8] = v8;
    _values[ 9] = v9;
    _values[10] = v10;
    _values[11] = v11;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12)
    : _values(NULL)
{
    Resize(13);
    _values[ 0] = v0;
    _values[ 1] = v1;
    _values[ 2] = v2;
    _values[ 3] = v3;
    _values[ 4] = v4;
    _values[ 5] = v5;
    _values[ 6] = v6;
    _values[ 7] = v7;
    _values[ 8] = v8;
    _values[ 9] = v9;
    _values[10] = v10;
    _values[11] = v11;
    _values[12] = v12;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13)
    : _values(NULL)
{
    Resize(14);
    _values[ 0] = v0;
    _values[ 1] = v1;
    _values[ 2] = v2;
    _values[ 3] = v3;
    _values[ 4] = v4;
    _values[ 5] = v5;
    _values[ 6] = v6;
    _values[ 7] = v7;
    _values[ 8] = v8;
    _values[ 9] = v9;
    _values[10] = v10;
    _values[11] = v11;
    _values[12] = v12;
    _values[13] = v13;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14)
    : _values(NULL)
{
    Resize(15);
    _values[ 0] = v0;
    _values[ 1] = v1;
    _values[ 2] = v2;
    _values[ 3] = v3;
    _values[ 4] = v4;
    _values[ 5] = v5;
    _values[ 6] = v6;
    _values[ 7] = v7;
    _values[ 8] = v8;
    _values[ 9] = v9;
    _values[10] = v10;
    _values[11] = v11;
    _values[12] = v12;
    _values[13] = v13;
    _values[14] = v14;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15)
    : _values(NULL)
{
    Resize(16);
    _values[ 0] = v0;
    _values[ 1] = v1;
    _values[ 2] = v2;
    _values[ 3] = v3;
    _values[ 4] = v4;
    _values[ 5] = v5;
    _values[ 6] = v6;
    _values[ 7] = v7;
    _values[ 8] = v8;
    _values[ 9] = v9;
    _values[10] = v10;
    _values[11] = v11;
    _values[12] = v12;
    _values[13] = v13;
    _values[14] = v14;
    _values[15] = v15;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15, Value_T const & v16)
    : _values(NULL)
{
    Resize(17);
    _values[ 0] = v0;
    _values[ 1] = v1;
    _values[ 2] = v2;
    _values[ 3] = v3;
    _values[ 4] = v4;
    _values[ 5] = v5;
    _values[ 6] = v6;
    _values[ 7] = v7;
    _values[ 8] = v8;
    _values[ 9] = v9;
    _values[10] = v10;
    _values[11] = v11;
    _values[12] = v12;
    _values[13] = v13;
    _values[14] = v14;
    _values[15] = v15;
    _values[16] = v16;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15, Value_T const & v16, Value_T const & v17)
    : _values(NULL)
{
    Resize(18);
    _values[ 0] = v0;
    _values[ 1] = v1;
    _values[ 2] = v2;
    _values[ 3] = v3;
    _values[ 4] = v4;
    _values[ 5] = v5;
    _values[ 6] = v6;
    _values[ 7] = v7;
    _values[ 8] = v8;
    _values[ 9] = v9;
    _values[10] = v10;
    _values[11] = v11;
    _values[12] = v12;
    _values[13] = v13;
    _values[14] = v14;
    _values[15] = v15;
    _values[16] = v16;
    _values[17] = v17;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15, Value_T const & v16, Value_T const & v17, Value_T const & v18)
    : _values(NULL)
{
    Resize(19);
    _values[ 0] = v0;
    _values[ 1] = v1;
    _values[ 2] = v2;
    _values[ 3] = v3;
    _values[ 4] = v4;
    _values[ 5] = v5;
    _values[ 6] = v6;
    _values[ 7] = v7;
    _values[ 8] = v8;
    _values[ 9] = v9;
    _values[10] = v10;
    _values[11] = v11;
    _values[12] = v12;
    _values[13] = v13;
    _values[14] = v14;
    _values[15] = v15;
    _values[16] = v16;
    _values[17] = v17;
    _values[18] = v18;
}

template<typename Value_T>
inline Array<Value_T>::Array (Value_T const & v0, Value_T const & v1, Value_T const & v2, Value_T const & v3, Value_T const & v4, Value_T const & v5, Value_T const & v6, Value_T const & v7, Value_T const & v8, Value_T const & v9, Value_T const & v10, Value_T const & v11, Value_T const & v12, Value_T const & v13, Value_T const & v14, Value_T const & v15, Value_T const & v16, Value_T const & v17, Value_T const & v18, Value_T const & v19)
    : _values(NULL)
{
    Resize(20);
    _values[ 0] = v0;
    _values[ 1] = v1;
    _values[ 2] = v2;
    _values[ 3] = v3;
    _values[ 4] = v4;
    _values[ 5] = v5;
    _values[ 6] = v6;
    _values[ 7] = v7;
    _values[ 8] = v8;
    _values[ 9] = v9;
    _values[10] = v10;
    _values[11] = v11;
    _values[12] = v12;
    _values[13] = v13;
    _values[14] = v14;
    _values[15] = v15;
    _values[16] = v16;
    _values[17] = v17;
    _values[18] = v18;
    _values[19] = v19;
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
