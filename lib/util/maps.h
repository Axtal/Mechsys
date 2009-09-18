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

#ifndef MECHSYS_MAPS_H
#define MECHSYS_MAPS_H

// Std lib
#include <cstring>  // for strcmp
#include <iostream> // for cout
#include <sstream>  // for istringstream, ostringstream
#include <cstdarg>  // for va_list, va_start, va_end
#include <cmath>    // for fabs
#include <map>

// MechSys
#include "util/array.h"
#include "util/string.h"
#include "util/numstreams.h"

// Constants
#define TRUE  1.0
#define FALSE 1.0


/////////////////////////////////////////////////////////////////////////////////////////// SDPair


typedef std::map<String,double> StrDbl_t;

/** String-Double Pair. */
class SDPair : public StrDbl_t
{
public:
    /** Ex:  pair("ux uy", 1.0,2.0)        OR   Case-A
     *       pair("ux=1.0  uy=2.0")        OR   Case-B
     *       pair("{'ux':1.0, 'uy':2.0}")       Case-C  (Python-Dict) */
    SDPair & Set (const char * Str, ...);
    SDPair & Set (const char * Str, va_list ArgList);

    // Operators
    double operator() (char   const * Key) const;
    double operator() (String const & Key) const;

    // Methods
    bool HasKey  (char const   * Key) const;
    bool HasKey  (String const & Key) const;
    void Val2Key (double Val, String & Key, double Tol=1.0e-15) const;
    void clear   () { Keys.Clear(); StrDbl_t::clear(); }

    // Data
    Array<String> Keys;
};


/////////////////////////////////////////////////////////////////////////////////////////// SIPair


typedef std::map<String,int> StrInt_t;

/** String-Int Pair. */
class SIPair : public StrInt_t
{
public:
    /** Ex:  pair("a b", 1,2)        OR   Case-A
     *       pair("a=1  b=2")        OR   Case-B
     *       pair("{'a':1, 'b':2}")       Case-C  (Python-Dict) */
    SIPair & Set (const char * Str, ...);

    // Operators
    double operator() (char const   * Key) const;
    double operator() (String const & Key) const;

    // Methods
    bool HasKey (char const   * Key) const;
    bool HasKey (String const & Key) const;
    void clear  () { Keys.Clear(); StrInt_t::clear(); }

    // Data
    Array<String> Keys;
};


/////////////////////////////////////////////////////////////////////////////////////////// ISPair


typedef std::map<int,char const *> IntStr_t;

/** Int-String Pair. */
class ISPair : public IntStr_t
{
public:
    /** Ex:  pair("-1 -2", "alpha","beta")    OR   Case-A
     *       pair("{-1:'alpha', -2:'beta'}")       Case-B  (Python-Dict) */
    ISPair & Set (const char * Str, ...);

    // Operators
    char const * operator() (int Key) const;

    // Methods
    bool HasKey (int Key) const;
    void clear  () { Keys.Clear(); IntStr_t::clear(); }

    // Data
    Array<int> Keys;
};


/////////////////////////////////////////////////////////////////////////////////////////// IDPair


typedef std::map<int,double> IntDbl_t;

/** Int-Double Pair. */
class IDPair : public IntDbl_t
{
public:
    /** Ex:  pair("-1 -2", 0.0,0.1)    OR   Case-A
     *       pair("{-1:0.0, -2:0.1}")       Case-B  (Python-Dict) */
    IDPair & Set (const char * Str, ...);

    // Operators
    double operator() (int Key) const;

    // Methods
    bool HasKey (int Key) const;
    void clear  () { Keys.Clear(); IntDbl_t::clear(); }

    // Data
    Array<int> Keys;
};


/////////////////////////////////////////////////////////////////////////////////////////// Dict


typedef std::map<int, SDPair> Dict_t;

/** Dictionary. */
class Dict : public Dict_t
{
public:
    /** Ex:  dict(-1, "ux uy", 1.0,2.0)        OR   Case-A
     *       dict(-1, "ux=1.0  uy=2.0")        OR   Case-B
     *       dict(-1, "{'ux':1.0, 'uy':2.0}")       Case-C  (Python-Dict) */
    Dict & Set (int Key, const char * Str, ...);

    // Operators
    SDPair const & operator() (int Key) const;

    // Methods
    bool HasKey (int Key) const;
    void clear  () { Keys.Clear(); Dict_t::clear(); }

    // Data
    Array<int> Keys;
};


/////////////////////////////////////////////////////////////////////////////////////////// Table


typedef std::map<String,Array<double> > Table_t;

class Table : public Table_t
{
public:
    /** Ex: Set("ux uy", 2, 1.0, 2.0,
     *                      3.0, 4.0);  */
    Table & Set     (const char * StrKeys, size_t NumRows, ...);
    Table & SetZero (const char * StrKeys, size_t NumRows);

    // Operators
    double operator() (String const & Key, size_t iRow) const;
    void   SetVal     (String const & Key, size_t iRow, double Value);

    // Data
    size_t        NRows;
    Array<String> Keys;
};


///////////////////////////////////////////////////////////////////////////////////// operator << //////////////


std::ostream & operator<< (std::ostream & os, SDPair const & P)
{
    int nkeys = P.size();
    int k     = 0;
    os << "{";
    for (size_t i=0; i<P.Keys.Size(); ++i)
    {
        String key = P.Keys[i];
        String val;  val.Printf("%g",P(key));
        os << "'" << key << "':" << val;
        if (k<nkeys-1) os << ", ";
        k++;
    }
    os << "}";
    return os;
}

std::ostream & operator<< (std::ostream & os, SIPair const & P)
{
    int nkeys = P.size();
    int k     = 0;
    os << "{";
    for (size_t i=0; i<P.Keys.Size(); ++i)
    {
        String key = P.Keys[i];
        String val;  val.Printf("%d",P(key));
        os << "'" << key << "':" << val;
        if (k<nkeys-1) os << ", ";
        k++;
    }
    os << "}";
    return os;
}

std::ostream & operator<< (std::ostream & os, ISPair const & P)
{
    int nkeys = P.size();
    int k     = 0;
    os << "{";
    for (size_t i=0; i<P.Keys.Size(); ++i)
    {
        int    key = P.Keys[i];
        String val = P(key);
        os << "" << key << ":" << "'" << val << "'";
        if (k<nkeys-1) os << ", ";
        k++;
    }
    os << "}";
    return os;
}

std::ostream & operator<< (std::ostream & os, IDPair const & P)
{
    int nkeys = P.size();
    int k     = 0;
    os << "{";
    for (size_t i=0; i<P.Keys.Size(); ++i)
    {
        int    key = P.Keys[i];
        double val = P(key);
        os << "" << key << ":" << val;
        if (k<nkeys-1) os << ", ";
        k++;
    }
    os << "}";
    return os;
}

std::ostream & operator<< (std::ostream & os, Dict const & D)
{
    int nkeys = D.size();
    int k = 0;
    os << "{";
    for (size_t i=0; i<D.Keys.Size(); ++i)
    {
        int key = D.Keys[i];
        os << key << ":" << D(key);
        if (k<nkeys-1) os << ",\n ";
        k++;
    }
    os << "}";
    return os;
}

std::ostream & operator<< (std::ostream & os, Table const & T)
{
    // keys
    for (size_t i=0; i<T.Keys.Size(); ++i) os << Util::_8s << T.Keys[i];
    os << "\n";

    // values
    for (size_t i=0; i<T.NRows; ++i)
    {
        for (size_t j=0; j<T.Keys.Size(); ++j) os << Util::_8s << T(T.Keys[j],i);
        os << "\n";
    }

    return os;
}


/////////////////////////////////////////////////////////////////////////////////// SDPair: Implementation /////


inline SDPair & SDPair::Set(const char * Str, ...)
{
    String str(Str);
    if (str.find(":")!=String::npos) // Case-C
    {
        // replace ":',{}" with spaces
        size_t pos = str.find_first_of(":',{}");
        while (pos!=String::npos)
        {
            str[pos] = ' ';
            pos      = str.find_first_of(":',{}",pos+1);
        }

        // fill map
        std::istringstream iss(str);
        String key;
        double val;
        while (iss>>key)
        {
            if (Keys.Find(key)<0) Keys.Push(key);
            if (iss>>val) (*this)[key] = val;
            else break;
        }
    }
    else if (str.find("=")!=String::npos) // Case-B
    {
        // replace "=" with spaces
        size_t pos = str.find_first_of("=");
        while (pos!=String::npos)
        {
            str[pos] = ' ';
            pos      = str.find_first_of("=",pos+1);
        }

        // fill map
        std::istringstream iss(str);
        String key;
        double val;
        while (iss>>key)
        {
            if (Keys.Find(key)<0) Keys.Push(key);
            if (iss>>val) (*this)[key] = val;
            else break;
        }
    }
    else // Case-A
    {
        std::istringstream iss(Str);
        String    key;
        va_list   arg_list;
        va_start (arg_list, Str);
        while (iss>>key)
        {
            if (Keys.Find(key)<0) Keys.Push(key);
            (*this)[key] = va_arg(arg_list,double);
        }
        va_end (arg_list);
    }
    return (*this);
}

inline SDPair & SDPair::Set(const char * Str, va_list ArgList)
{
    std::istringstream iss(Str);
    String key;
    while (iss>>key)
    {
        if (Keys.Find(key)<0) Keys.Push(key);
        (*this)[key] = va_arg(ArgList,double);
    }
    return (*this);
}

inline double SDPair::operator() (char const * Key) const
{
    StrDbl_t::const_iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("SDPair::operator(): String-Double pair: %s does not have a key = '%s'",oss.str().c_str(),Key);
    }
    return p->second;
}

inline double SDPair::operator() (String const & Key) const
{
    StrDbl_t::const_iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("SDPair::operator(): String-Double pair: %s does not have a key = '%s'",oss.str().c_str(),Key.CStr());
    }
    return p->second;
}

inline bool SDPair::HasKey (char const * Key) const
{
    StrDbl_t::const_iterator p = this->find(Key);
    return (p!=this->end());
}

inline bool SDPair::HasKey (String const & Key) const
{
    StrDbl_t::const_iterator p = this->find(Key);
    return (p!=this->end());
}

inline void SDPair::Val2Key (double Val, String & Key, double Tol) const
{
    bool found = false;
    for (StrDbl_t::const_iterator p=this->begin(); p!=this->end(); ++p)
    {
        if (fabs(Val-p->second)<Tol)
        {
            Key   = p->first;
            found = true;
            break;
        }
    }
    if (!found)
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("SDPair::Val2Key: Could not find Val=%g in map: %s",Val,oss.str().c_str());
    }
}


/////////////////////////////////////////////////////////////////////////////////// SIPair: Implementation /////


inline SIPair & SIPair::Set(const char * Str, ...)
{
    String str(Str);
    if (str.find(":")!=String::npos) // Case-C
    {
        // replace ":',{}" with spaces
        size_t pos = str.find_first_of(":',{}");
        while (pos!=String::npos)
        {
            str[pos] = ' ';
            pos      = str.find_first_of(":',{}",pos+1);
        }

        // fill map
        std::istringstream iss(str);
        String key;
        int    val;
        while (iss>>key)
        {
            if (Keys.Find(key)<0) Keys.Push(key);
            if (iss>>val) (*this)[key] = val;
            else break;
        }
    }
    else if (str.find("=")!=String::npos) // Case-B
    {
        // replace "=" with spaces
        size_t pos = str.find_first_of("=");
        while (pos!=String::npos)
        {
            str[pos] = ' ';
            pos      = str.find_first_of("=",pos+1);
        }

        // fill map
        std::istringstream iss(str);
        String key;
        int    val;
        while (iss>>key)
        {
            if (Keys.Find(key)<0) Keys.Push(key);
            if (iss>>val) (*this)[key] = val;
            else break;
        }
    }
    else // Case-A
    {
        std::istringstream iss(Str);
        String    key;
        va_list   arg_list;
        va_start (arg_list, Str);
        while (iss>>key)
        {
            if (Keys.Find(key)<0) Keys.Push(key);
            (*this)[key] = va_arg(arg_list,int);
        }
        va_end (arg_list);
    }
    return (*this);
}

inline double SIPair::operator() (char const * Key) const
{
    StrInt_t::const_iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("SIPair::operator(): String-Int pair: %s does not have a key = '%s'",oss.str().c_str(),Key);
    }
    return p->second;
}

inline double SIPair::operator() (String const & Key) const
{
    StrInt_t::const_iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("SIPair::operator(): String-Int pair: %s does not have a key = '%s'",oss.str().c_str(),Key.CStr());
    }
    return p->second;
}

inline bool SIPair::HasKey (char const * Key) const
{
    StrInt_t::const_iterator p = this->find(Key);
    return (p!=this->end());
}

inline bool SIPair::HasKey (String const & Key) const
{
    StrInt_t::const_iterator p = this->find(Key);
    return (p!=this->end());
}


/////////////////////////////////////////////////////////////////////////////////// ISPair: Implementation /////


inline ISPair & ISPair::Set(const char * Str, ...)
{
    String str(Str);
    if (str.find(":")!=String::npos) // Case-B
    {
        // replace ":',{}" with spaces
        size_t pos = str.find_first_of(":',{}");
        while (pos!=String::npos)
        {
            str[pos] = ' ';
            pos      = str.find_first_of(":',{}",pos+1);
        }

        // fill map
        std::istringstream iss(str);
        int    key;
        String val;
        while (iss>>key)
        {
            if (Keys.Find(key)<0) Keys.Push(key);
            if (iss>>val) (*this)[key] = val.CStr();
            else break;
        }
    }
    else // Case-A
    {
        std::istringstream iss(Str);
        int       key;
        va_list   arg_list;
        va_start (arg_list, Str);
        while (iss>>key)
        {
            if (Keys.Find(key)<0) Keys.Push(key);
            (*this)[key] = va_arg(arg_list,char const *);
        }
        va_end (arg_list);
    }
    return (*this);
}

inline char const * ISPair::operator() (int Key) const
{
    IntStr_t::const_iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("ISPair::operator(): Int-String pair: %s does not have a key = %d",oss.str().c_str(),Key);
    }
    return p->second;
}

inline bool ISPair::HasKey (int Key) const
{
    if (Keys.Find(Key)<0) return false;
    else                  return true;
}


/////////////////////////////////////////////////////////////////////////////////// IDPair: Implementation /////


inline IDPair & IDPair::Set(const char * Str, ...)
{
    String str(Str);
    if (str.find(":")!=String::npos) // Case-B
    {
        // replace ":',{}" with spaces
        size_t pos = str.find_first_of(":',{}");
        while (pos!=String::npos)
        {
            str[pos] = ' ';
            pos      = str.find_first_of(":',{}",pos+1);
        }

        // fill map
        std::istringstream iss(str);
        int    key;
        double val;
        while (iss>>key)
        {
            if (Keys.Find(key)<0) Keys.Push(key);
            if (iss>>val) (*this)[key] = val;
            else break;
        }
    }
    else // Case-A
    {
        std::istringstream iss(Str);
        int       key;
        va_list   arg_list;
        va_start (arg_list, Str);
        while (iss>>key)
        {
            if (Keys.Find(key)<0) Keys.Push(key);
            (*this)[key] = va_arg(arg_list,double);
        }
        va_end (arg_list);
    }
    return (*this);
}

inline double IDPair::operator() (int Key) const
{
    IntDbl_t::const_iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("IDPair::operator(): Int-Double pair: %s does not have a key = %d",oss.str().c_str(),Key);
    }
    return p->second;
}

inline bool IDPair::HasKey (int Key) const
{
    if (Keys.Find(Key)<0) return false;
    else                  return true;
}


///////////////////////////////////////////////////////////////////////////////////// Dict: Implementation /////


inline Dict & Dict::Set (int Key, const char * Str, ...)
{
    String str(Str);
         if (str.find(":")!=String::npos) { (*this)[Key].Set(Str); } // Case-C
    else if (str.find("=")!=String::npos) { (*this)[Key].Set(Str); } // Case-B
    else // Case-A
    {
        va_list   arg_list;
        va_start (arg_list, Str);
        (*this)[Key].Set(Str,arg_list);
        va_end (arg_list);
    }
    if (Keys.Find(Key)<0)
    {
        Keys.Push(Key);
    }
    return (*this);
}

inline SDPair const & Dict::operator() (int Key) const
{
    Dict_t::const_iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("Dict::operator(): Dictionary: %s does not have a key = %d",oss.str().c_str(),Key);
    }
    return p->second;
}

inline bool Dict::HasKey (int Key) const
{
    if (Keys.Find(Key)<0) return false;
    else                  return true;
}


//////////////////////////////////////////////////////////////////////////////////// Table: Implementation /////


inline Table & Table::Set (const char * StrKeys, size_t NumRows, ...)
{
    NRows = NumRows;

    // retrieve keys and initialize table
    String             key;
    std::istringstream iss(StrKeys);
    while (iss>>key)
    {
        if (Keys.Find(key)<0) Keys.Push(key);
        (*this)[key].Resize (NRows);
    }

    // read values
    va_list   arg_list;
    va_start (arg_list, NumRows);
    for (size_t i=0; i<NumRows; ++i)
    {
        for (size_t j=0; j<Keys.Size(); ++j)
        {
            (*this)[Keys[j]][i] = va_arg(arg_list,double);
        }
    }
    va_end (arg_list);
    return (*this);
}

inline Table & Table::SetZero (const char * StrKeys, size_t NumRows)
{
    NRows = NumRows;

    // retrieve keys and initialize table
    String             key;
    std::istringstream iss(StrKeys);
    while (iss>>key)
    {
        if (Keys.Find(key)<0) Keys.Push(key);
        (*this)[key].Resize    (NRows);
        (*this)[key].SetValues (0.0);
    }
    return (*this);
}

inline double Table::operator() (String const & Key, size_t iRow) const
{
    Table_t::const_iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("Table::operator(): Table: %s does not have a key = '%s'",oss.str().c_str(),Key.CStr());
    }
    return p->second[iRow];
}

inline void Table::SetVal (String const & Key, size_t iRow, double Value)
{
    Table_t::iterator p = this->find(Key);
    if (p==this->end())
    {
        std::ostringstream oss;
        oss << (*this);
        throw new Fatal("Table::operator(): Table: %s does not have a key = '%s'",oss.str().c_str(),Key.CStr());
    }
    p->second[iRow] = Value;
}


#endif // MECHSYS_MAPS_H