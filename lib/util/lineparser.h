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

#ifndef MECHSYS_LINEPARSER_H
#define MECHSYS_LINEPARSER_H

#include <cstdlib> // for getenv
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "util/array.h"
#include "util/string.h"
#include "util/exception.h"

class LineParser : public std::istringstream
{
public:
	LineParser(std::string const & Line) : std::istringstream(Line)        {}
	LineParser(     String const & Line) : std::istringstream(Line.CStr()) {}
	LineParser(       char const * Line) : std::istringstream(Line)        {}
	void Reset(     String const & Line) { Reset(Line.CStr()); }
	void Reset(       char const * Line)
	{
		(*this).clear();
		(*this).str(Line);
	}
	void Reset(std::string const & Line)
	{
		(*this).clear();
		(*this).str(Line);
	}
	void ReplaceAllChars(char OldChar, char NewChar)
	{
		std::string tmp(std::istringstream::str());
		std::replace(tmp.begin(), tmp.end(), OldChar, NewChar);
		std::istringstream::str(tmp);
	}

	void Set(String const & S)
	{
		this->str(S.CStr());
	}

	template<typename Type>
	void ToArray(Array<Type> & A);

	template<typename Type>
	void StructuredLine(int & A, int & B, Array< Array<Type> > & C);

	void SplitLine(std::string const & Separator, Array<std::string> & Res);
	void SplitLine(     String const & Separator, Array<     String> & Res);

	template<typename Type1, typename Type2>
	bool BreakExpressions(Array<Type1> & lvalue, Array<Type2> & rvalue);

	void PathSubstituteEnv();
	void FileBasename(String const & ExtensionToRemove, String & Basename);

private:
}; // class LineParser


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


template<typename Type>
inline void LineParser::ToArray(Array<Type> & A)
{
	// Fill A array with iss.str() values
	A.Resize(0);
	Type tmp;
	while ((*this)>>tmp)
	{
		A.Push(tmp);
	}
}

template<typename Type>
inline void LineParser::StructuredLine(int & A, int & B, Array< Array<Type> > & C)
{
	// Ex.:
	//      7  5  { 1 2 3 } { 4 5 6 8 0 1 }
	//      A  B     C[0]        C[1]

	C.Resize(0);
	std::string str_open_key;
	(*this) >> A >> B >> str_open_key;
	if (str_open_key[0]!='{')
		throw new Fatal(_("LineParser::StructuredLine: Line is not corrected formatted. Line=< %s >"),this->str().c_str());
	
	Array<Type> a_inner;
	std::string str_tmp;
	while ((*this)>>str_tmp)
	{
		if (str_tmp[0]=='}')
		{
			(*this)>>str_open_key; // read next '{'
			C.Push(a_inner);
			a_inner.Resize(0);
		}
		else
		{
			std::istringstream iss(str_tmp);
			Type tmp; iss>>tmp;
			a_inner.Push(tmp);
		}
	}
}

inline void LineParser::SplitLine(std::string const & Separator, Array<std::string> & Res)
{
	/* Ex.:
	 *         /home/dorival/teste/An File.txt 
	 *          R[0]  R[1]   R[2]     R[3]
	 *
	 *         ->1234->456a->789as
	 *           R[0]  R[1]  R[2]
	 */

	// Clear result array
	Res.Resize(0);

	// Fill result array
	if (Separator.empty())
		Res.Push(this->str());
	else if (!this->str().empty())
	{
		// Loop along the line
		size_t start = 0;
		size_t end   = this->str().find(Separator, start);
		while (end!=std::string::npos)
		{
			Res.Push(this->str().substr(start, end - start));
			start = end + Separator.size();
			  end = this->str().find(Separator, start);
		}
		Res.Push(this->str().substr(start));
		// Clear the line (all extracted)
		this->clear(); this->str("");
	}
}

inline void LineParser::SplitLine(String const & Separator, Array<String> & Res)
{
	std::string        sep = Separator.CStr();
	Array<std::string> arr;
	SplitLine(sep,arr);
	size_t len = arr.Size();
	Res.Resize(len);
	for (size_t i=0; i<len; ++i)
		Res[i] = String(arr[i]);
}

template<typename Type1, typename Type2>
inline bool LineParser::BreakExpressions(Array<Type1> & Lvalue, Array<Type2> & Rvalue)
{
	/* For:  this->str() == "var1=value1 var2=value3 var3=value4"
	 *
	 * with:  var1,   var2   and var3   of Type1
	 * and    value1, value2 and value3 of Type2
	 *
	 * result:
	 *           Lvalue={var1,var2,var3}
	 *      and  Rvalue={value1,value2,value3}
	 * return:
	 *           false if some problem occurred with istringstream, i.e: format is invalid
	 */

	Lvalue.Resize(0);
	Rvalue.Resize(0);
	ReplaceAllChars('=',' ');

	// Parse bry marks
	Type1 lval;
	Type2 rval;
	while ((*this)>>lval)
	{
		if (((*this)>>rval))
		{
			Lvalue.Push(lval);
			Rvalue.Push(rval);
		}
		else return false;
	}
	if (Lvalue.Size()==Rvalue.Size()) return true;
	else                              return false;
}

inline void LineParser::PathSubstituteEnv()
{
	/* This function substitute any environmental variable from a path string
	 * Example:
	 *  input:    /one/two/$THREE/four
	 *  output:   /one/two/three/four    if $THREE=="three"
	 */

	// Find if the string start with a slash "/" => fullpath
	String path;
	bool is_fullpath = ( this->str().c_str()[0]=='/' ? true : false );
	if (is_fullpath) path.append("/");

	// Substitute environmental variables
	Array<String> pieces;
	this->SplitLine("/",pieces);


	for (size_t k=0; k<pieces.Size(); ++k)
	{
		char const * piece = pieces[k].CStr();
		if (piece[0]=='$')
		{
			String envvar_name(pieces[k].substr(1).c_str()); // substr(1) to remove the "$"
			char * envvar_res = getenv(envvar_name.CStr());
			if (envvar_res==NULL) throw new Fatal(_("Could not find \"%s\" enviroment variable in this system."),envvar_name.CStr());
			String envvar_value(envvar_res);
			//std::cout << "envvar_name  = " << envvar_name << std::endl;
			//std::cout << "envvar_value = " << envvar_value << std::endl;
			if (k==0)
			{
				char const * tmp_piece = envvar_value.CStr();
				if (tmp_piece[0]!='/' && is_fullpath) path.append("/");
			}
			else path.append("/");
			path.append(envvar_value);
		}
		else
		{
			if (k==0)
			{
				if (piece[0]!='/' && is_fullpath) path.append("/");
			}
			else path.append("/");
			path.append(pieces[k]);
		}
	}

	this->Reset(path);

}

inline void LineParser::FileBasename(String const & ExtensionToRemove, String & Basename)
{
	/* Returns the basename of a file, ex.:
	 *  In:   /home/dorival/test.dat
	 *  Out:  test.dat   if ExtensionToRemove==""     or
	 *        test       if ExtensionToRemove=="dat"
	 */

	Array<String> pieces;
	this->SplitLine("/",pieces);
	Basename = pieces[pieces.Size()-1];
	if (ExtensionToRemove!=String(""))
	{
		if (ExtensionToRemove.size()!=3)
			throw new Fatal(_("LineParser::FileBasename: The ExtensionToRemove (=%s) must have three characters, ex.: \"dat\""),ExtensionToRemove.CStr());
		String extension(Basename.substr(Basename.size()-3).c_str());
		if (extension==ExtensionToRemove)
			Basename = Basename.substr(0,Basename.size()-4).c_str();
	}

}

#endif //MECHSYS_FILEPARSER
