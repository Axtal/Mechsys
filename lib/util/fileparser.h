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

#ifndef MECHSYS_FILEPARSER_H
#define MECHSYS_FILEPARSER_H

// STL
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>

// MechSys
#include "util/array.h"
#include "util/string.h"
#include "util/lineparser.h"
#include "util/fatal.h"

class FileParser
{
public:
	// Constructor
	FileParser(String const & Filename) : _comm_char('#'), _filename(Filename) { _open(_filename); }
	// Destructor
	~FileParser() { _file.close(); }
	// Structures
	struct Tags
	{
		String L1;
		String L2;
		String Comment;
		bool   Required;
		String Default;
	};
	typedef std::map< String,Array<double> > Table;
	typedef std::map< String,String >        StructVals;
	// Methods
	int         CurrentLineNumber    ();
	String      GetCurrentLine       ();
	bool        TryGetCurrentLine    (String & line);
	void        Advance              ();
	void        JumpCommentsOrBlanks ();
	template<typename Type>
	void        ReadFirstWord        (Type & Val);
	bool        IsEOF                () { return _is_EOF; }
	void        Reset                () { _file.close(); _open(_filename); }
	template<typename Type>
	void        ReadNextValue        (Type & Val) { JumpCommentsOrBlanks(); ReadFirstWord(Val); Advance(); }
	template<typename Type>
	void        FindKeyAndFillArray  (String const & StrKey, Array<Type> & A);
	void        SetCommChar          (char CommChar) { _comm_char = CommChar; }
	void        StructFile           (Tags const * tags, int num_tags, StructVals & values);
	void        ReadTable            (Table & t);
	// Static
	static bool CheckForFile(String const & Filename)
	{
		// returns true if file exists and false otherwise
		std::ifstream file(Filename.CStr(),std::ios::in);
		if     (file.fail()) { return false; }
		else { file.close();   return true;  }
	}
	static String GetBackupFilename(String const & originalFN)
	{
		// ex.:
		//         originalFN                    |   newFN
		//        -------------------------------------------
		//        name.txt                       |  name.txt~
		//        name.txt~                      |  name.txt~
		//        name.txt<file does not exist>  |  name.txt
		//
		// OBS.: return originalFN if file does not exist
		
		// Check if file exists
		if (FileParser::CheckForFile(originalFN)) // file exists
		{
			if ( originalFN[ originalFN.size()-1 ] == '~' ) // already has a tilde at the end
				return originalFN;
			else
			{
				String tmp_fn = originalFN;
				tmp_fn.append(_("~"));
				return tmp_fn;
			}
		}
		else
			return originalFN;
	}
private:
	// Data
	char          _comm_char; // Comments char, ex.: '#'
	bool          _is_EOF;
	String        _filename;
	std::ifstream _file;
	String        _str_current_line;
	int           _current_line_num;
	// Functions
	void _open(String const & Filename);
}; // 


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline int FileParser::CurrentLineNumber()
{
	if (_is_EOF)
		throw new Fatal(_("FileParser::GetCurrentLine: Can not return current line number because EOF. Filename = < %s >"), _filename.c_str());
	return _current_line_num; 
}

inline String FileParser::GetCurrentLine()
{
	if (_is_EOF)
		throw new Fatal(_("FileParser::GetCurrentLine: Can not return current line number because EOF. Filename = < %s >"), _filename.c_str());
	return _str_current_line; 
}

inline bool FileParser::TryGetCurrentLine(String & line)
{
	if (_is_EOF) return false;
	else
	{
		line = _str_current_line;
		return true;
	}
}

inline void FileParser::Advance()
{
	// Check EOF
	if (_is_EOF) 
		throw new Fatal(_("FileParser::Advance: < %s > Can not advance because EOF"),_filename.c_str());

	// If not is EOF read line
	std::string buf;
	if (!std::getline(_file, buf)) _is_EOF = true;
	else
	{
		_str_current_line = buf;
		_current_line_num++;
	}
}

inline void FileParser::JumpCommentsOrBlanks()
{
	do
	{
		std::istringstream iss(_str_current_line.CStr());
		String word;
		if (iss>>word) // Comments or data
		{
			if (word[0]==_comm_char)
			{
				if (!TryGetCurrentLine(_str_current_line)) break;
				Advance();
			}
			else
				break;
		}
		else // Blank lines (?)
		{
			if (!TryGetCurrentLine(_str_current_line)) break;
			Advance();
		}
	} while (true);
}

template<typename Type>
inline void FileParser::ReadFirstWord(Type & Val)
{
	std::istringstream iss(_str_current_line.CStr());
	if (!(iss >> Val))
		throw new Fatal(_("FileParser::ReadFirstWord: Invalid data value. _current_line_num = %d"), _current_line_num);
}

inline void FileParser::_open(String const & Filename)
{
	// EOF
	_is_EOF = false;

	// Open file
	_file.open(Filename.CStr(), std::ios::in);
	if (_file.fail())
		throw new Warning(_("FileParser::Constructor: Could not open file < %s >"), Filename.c_str());
	
	// Read current (first) line
	std::string buf;
	if (!std::getline(_file, buf))
		throw new Warning(_("FileParser::Constructor: Could not read first line of file < %s >"), Filename.c_str());
	_str_current_line = buf;
	_current_line_num = 1;
}

template<typename Type>
inline void FileParser::FindKeyAndFillArray(String const & StrKey, Array<Type> & A)
{
	/* filename.example
	################################################ comments

	WrongKey 1.1 2.2 3.3 4.4 5.5 6.6
	
	#-------------------------------------------

	StrKey 1.1 2.2 3.3 4.4 5.5 6.6

	4 123 1234 

	*/

	// Find StrKey
	bool is_key_found = false;
	A.Resize(0);
	do
	{
		// Check EOF
		if (_is_EOF) 
			throw new Fatal(_("FileParser::FindKeyAndFillArray: Could not find StrKey = < %s > inside file < %s >"), StrKey.c_str(), _filename.c_str());

		// Find str_key == ModelName
		JumpCommentsOrBlanks();
		LineParser LP(GetCurrentLine());
		String str_key;
		LP >> str_key;
		if (str_key==StrKey)
		{
			Type value;
			while (LP>>value)
			{
				A.Push(value);
			}
			is_key_found = true;
		}
		else
		{
			try { Advance(); }
			catch (Exception * e)
			{
				delete e;
				throw new Fatal(_("FileParser::FindKeyAndFillArray: Could not find StrKey = < %s > inside file < %s >"), StrKey.c_str(), _filename.c_str());
			}
		}
	} while (!is_key_found);
}

inline void FileParser::StructFile(Tags const * tags, int num_tags, StructVals & values)
{
	// Check if file wasn't read
	if (IsEOF()) throw new Fatal(_("FileParser::StructFile: Could not read file < %s >. Probably file was already read."), _filename.c_str());

	// Initialize values
	values.clear();
	for (int i=0; i<num_tags; ++i)
	{
		String fullkey(tags[i].L1); fullkey.append(_(".")); fullkey.append(tags[i].L2);
		values[fullkey] = "";
	}
	
	// Read file
	while (!IsEOF())
	{
		JumpCommentsOrBlanks();
		String line;
		if (!TryGetCurrentLine(line)) break;
		LineParser LP(line);
		Advance();
		String key;    LP>>key;
		String subkey; LP>>subkey;
		String fullkey(key); fullkey.append(_(".")); fullkey.append(subkey);
		bool found = false;
		for (int i=0; i<num_tags; ++i)
		{
			if (key==tags[i].L1 && subkey==tags[i].L2)
			{
				String val;
				while (LP>>val)
				{
					if (val[0]==_comm_char) break;
					else
					{
						if (values[fullkey].size()>0) values[fullkey].append(_(" "));
						values[fullkey].append(val);
					}
				}
				found = true;
				break;
			}
		}
		if (!found)
			throw new Fatal(_("FileParser::StructFile: File < %s >: Pair < %s %s > is not valid."), _filename.c_str(), key.c_str(), subkey.c_str());
	}

	// Additional check and configuration
	for (int i=0; i<num_tags; ++i)
	{
		String fullkey(tags[i].L1); fullkey.append(_(".")); fullkey.append(tags[i].L2);

		// Check required values
		if (tags[i].Required && values[fullkey]=="")
			throw new Fatal(_("FileParser::StructFile: The pair < %s %s > must be provided by file < %s >"), tags[i].L1.c_str(), tags[i].L2.c_str(), _filename.c_str());

		// Set default values
		if (!tags[i].Required && tags[i].Default!=_("") && values[fullkey]==_(""))
			values[fullkey] = tags[i].Default;
	}
	
}

inline void FileParser::ReadTable(Table & T)
{
	// Check if file wasn't read
	if (IsEOF()) throw new Fatal(_("FileParser::ReadTable: Could not read file < %s >. Probably file was already read."), _filename.c_str());

	// Read Header ==> num of columns
	Array<String> header;
	LineParser LP(GetCurrentLine());
	String key;
	while (LP>>key) header.Push(key);
	if (header.Size()<1) throw new Fatal(_("FileParser::ReadTable: File format invalid < %s >. Header must have at least one tag (one column)"), _filename.c_str());

	// Prepare table (map)
	T.clear(); // clear current values
	Array<double> a_null;
	for (size_t i=0; i<header.Size(); ++i)
		T[header[i]] = a_null;

	// Read values
	size_t num_lines=0;
	while (!IsEOF())
	{
		// read a line and distribute the values to all columns
		Advance();
		String line;
		if (!TryGetCurrentLine(line)) break;
		LP.Reset(line);
		size_t i_col=0; double col_val;
		while (LP>>col_val)
		{
			if (i_col>header.Size()-1)
				throw new Fatal(_("FileParser::ReadTable: File < %s >: Too many values: Number of columns (%d) must be equal to the number (%d) of header tags"), _filename.c_str(), i_col+1, header.Size());
			T[header[i_col]].Push(col_val);
			i_col++;
		}
		if (i_col!=header.Size())
			throw new Fatal(_("FileParser::ReadTable: File < %s >: Not enough values: Number of columns (%d) must be equal to the number (%d) of header tags"), _filename.c_str(), i_col, header.Size());
		num_lines++;
	}

	// Check-up
	Table::const_iterator it = T.begin();
	while (it!=T.end())
	{
		if (it->second.Size()!=num_lines)
			throw new Fatal(_("FileParser::ReadTable:  INTERNAL ERROR:  There is an inconsistency with the number of lines of each column (%d != %d)"), it->second.Size(), num_lines);
		it++;
	}
	
}


std::ostream & operator<< (std::ostream & os, FileParser::Table const & T)
{
	FileParser::Table::const_iterator it = T.begin();
	while (it!=T.end())
	{
		os << "Column " << it->first << " has " << it->second.Size() << " values:\n";
		for (size_t i=0; i<it->second.Size(); ++i)
			os << it->second[i] << " ";
		os << "\n\n";
		it++;
	}
	return os;
}


#endif //MECHSYS_FILEPARSER
