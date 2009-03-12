/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo, Fernando Alonso                   *
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

#ifndef MECHSYS_FILEUTILS_H
#define MECHSYS_FILEUTILS_H

// External
#include <sys/stat.h> // for mkdir()

// MechSys
#include "util/string.h"
#include "util/lineparser.h"

namespace Util
{

bool MkdirP(String const & Path)
{
	// tree
	LineParser    lp(Path);
	Array<String> tree;
	lp.SeparatedLine("/",tree);

	int res = -1;
	if (tree.Size()>0)
	{
		String left = _T("");
		for (size_t i=0; i<tree.Size(); ++i)
		{
			left += tree[i] + _T("/");
			// Upon successful completion, mkdir() shall return 0.
			// Otherwise, -1 shall be returned, no directory shall be created, and errno shall be set to indicate the error.
			res = mkdir(left.GetSTL().c_str(),S_IRWXU); //if (res!=) // directory may be existent, so, it is ok...
		}
	}

	return res==0 ? true : false;
}

}; // namespace Util

#endif // MECHSYS_FILEUTILS_H
