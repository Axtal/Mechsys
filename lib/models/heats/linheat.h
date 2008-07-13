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

#ifndef MECHSYS_LINHEAT_H
#define MECHSYS_LINHEAT_H

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "models/heatmodel.h"
#include "util/string.h"
#include "util/util.h"
#include "util/lineparser.h"

using LinAlg::Vector;
using LinAlg::Matrix;

class LinHeat : public HeatModel
{
public:
	// Destructor
	virtual ~LinHeat () {}

	// Derived Methods
	void         SetPrms      (char const * Prms);
	void         SetInis      (char const * Inis);
	void         TgStiffness  (Matrix<double> & D) const { D = _K; };
	int          UpdateState  (Vector<double> const & DEps, Vector<double> & DSig);
	void         BackupState  ();
	void         RestoreState ();
	char const * Name         () const { return "LinHeat"; }

private:
	// Data
	Matrix<double> _K;

	// Private methods
	double _val (char const * Name) const { throw new Fatal("LinHeat::_val The Name==%s is invalid",Name); }

}; // class LinHeat


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void LinHeat::SetPrms(char const * Prms)
{
	//if (_geom<0) throw new Fatal("LinHeat::SetPrms: Geometry type:\n\t[1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)] must be set via SetGeom before calling this method");

	/* "E=20000.0 nu=0.2" */
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Set
	double k  = 0.0;
	for (size_t i=0; i<names.Size(); ++i)
	{
			 if (names[i]=="k" ) k  = values[i];
	}

	//_K.Resize(3,3);
	//_K =   k, 0.0, 0.0,
	//     0.0,   k, 0.0,
	//     0.0, 0.0,   k;
	_K.Resize(2,2);
	_K =   k, 0.0,
	     0.0,   k;
}

inline void LinHeat::SetInis(char const * Inis)
{
	/* "Sx=0.0 Sy=0.0 Sxy=0.0" */
	LineParser lp(Inis);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Check
	for (size_t i=0; i<names.Size(); i++)
	{
	}
}

inline int LinHeat::UpdateState(Vector<double> const & DGrad, Vector<double> & DFlow)
{
	DFlow = _K*DGrad;
	return 1;
}

inline void LinHeat::BackupState()
{
}

inline void LinHeat::RestoreState()
{
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new LinHeat model
Model * LinHeatMaker()
{
	return new LinHeat();
}

// Register an LinHeat model into ModelFactory array map
int LinearHeatRegister()
{
	ModelFactory["LinHeat"] = LinHeatMaker;
	return 0;
}

// Execute the autoregistration
int __LinearHeat_dummy_int = LinearHeatRegister();


#endif // MECHSYS_LINHEAT_H
