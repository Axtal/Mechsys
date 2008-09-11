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
	void         SetPrms        (char const * Prms);
	void         SetInis        (char const * Inis);
	void         TgConductivity (Matrix<double> & D) const { D = _K; };
	int          UpdateState    (Vector<double> const & DEps, Vector<double> & DSig);
	void         BackupState    ();
	void         RestoreState   ();
	char const * Name           () const { return "LinHeat"; }

private:
	// Data
	Matrix<double> _K; // Conductivity

	// Private methods
	double _val (char const * Name) const { throw new Fatal("LinHeat::_val The Name==%s is invalid",Name); }

}; // class LinHeat


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void LinHeat::SetPrms(char const * Prms)
{
	if (_geom<0) throw new Fatal("LinHeat::SetPrms: Geometry type:\n\t[1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)] must be set via SetGeom before calling this method");

	/* "kxx=1.0 kxy=0.0 kxz=0.0
	            kyy=1.0 kxz=0.0
	                    kzz=1.0"
	    or "k=1.0" => isotropic
	*/
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Conductivity matrix
	     if (_geom==1)                         _k.Resize (1,1);
	else if (_geom==2 || _geom==4 || _geom==5) _k.Resize (2,2);
	else if (_geom==3)                         _k.Resize (3,3);
	_k.SetValues (0.0);

	// Set
	if (names.Size()==1)
	{
		if (names[0]=="k")
		{
			if (_geom==1)
			{
				_k(0,0) = values[0];
			}
			else if (_geom==2 || _geom==4 || _geom==5)
			{
				_k = values[0],       0.0,
						   0.0, values[0];
			}
			else if (_geom==3)
			{
				_k = values[0],       0.0,       0.0,
						   0.0, values[0],       0.0,
						   0.0,       0.0, values[0];
			}
		}
		else throw new Fatal("LinHeat::SetPrms: Parameter key==%s for isotropic models is invalid. It must be equal to 'k'. Ex.: k=1.0",names[0].CStr());
	}
	else
	{
		if (_geom==1) throw new Fatal("LinHeat::SetPrms: For unidimensional problems, only one parameter key (equal to 'k') must be given. Ex.: k=1.0 (%s is invalid)");
		for (size_t i=0; i<names.Size(); ++i)
		{
			      if (names[i]=="kxx"                    && _geom> 1) _k(0,0) = values[i];
			 else if (names[i]=="kxy" || names[i]=="kyx" && _geom> 1) _k(0,1) = values[i];
			 else if (names[i]=="kxz" || names[i]=="kzx" && _geom==3) _k(0,2) = values[i];
			 else if (names[i]=="kyy"                    && _geom> 1) _k(1,1) = values[i];
			 else if (names[i]=="kyz" || names[i]=="kzy" && _geom==3) _k(1,2) = values[i];
			 else if (names[i]=="kzz"                    && _geom==3) _k(2,2) = values[i];
			 else throw new Fatal("LinHeat::SetPrms: Parameter key==%s is invalid. It must be: kxx, kxy, kxz,  kyy, kyz,  kzz  (or kyx, kzx, kzy).",names[i].CStr());
		}
		_k(1,0) = _k(0,1);
		_k(2,0) = _k(0,2);
		_k(2,1) = _k(1,2);
	}
}

inline void LinHeat::SetInis(char const * Inis)
{
	/* "T=0.0" */
	LineParser lp(Inis);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Check
	for (size_t i=0; i<names.Size(); i++)
	{
		if (names[i]=="T") _T = values[i];
		else throw new Fatal("LinHeat::SetInis: Initial value key==%s is invalid. It must be 'T'. Ex.: T=1.0",names[i]);
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
