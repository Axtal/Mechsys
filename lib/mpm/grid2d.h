/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso                                *
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

/* Grid2D - Copyright (C) 2007 Dorival de Moraes Pedroso */

#ifndef MPM_GRID2D_H
#define MPM_GRID2D_H

// Local
#include <mechsys/mpm/defs.h>

namespace MPM {

class Grid2D
{
public:
	/* Constructor. */
	Grid2D (double xMin, double yMin, size_t nRow, size_t nCol, double Lx, double Ly, ptIsFixed pIsFixed);

	// Methods
	size_t nNodes        () const { return _X.Size(); }               ///< Return the number of nodes = nRow*nCol
	size_t nCells        () const { return _c.Size(); }               ///< Return the number of cells = (nRow-1)*(nCol-1)
	bool   IsInside      (Vector3D const & P, bool WithBorder=false); ///< Check if a given point (P) is inside the grid
	void   ClearState    (bool CurrentCoords=true);                   ///< Reset all state values in the grid, such as mass, velocity, momentum, etc.
	void   SetFixed      (size_t i, FixType Type);                    ///< Set fixed node # i
	void   SetLoading    (size_t i, Vector3D const & Fe);             ///< Set external force (loading) on node # i
	void   FixNodes      (Array<Vector3D> * VelOrDqDt);               ///< Fix nodes by clearing v (velocities) or dqdt (rate of momentum)
	void   ApplyLoadings (double M);                                  ///< Apply loading using data on loaded nodes and scaling by the multiplier M
	void   Refine        ();                                          ///< Refine the grid (bi-section)
	size_t LBN           (Vector3D const & P, bool MPM=false);        ///< Left-bottom node
	bool   IsFixed       (size_t n, size_t iD) const;                 ///< Check if node n is fixed along iD dimension

	// Access methods
	double                  xMin ()             const { return _xmin;          } ///< Returns Min x
	double                  yMin ()             const { return _ymin;          } ///< Returns Min y
	double                  xMax ()             const { return _xmax;          } ///< Returns Max x
	double                  yMax ()             const { return _ymax;          } ///< Returns Max y
	size_t                  nRow ()             const { return _nrow;          } ///< Returns Num of rows
	size_t                  nCol ()             const { return _ncol;          } ///< Returns Num of columns
	double                  L    (size_t iComp) const { return _L(iComp);      } ///< Returns Cell-length 0=>x, 1=>y
	Vector3D const &        L    ()             const { return _L;             } ///< Returns Cell-length 0=>x, 1=>y
	Array<Vector3D> const & X    ()             const { return _X;             } ///< Access all nodes (Undeformed coordinates)
	Array<Vector3D> const & x    ()             const { return _x;             } ///< Access all nodes (Current coordinates)
	Vector3D const &        X    (size_t i)     const { return _X[i];          } ///< Returns access to a node position - Undeformed (read)
	Vector3D const &        x    (size_t i)     const { return _x[i];          } ///< Returns access to a node position - Current (read)
	Vector3D       &        x    (size_t i)           { return _x[i];          } ///< Returns access to a node position - Current (write)
	Connec2D const &        C    (size_t i)     const { return _c[i];          } ///< Returns access to a connectivity (read)
	size_t                  nFix ()             const { return _fixed.Size();  } ///< Returns the number of fixed nodes
	Vector3D const &        FixN (size_t i)     const { return _X[_fixed[i]];  } ///< Access a fixed node (Undeformed coordinates)
	FixType                 FixT (size_t i)     const { return _fixtype[i];    } ///< Access the type of a fixed node
	size_t                  nLd  ()             const { return _loaded.Size(); } ///< Returns the number of loaded nodes
	Vector3D const &        LdN  (size_t i)     const { return _X[_loaded[i]]; } ///< Access a loaded node (Undeformed coordinates)
	Vector3D const &        Ld   (size_t i)     const { return _loads[i];      } ///< Access a loading

	// Set methods
	double   & m    (size_t i) { return _m   [i]; } ///< Nodes masses
	Vector3D & q    (size_t i) { return _q   [i]; } ///< Nodes momenta (x,y comps)
	Vector3D & v    (size_t i) { return _v   [i]; } ///< Nodes velocities (x,y comps)
	Vector3D & fe   (size_t i) { return _fe  [i]; } ///< Nodes external forces (x,y comps)
	Vector3D & fi   (size_t i) { return _fi  [i]; } ///< Nodes internal forces (x,y comps)
	Vector3D & dqdt (size_t i) { return _dqdt[i]; } ///< Rate of nodes momenta (x,y comps)
	Array<Vector3D> & q    () { return _q   ; } ///< Nodes momenta (x,y comps)
	Array<Vector3D> & v    () { return _v   ; } ///< Nodes velocities (x,y comps)
	Array<Vector3D> & dqdt () { return _dqdt; } ///< Rate of nodes momenta (x,y comps)

	// Get methods
	double   const & m    (size_t i) const { return _m   [i]; } ///< Nodes masses
	Vector3D const & q    (size_t i) const { return _q   [i]; } ///< Nodes momenta (x,y comps)
	Vector3D const & v    (size_t i) const { return _v   [i]; } ///< Nodes velocities (x,y comps)
	Vector3D const & fe   (size_t i) const { return _fe  [i]; } ///< Nodes external forces (x,y comps)
	Vector3D const & fi   (size_t i) const { return _fi  [i]; } ///< Nodes internal forces (x,y comps)
	Vector3D const & dqdt (size_t i) const { return _dqdt[i]; } ///< Rate of nodes momenta (x,y comps)

private:
	// Input
	double     _xmin;    ///< Min x
	double     _ymin;    ///< Min y
	size_t     _nrow;    ///< Num of rows
	size_t     _ncol;    ///< Num of columns
	Vector3D  _L;        ///< Cell-length 0=>x, 1=>y
	ptIsFixed _is_fixed; ///< Point to function which detects if a node is fixed

	// Geometry
	double          _xmax; ///< Max x
	double          _ymax; ///< Max y
	Array<Vector3D> _X;    ///< Undeformed (x,y coords) (size=nRow*nCol)
	Array<Vector3D> _x;    ///< Current (x,y coords) (size=nRow*nCol)
	Array<Connec2D> _c;    ///< Connectivities (size=(nRow-1)*(nCol-1))

	// State values
	Array<double>   _m;    ///< Nodes masses
	Array<Vector3D> _q;    ///< Nodes momenta (x,y comps)
	Array<Vector3D> _v;    ///< Nodes velocities (x,y comps)
	Array<Vector3D> _fe;   ///< Nodes external forces (x,y comps)
	Array<Vector3D> _fi;   ///< Nodes internal forces (x,y comps)
	Array<Vector3D> _dqdt; ///< Rate of nodes momenta (x,y comps)

	// Fixed nodes
	Array<size_t>  _fixed;   ///< Indexes of the fixed nodes
	Array<FixType> _fixtype; ///< Type of the fixing

	// Loaded nodes
	Array<size_t>   _loaded; ///< Indexes of the loaded nodes
	Array<Vector3D> _loads;  ///< Loading values (applied/external forces)

	// Private methods
	void _refine(); ///< Initialize/refine the grid

}; // class Grid2D


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline Grid2D::Grid2D(double xMin, double yMin, size_t nRow, size_t nCol, double Lx, double Ly, ptIsFixed pIsFixed)
	: _xmin     (xMin),
	  _ymin     (yMin),
	  _nrow     (nRow),
	  _ncol     (nCol),
	  _is_fixed (pIsFixed)
{
	// Save lengths
	_L = Lx, Ly, 0.0;

	// Do generate grid
	_refine ();
}

inline bool Grid2D::IsInside(Vector3D const & P, bool WithBorder)
{
	if (WithBorder)
	{
		if (P(0)>=_xmin+_L(0) && P(0)<=_xmax-_L(0) && P(1)>=_ymin+_L(1) && P(1)<=_ymax-_L(1)) return true;
		else                                                                                  return false;
	}
	else
	{
		if (P(0)>=_xmin && P(0)<=_xmax && P(1)>=_ymin && P(1)<=_ymax) return true;
		else                                                          return false;
	}		
}

inline void Grid2D::ClearState(bool CurrentCoords)
{
	for (size_t i=0; i<_X.Size(); ++i)
	{
		if (CurrentCoords)
		_x   [i] = _X[i]; // nodes undeformed position
		_m   [i] = 0.0;   // nodes masses
		_q   [i] = 0.0;   // nodes momenta
		_v   [i] = 0.0;   // nodes velocities
		_fe  [i] = 0.0;   // nodes external forces
		_fi  [i] = 0.0;   // nodes internal forces
		_dqdt[i] = 0.0;   // nodes momenta rate
	}
}

inline void Grid2D::SetFixed(size_t i, FixType Type)
{
	long k = _fixed.Find(i);
	if (k<0) // not added yet
	{
		_fixed  .Push(i);
		_fixtype.Push(Type);
	}
	else
	{
		_fixed  [k] = i;
		_fixtype[k] = Type;
	}
}

inline void Grid2D::SetLoading(size_t i, Vector3D const & Fe)
{
	long k = _loaded.Find(i);
	if (k<0) // not added yet
	{
		_loaded.Push(i);
		_loads .Push(Fe);
	}
	else
	{
		_loaded[k] = i;
		_loads [k] = Fe;
	}
}

inline void Grid2D::FixNodes(Array<Vector3D> * VelOrDqDt)
{
	for (size_t i=0; i<_fixed.Size(); ++i)
	{
		switch (_fixtype[i])
		{
			case FIX_X   : { (*VelOrDqDt)[_fixed[i]](0)=0.0;                                                                 break; }
			case FIX_Y   : { (*VelOrDqDt)[_fixed[i]](1)=0.0;                                                                 break; }
			case FIX_Z   : { (*VelOrDqDt)[_fixed[i]](2)=0.0;                                                                 break; }
			case FIX_XY  : { (*VelOrDqDt)[_fixed[i]](0)=0.0; (*VelOrDqDt)[_fixed[i]](1)=0.0;                                 break; }
			case FIX_YZ  : { (*VelOrDqDt)[_fixed[i]](1)=0.0; (*VelOrDqDt)[_fixed[i]](2)=0.0;                                 break; }
			case FIX_ZX  : { (*VelOrDqDt)[_fixed[i]](2)=0.0; (*VelOrDqDt)[_fixed[i]](0)=0.0;                                 break; }
			case FIX_XYZ : { (*VelOrDqDt)[_fixed[i]](0)=0.0; (*VelOrDqDt)[_fixed[i]](1)=0.0; (*VelOrDqDt)[_fixed[i]](2)=0.0; break; }
		}
	}
}

inline bool Grid2D::IsFixed(size_t n, size_t iD) const
{
	long k = _fixed.Find(n);
	if (k<0) return false; // not fixed
	else
	{
		switch (_fixtype[k])
		{
			case FIX_X   : { if (iD==0)          return true; else return false; }
			case FIX_Y   : { if (iD==1)          return true; else return false; }
			case FIX_Z   : { if (iD==2)          return true; else return false; }
			case FIX_XY  : { if (iD==0 || iD==1) return true; else return false; }
			case FIX_YZ  : { if (iD==1 || iD==2) return true; else return false; }
			case FIX_ZX  : { if (iD==2 || iD==0) return true; else return false; }
			case FIX_XYZ : {                     return true;                    }
		}
	}
	return false;
}

inline void Grid2D::ApplyLoadings(double M)
{
	for (size_t i=0; i<_loaded.Size(); ++i)
		_fe[_loaded[i]] = M*_loads[i];
}

inline void Grid2D::Refine()
{
	// New geometry
	_nrow = 2*_nrow-1;
	_ncol = 2*_ncol-1;
	_L    = _L(0)/2.0, _L(1)/2.0, 0.0;

	// Do refine
	_refine();
}

inline size_t Grid2D::LBN(Vector3D const & P, bool MPM)
{
	if (MPM) return static_cast<size_t>((P(0)-_xmin)/_L(0))   + (static_cast<size_t>((P(1)-_ymin)/_L(1))  )*_ncol;
	else     return static_cast<size_t>((P(0)-_xmin)/_L(0))-1 + (static_cast<size_t>((P(1)-_ymin)/_L(1))-1)*_ncol;
}

/* private */

inline void Grid2D::_refine()
{
	// Allocate nodes
	_X.Resize(_nrow*_ncol);

	// Set nodes coordinates
	size_t k = 0;
	for (size_t i=0; i<_nrow; ++i)
	for (size_t j=0; j<_ncol; ++j)
	{
		_X[k] = _xmin+j*_L(0),
		        _ymin+i*_L(1),
				0.0;
		k++;
	}
	_xmax = _xmin+(_ncol-1)*_L(0);
	_ymax = _ymin+(_nrow-1)*_L(1);

	// Allocate and initialize (zero) nodes data
	_x   .Resize(_X.Size()); // nodes current coordinates
	_m   .Resize(_X.Size()); // nodes masses
	_q   .Resize(_X.Size()); // nodes momenta
	_v   .Resize(_X.Size()); // nodes velocities
	_fe  .Resize(_X.Size()); // nodes external forces
	_fi  .Resize(_X.Size()); // nodes internal forces
	_dqdt.Resize(_X.Size()); // nodes momenta rate
	ClearState();

	// Allocate cells connectivities
	_c.Resize((_nrow-1)*(_ncol-1));

	// Set connectivities
	for (size_t i=0; i<_nrow-1; ++i)
	{
		for (size_t j=0; j<_ncol-1; ++j)
		{
			_c[i*(_ncol-1)+j] = i   *_ncol+j  ,
			                    i   *_ncol+j+1,
			                   (i+1)*_ncol+j+1,
			                   (i+1)*_ncol+j  ;
		}
	}

	// Set current coordinates and boundary conditions
	_fixed  .Resize(0);
	_fixtype.Resize(0);
	FixType type;
	for (size_t i=0; i<_X.Size(); ++i)
	{
		_x[i] = _X[i];
		if ((*_is_fixed) (_X[i], type))
		{
			long k = _fixed.Find(i);
			if (k<0) // not added yet
			{
				_fixed  .Push(i);
				_fixtype.Push(type);
			}
			else
			{
				_fixed  [k] = i;
				_fixtype[k] = type;
			}
		}
	}
}

}; // namespace MPM

#endif // MPM_GRID2D_H
