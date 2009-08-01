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

#ifndef DEM_DOMAIN2D_H
#define DEM_DOMAIN2D_H

// MechSys
#include "dem/disk.h"
#include "util/array.h"

class Domain2D
{
public:
	// Constructor
	Domain2D (double rho=1.0, double dt=1.0, size_t np=100);

	// Methods
	void Forces (double k=1.0);
	void Solve  (double tini, double tfin, double dtout);

protected:
	Array<Disk> _P;
	double      _dt;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Domain2D::Domain2D(double rho, double dt, size_t np)
{
	_P.Resize(np);
	_dt = dt;
	for (size_t i=0; i<10; i++)
	{
		for (size_t j=0; j<10; j++)
		{
			Vec2_t x0;
			x0 = i,j;
			double vx = (rand()/RAND_MAX),vy = (rand()/RAND_MAX);
			Vec2_t v0; 
			v0 = vx,vy;
			_P[j+i*10] = Disk(rho,0.5,x0,v0,1.0);
		}
	}
}

inline void Domain2D::Forces(double k)
{
	for (size_t i=0;i<100;i++)
	{
		for (size_t j=i+1;j<100;j++)
		{
			Vec2_t dx = _P[i].x() - _P[j].x();
			double dist = sqrt(dx(1)*dx(1)+dx(0)*dx(0));
			double delta = _P[i].r() + _P[j].r() - dist;
			if (delta>0) {
				Vec2_t dF = (k*delta/dist)*dx;
				_P[i].F() = _P[i].F() + dF;
				_P[j].F() = _P[j].F() - dF;
			}
		}
	}
}

inline void Domain2D::Solve(double tini,double tfin,double dtout)
{
	double t=tini;
	double tout=t+dtout;
	while (t<tfin)
	{
		for (size_t i=0; i<100; i++) _P[i].F() = 0,0;
		Forces();
		for (size_t i=0; i<100; i++) _P[i].Move();
		t += _dt;
	}
}

#endif // DEM_DOMAIN2D_H